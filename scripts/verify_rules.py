#!/usr/bin/env python3
"""Rule-response checks from the Handover (Batches 1-5).

For each rule: run baseline / forced-off / forced-on over several seeds and
report the observed column at days 3, 6, 9. A rule passes if 'off' and 'on'
both move the column away from baseline in the expected direction.
"""
import argparse, csv, json, os, re, statistics as st, subprocess, sys
from concurrent.futures import ProcessPoolExecutor

DAYS = {3: 144, 6: 288, 9: 432}
B = "biology.fibroblast."
M = "biology.biomaterial."

# Cytokine base diffusivities in gh*.json (TNF, TGF, FGF, IL6, IL8, IL10), mm^2/min.
# Q no longer scales D in the code; the D tests scale all six together (Nicole, Nov 11).
D_BASE = [0.00018, 0.000156, 0.000177, 0.000168, 0.000228, 0.000139]
Q_REF = 15.71   # GH5 Q at 1400 min; config D x Q_REF = the D the code used before Q was removed
def d_scale(f):
    return {f"chemistry.species.{i}.base_diffusivity_mm2_per_min": f * d
            for i, d in enumerate(D_BASE)}



# name: (column prefix, off overrides, on overrides, expected off, expected on)
TESTS = {
    "proliferation": ("Total Cells",
        {B+"proliferation.baseline": 0, B+"proliferation.activating_factor_effect": 0, B+"proliferation.baseline_offset": 0},
        {B+"proliferation.baseline": 250}, "lower", "higher"),
    "migration": ("Mean Displacement",
        {B+"migration.baseline_speed": 0, B+"ha_synthesis.move_probability": 0},
        {B+"migration.baseline_speed": 5}, "lower", "higher"),
    "migration2x": ("Mean Displacement",
        {B+"migration.baseline_speed": 0, B+"ha_synthesis.move_probability": 0},
        {B+"migration.baseline_speed": 1.73}, "lower", "higher"),
    "viability": ("Live Cells",
        {B+"viability.time_effect": 0, B+"viability.biomaterial_effect": 100},
        {B+"viability.time_effect": 0, B+"viability.biomaterial_effect": 10}, "higher", "lower"),
    "tgf": ("Total TGF",
        {B+"tgf_synthesis.baseline": 0, B+"tgf_synthesis.feedback_effect": 0},
        {B+"tgf_synthesis.baseline": 100}, "lower", "higher"),
    "fgf": ("Total FGF", {B+"fgf_synthesis.rate": 0}, {B+"fgf_synthesis.rate": 50}, "lower", "higher"),
    "tnf": ("Total TNF", {B+"tnf_synthesis.baseline": 0}, {B+"tnf_synthesis.baseline": 100}, "lower", "higher"),
    "il6": ("Total IL6", {B+"il6_synthesis.baseline": 0}, {B+"il6_synthesis.baseline": 0.1}, "lower", "higher"),
    "il8": ("Total IL8", {B+"il8_synthesis.baseline": 0}, {B+"il8_synthesis.baseline": 10}, "lower", "higher"),
    "activation": ("Activated Fibroblast",
        {B+"activation.independent_probability": 0,
         B+"activation.tgf_intermediate_threshold": 1e12, B+"activation.tgf_high_threshold": 1e12},
        {B+"activation.independent_probability": 0.05}, "lower", "higher"),
    "deactivation": ("Activated Fibroblast",
        {B+"activation.deactivation_probability": 0},
        {B+"activation.deactivation_probability": 1.0}, "higher", "lower"),
    "collagen": ("Collagen",
        {B+"collagen_synthesis.baseline": 0, B+"collagen_synthesis.equal_frag_offset": 0,
         B+"collagen_synthesis.baseline_rate": 0},
        {B+"collagen_synthesis.baseline": 8922, B+"collagen_synthesis.baseline_rate": 100}, "lower", "higher"),
    "elastin": ("Elastin", {B+"elastin_synthesis.baseline": 0},
        {B+"elastin_synthesis.baseline": 1834}, "lower", "higher"),
    "ha": ("HA (", {B+"ha_synthesis.baseline": 0},
        {B+"ha_synthesis.baseline": 166.8}, "lower", "higher"),
    "act_multiplier": ("Collagen",
        {B+"activation.synthesis_multiplier": 0},
        {B+"activation.synthesis_multiplier": 10, B+"activation.independent_probability": 0.05},
        "lower", "higher"),
    "gate_removed": ("Collagen",
        {B+"activation.independent_probability": 0,
         B+"activation.tgf_intermediate_threshold": 1e12, B+"activation.tgf_high_threshold": 1e12},
        {B+"activation.independent_probability": 0.05},
        "lower", "higher"),
    "growth_cap": ("Total Cells",
        {B+"proliferation.carrying_capacity_per_mL": 1},      # cap ~0: no divisions
        {B+"proliferation.carrying_capacity_per_mL": 1e30},   # cap ~1: no limit
        "lower", "higher"),
    "collagen_xb10": ("Collagen",
        {B+"collagen_synthesis.baseline": 0, B+"collagen_synthesis.equal_frag_offset": 0,
         B+"collagen_synthesis.baseline_rate": 0},
        {B+"collagen_synthesis.baseline": 8922, B+"collagen_synthesis.mesh_effect": 255.4,
         B+"collagen_synthesis.elasticity_effect": 0.02429}, "lower", "higher"),
    "elastin_xb10": ("Elastin", {B+"elastin_synthesis.baseline": 0},
        {B+"elastin_synthesis.baseline": 1834, B+"elastin_synthesis.mesh_effect": 47.02,
         B+"elastin_synthesis.elasticity_effect": 0.009759}, "lower", "higher"),
        # Table 4 (biomaterial) rules: deterministic, so SD = 0 across seeds.
    "g_prime": ("Shear Storage Modulus",
        {M+"elastic_modulus.ha_crosslinker_interaction": 0},         # c4 = 0
        {M+"elastic_modulus.ha_crosslinker_interaction": 2666.4},    # c4 x2
        "lower", "higher"),
    "mesh": ("Mesh Size",
        {M+"elastic_modulus.ha_crosslinker_interaction": 2666.4},    # stiffer -> smaller xi
        {M+"elastic_modulus.ha_crosslinker_interaction": 0},         # softer  -> larger xi
        "lower", "higher"),
    "rho_xl": ("Crosslink Density",
        {M+"elastic_modulus.ha_crosslinker_interaction": 0},
        {M+"elastic_modulus.ha_crosslinker_interaction": 2666.4},
        "lower", "higher"),
    "wl": ("Dry-weight Loss",
        {M+"mass_loss.ha_time_effect": 0, M+"mass_loss.ha_effect": 0, M+"mass_loss.baseline": 0},
        {M+"mass_loss.baseline": 53},                                # c16 x10
        "lower", "higher"),

    # D test: off = D x 0.1, base = D x 1, on = D x 3. No expected direction;
    # read Total Cells and Collagen with compare_to_data.py, ignore the verdict.
    "diffusion": ("Total Cells", d_scale(0.1), d_scale(3), "lower", "higher"),
    # D test on the pre-removal D (config x Q_REF): off = x0.1, base = x1, on = x3.
    # Run with --base-dscale 15.71 so the baseline is config x Q_REF too.
    "diffusion_q": ("Total Cells", d_scale(0.1 * Q_REF), d_scale(3 * Q_REF), "lower", "higher"),

}

def load_jsonc(path):
    text = "\n".join(re.sub(r"//.*$", "", l) for l in open(path).read().splitlines())
    return json.loads(text)

def set_key(cfg, dotted, value):
    node = cfg
    *parents, leaf = dotted.split(".")
    for p in parents:
        node = node[int(p)] if isinstance(node, list) else node[p]
    if isinstance(node, list):
        node[int(leaf)] = value
        return
    if leaf not in node:
        sys.exit(f"unknown config key: {dotted}")
    node[leaf] = value


def run(job):
    a, name, case, overrides, seed = job
    out = os.path.join(a.out, f"{name}_{case}_s{seed}")
    csv_path = os.path.join(out, "Output_Biomarkers.csv")
    if not os.path.exists(csv_path):
        os.makedirs(out, exist_ok=True)
        cfg = load_jsonc(a.config)
        # wi = cfg["world_init"]
        # if wi["fibroblast_count"] > 0:   # counts in gh*.json are for the 3 mm world
        #     wi["fibroblast_count"] = max(1, round(wi["fibroblast_count"] * (float(a.world) / 3) ** 3))
        if a.base_dscale != 1.0:                       # applies to base, off and on
            for k, v in d_scale(a.base_dscale).items():
                set_key(cfg, k, v)

        for k, v in overrides.items():
            set_key(cfg, k, v)
        cfg_path = os.path.join(out, "config.json")
        json.dump(cfg, open(cfg_path, "w"), indent=1)
        env = dict(os.environ, LD_LIBRARY_PATH="build/lib:" + os.environ.get("LD_LIBRARY_PATH", ""))
        with open(os.path.join(out, "stdout.log"), "w") as log:
            subprocess.run([a.bin, "--config", cfg_path, "--seed", str(seed),
                            "--numticks", str(max(DAYS.values())),
                            "--wxw", a.world, "--wyw", a.world, "--wzw", a.world,
                            "--output-dir", out], stdout=log, stderr=subprocess.STDOUT,
                           env=env, check=True)
    return (name, case, seed, csv_path)

def read_col(csv_path, prefix):
    rows = list(csv.reader(open(csv_path)))
    head = rows[0]
    idx = next((i for i, h in enumerate(head) if h.startswith(prefix)), None)
    if idx is None:
        sys.exit(f"column starting with '{prefix}' not found in {csv_path}")
    return {d: float(rows[t + 1][idx]) for d, t in DAYS.items()}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default="configFiles/gh10.json")
    p.add_argument("--bin", default="build/bin/testRun")
    p.add_argument("--world", default="0.6")
    p.add_argument("--seeds", type=int, default=3)
    p.add_argument("--jobs", type=int, default=os.cpu_count())
    p.add_argument("--out", default="output/verify")
    p.add_argument("--only", nargs="*", help="subset of test names")
    p.add_argument("--base-dscale", type=float, default=1.0,
                help="multiply all six cytokine D in every case, including baseline")

    a = p.parse_args()

    names = a.only or list(TESTS)
    jobs = []
    for s in range(1, a.seeds + 1):
        jobs.append((a, "baseline", "base", {}, s))
        for n in names:
            jobs.append((a, n, "off", TESTS[n][1], s))
            jobs.append((a, n, "on", TESTS[n][2], s))
    with ProcessPoolExecutor(a.jobs) as ex:
        done = list(ex.map(run, jobs))
    paths = {(n, c, s): path for n, c, s, path in done}

    print(f"{'rule':14s} {'case':5s} " + " ".join(f"day{d:<12d}" for d in DAYS) + " verdict")
    for n in names:
        col, _, _, exp_off, exp_on = TESTS[n]
        stats = {}
        for case, key in (("base", "baseline"), ("off", n), ("on", n)):
            vals = [read_col(paths[(key, case, s)], col) for s in range(1, a.seeds + 1)]
            stats[case] = {d: (st.mean(v[d] for v in vals),
                               st.pstdev([v[d] for v in vals])) for d in DAYS}
        def moved(case, direction):
            m, sd = stats[case][9]; mb, sdb = stats["base"][9]
            gap = 2 * max(sd, sdb, 1e-9)
            return m < mb - gap if direction == "lower" else m > mb + gap
        ok = moved("off", exp_off) and moved("on", exp_on)
        for case in ("base", "off", "on"):
            cells = " ".join(f"{stats[case][d][0]:10.3g}±{stats[case][d][1]:<5.2g}" for d in DAYS)
            print(f"{n:14s} {case:5s} {cells} {'PASS' if ok else 'CHECK'}" if case == "on"
                  else f"{n:14s} {case:5s} {cells}")
        print()

if __name__ == "__main__":
    main()
