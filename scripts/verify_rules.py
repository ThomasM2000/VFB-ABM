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

# name: (column prefix, off overrides, on overrides, expected off, expected on)
TESTS = {
    "proliferation": ("Total Cells",
        {B+"proliferation.baseline": 0, B+"proliferation.activating_factor_effect": 0,
         B+"proliferation.baseline_offset": 0},
        {B+"proliferation.baseline": 250}, "lower", "higher"),
    "migration": ("Mean Displacement",
        {B+"migration.baseline_speed": 0, B+"ha_synthesis.move_probability": 0},
        {B+"migration.baseline_speed": 5}, "lower", "higher"),
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
        {B+"activation.independent_probability": 5}, "lower", "higher"),
    "deactivation": ("Activated Fibroblast",
        {B+"activation.deactivation_probability": 0},
        {B+"activation.deactivation_probability": 100}, "higher", "lower"),
    "collagen": ("Collagen",
        {B+"collagen_synthesis.baseline": 0, B+"collagen_synthesis.equal_frag_offset": 0,
         B+"collagen_synthesis.baseline_rate": 0},
        {B+"collagen_synthesis.baseline": 89220, B+"collagen_synthesis.baseline_rate": 100}, "lower", "higher"),
    "elastin": ("Elastin", {B+"elastin_synthesis.baseline": 0},
        {B+"elastin_synthesis.baseline": 1834}, "lower", "higher"),
    "ha": ("HA (", {B+"ha_synthesis.baseline": 0},
        {B+"ha_synthesis.baseline": 166.8}, "lower", "higher"),
}

def load_jsonc(path):
    text = "\n".join(re.sub(r"//.*$", "", l) for l in open(path).read().splitlines())
    return json.loads(text)

def set_key(cfg, dotted, value):
    node = cfg
    *parents, leaf = dotted.split(".")
    for p in parents:
        node = node[p]
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
