#!/usr/bin/env python3
"""Write a Notion-ready Markdown record of a verify_rules.py output folder.

For one folder (e.g. output/verify_<jobid>) it prints:
  1. a run-configuration row (job, commit, gel, world, ticks, seeding, day-0
     viability/collagen, k3/k4, seeds, folder),
  2. a results table in the Oct 14 verification format
     (Rule | Output | Base | Off | On | Result) at the judging day,
  3. the same numbers at every reported day,
  4. the Off/On overrides of each rule, as a footnote draft.
Paste the output straight into Notion (Markdown tables paste as tables).

The verdict uses the same rule as verify_rules.py: Off and On must both move
the output away from Base in the expected direction by more than
2 x max(SD_case, SD_base). SD is the population SD over seeds.

Usage:
  python3 scripts/report_run.py output/verify_24631305
  python3 scripts/report_run.py output/verify_24631305 --exp E4 --judge-tick 24 --ticks 24 48
  python3 scripts/report_run.py output/verify_24631305 > output/verify_24631305/notion.md
"""
import argparse, csv, glob, json, math, os, re, statistics as st, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_rules import TESTS, DAYS  # noqa: E402

B = "biology.fibroblast."
SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def fmt(x):
    """3 significant figures; x10^n for very large or very small values."""
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    if -2 <= e < 4:
        return f"{x:.{max(0, 2 - e)}f}"
    return f"{x / 10 ** e:.2f}×10{str(e).translate(SUP)}"


def flat(d, prefix=""):
    out = {}
    for k, v in d.items():
        if isinstance(v, dict):
            out.update(flat(v, prefix + k + "."))
        elif k != "description":
            out[prefix + k] = v
    return out


def series(csv_path, col_prefix):
    rows = list(csv.reader(open(csv_path)))
    idx = next((i for i, h in enumerate(rows[0]) if h.startswith(col_prefix)), None)
    if idx is None:
        sys.exit(f"column starting with '{col_prefix}' not found in {csv_path}")
    return [float(r[idx]) for r in rows[1:]]


def find_commit(root):
    """verify.sbatch logs 'commit <hash>'; fall back to a kept manifest, then git."""
    m = re.search(r"verify_(\d+)", os.path.basename(os.path.normpath(root)))
    if m:
        for log in glob.glob(f"logs/verify_{m.group(1)}.out"):
            for line in open(log):
                if line.startswith("commit "):
                    return line.split()[1], m.group(1)
        for man in glob.glob(os.path.expanduser(
                f"~/*projects/def-nicoleli/results/vfb-abm-results/*_{m.group(1)}_*/manifest.txt")) + \
                glob.glob(os.path.expanduser(
                f"~/links/projects/def-nicoleli/results/vfb-abm-results/*_{m.group(1)}_*/manifest.txt")):
            for line in open(man):
                if line.startswith("commit:"):
                    return line.split(":", 1)[1].strip(), m.group(1)
    try:
        h = subprocess.run(["git", "describe", "--always", "--dirty"], capture_output=True,
                           text=True).stdout.strip()
        if h:
            return f"{h} (current checkout, not read from the job log)", m.group(1) if m else "local"
        return "unknown", m.group(1) if m else "local"
    except OSError:
        return "unknown", m.group(1) if m else "local"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("root", help="verify_rules.py output folder")
    p.add_argument("--exp", default="E?", help="experiment label used in the Notion tables")
    p.add_argument("--judge-tick", type=int, default=DAYS[9],
                   help="tick used for the verdict (default 432 = day 9)")
    p.add_argument("--ticks", type=int, nargs="*", default=[],
                   help="extra ticks to report besides days 3/6/9")
    a = p.parse_args()

    root = a.root.rstrip("/")
    cases = {}
    for d in sorted(os.listdir(root)):
        m = re.match(r"(.+)_(base|off|on)_s(\d+)$", d)
        csv_path = os.path.join(root, d, "Output_Biomarkers.csv")
        if m and os.path.exists(csv_path):
            cases.setdefault((m.group(1), m.group(2)), []).append(int(m.group(3)))
    if ("baseline", "base") not in cases:
        sys.exit(f"no complete baseline_base_s* runs in {root}")
    rules = sorted({r for r, c in cases if r != "baseline"},
                   key=lambda r: list(TESTS).index(r) if r in TESTS else 99)

    # --- run configuration (from the first baseline seed) -------------------
    seeds = sorted(cases[("baseline", "base")])
    base_dir = os.path.join(root, f"baseline_base_s{seeds[0]}")
    cfg = flat(json.load(open(os.path.join(base_dir, "config.json"))))
    log = open(os.path.join(base_dir, "stdout.log")).read()
    seeding = re.search(r"Seeding (\d+) fibroblasts \(([\d.e+]+) cells/mL\)", log)
    nx = re.search(r"nx, ny, nz: (\d+)", log)
    n_ticks = sum(1 for _ in open(os.path.join(base_dir, "Output_Biomarkers.csv"))) - 2
    commit, job = find_commit(root)
    gel = f"GH{cfg['world_init.biomaterial.ha_gtn_ratio']:g}"
    world = f"{int(nx.group(1)) * 0.015:.2g} mm" if nx else "?"
    seed_txt = (f"{seeding.group(1)} live / {fmt(float(seeding.group(2)))} cells/mL"
                if seeding else "?")

    print(f"## {a.exp}: {', '.join(rules)} (job {job})\n")
    print("**Run configuration**\n")
    print("| Exp | Job | Commit | Gel | World | Ticks | Seeding | Day-0 viability / collagen "
          "| k3 / k4 | Seeds | Folder |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    print(f"| {a.exp} | {job} | `{commit}` | {gel} | {world} | {n_ticks} | {seed_txt} "
          f"| {cfg['world_init.initial_viability']} / "
          f"{cfg['world_init.initial_ecm.collagen_per_patch']} pg/patch "
          f"| {cfg[B + 'viability.time_effect']} / {cfg[B + 'viability.biomaterial_effect']} "
          f"| {seeds[0]}–{seeds[-1]} | `{root}` |\n")

    # --- results ----------------------------------------------------------------
    ticks = sorted(set(list(DAYS.values()) + a.ticks + [a.judge_tick]))
    tick_name = {t: (f"day {t / 48:g}" if t % 48 == 0 else f"tick {t} ({t / 2:g} h)") for t in ticks}

    def stats(rule, case, col):
        runs = [series(os.path.join(root, f"{rule}_{case}_s{s}", "Output_Biomarkers.csv"), col)
                for s in sorted(cases[(rule, case)])]
        return {t: (st.mean(r[t] for r in runs), st.pstdev([r[t] for r in runs]), len(runs))
                for t in ticks}

    summary, detail, notes = [], [], []
    for rule in rules:
        if (rule, "off") not in cases or (rule, "on") not in cases:
            notes.append(f"- {rule}: Off or On runs missing, skipped.")
            continue
        col, _, _, exp_off, exp_on = TESTS[rule]
        s = {"Base": stats("baseline", "base", col),
             "Off": stats(rule, "off", col), "On": stats(rule, "on", col)}

        def moved(case, direction):
            m, sd, _ = s[case][a.judge_tick]
            mb, sdb, _ = s["Base"][a.judge_tick]
            gap = 2 * max(sd, sdb, 1e-9)
            return m < mb - gap if direction == "lower" else m > mb + gap

        verdict = "PASS" if moved("Off", exp_off) and moved("On", exp_on) else "CHECK"
        cell = lambda c, t: f"{fmt(s[c][t][0])} ± {fmt(s[c][t][1])}"
        summary.append(f"| {rule} | {col.strip(' (')} | {cell('Base', a.judge_tick)} "
                       f"| {cell('Off', a.judge_tick)} | {cell('On', a.judge_tick)} | {verdict} |")
        for t in ticks:
            detail.append(f"| {rule} | {tick_name[t]} | {cell('Base', t)} | {cell('Off', t)} "
                          f"| {cell('On', t)} |")
        # Overrides are read from the runs' own config.json, so the footnote
        # matches what was actually run even if TESTS changes later.
        def diff(case):
            run_cfg = flat(json.load(open(os.path.join(
                root, f"{rule}_{case}_s{sorted(cases[(rule, case)])[0]}", "config.json"))))
            return ", ".join(f"{k.replace(B, '')} = {fmt(v) if isinstance(v, (int, float)) else v}"
                             for k, v in run_cfg.items() if cfg.get(k) != v) or "(none)"
        notes.append(f"- **{rule}**: Off = {diff('off')} (expected {exp_off}); "
                     f"On = {diff('on')} (expected {exp_on}).")

    judge = tick_name[a.judge_tick]
    print(f"**Results** ({a.exp}; values at {judge}, mean ± SD over {len(seeds)} seeds)\n")
    print(f"| Rule | Output | Base ({judge}) | Off ({judge}) | On ({judge}) | Result |")
    print("|---|---|---|---|---|---|")
    print("\n".join(summary) + "\n")
    print("**All reported time points**\n")
    print("| Rule | Time | Base | Off | On |")
    print("|---|---|---|---|---|")
    print("\n".join(detail) + "\n")
    print("**Off / On settings** (footnote draft)\n")
    print("\n".join(notes))
    print(f"\nVerdict rule: Off and On each differ from Base in the expected direction by more "
          f"than 2 × max(SD_case, SD_base) at {judge}.")


if __name__ == "__main__":
    main()
