#!/usr/bin/env python3
"""Print the run configuration of a verify_rules.py output folder.

One block per case (baseline_base, <rule>_off, <rule>_on): gel, world size,
ticks, seeding, day-0 viability and collagen, k3/k4, seeds, and every config
key that differs from the baseline case.

Usage: python3 scripts/describe_run.py output/verify_<jobid>
"""
import csv, json, os, re, sys
from collections import defaultdict

B = "biology.fibroblast."

def flat(d, prefix=""):
    out = {}
    for k, v in d.items():
        if isinstance(v, dict):
            out.update(flat(v, prefix + k + "."))
        elif k != "description":
            out[prefix + k] = v
    return out

root = sys.argv[1]
cases = defaultdict(list)
for d in sorted(os.listdir(root)):
    m = re.match(r"(.+)_s(\d+)$", d)
    if m and os.path.exists(os.path.join(root, d, "config.json")):
        cases[m.group(1)].append(int(m.group(2)))

base_cfg = None
for case in sorted(cases, key=lambda c: c != "baseline_base"):
    run = os.path.join(root, f"{case}_s{cases[case][0]}")
    cfg = flat(json.load(open(os.path.join(run, "config.json"))))
    log = open(os.path.join(run, "stdout.log")).read()
    seeding = re.search(r"Seeding (\d+) fibroblasts \(([\d.e+]+) cells/mL\)", log)
    nx = re.search(r"nx, ny, nz: (\d+)", log)
    csv_path = os.path.join(run, "Output_Biomarkers.csv")
    ticks = sum(1 for _ in open(csv_path)) - 2 if os.path.exists(csv_path) else "n/a"
    if base_cfg is None:
        base_cfg = cfg
    diff = {k: v for k, v in cfg.items() if base_cfg.get(k) != v}

    print(f"== {case}")
    print(f"  gel        GH{cfg['world_init.biomaterial.ha_gtn_ratio']:g}")
    print(f"  world      {int(nx.group(1)) * 0.015:.2g} mm" if nx else "  world      ?")
    print(f"  ticks      {ticks}")
    print(f"  seeding    {seeding.group(1)} live cells ({float(seeding.group(2)):.3g} cells/mL)"
          if seeding else "  seeding    ?")
    print(f"  day 0      viability {cfg['world_init.initial_viability']}, "
          f"collagen {cfg['world_init.initial_ecm.collagen_per_patch']} pg/patch")
    print(f"  k3 / k4    {cfg[B + 'viability.time_effect']} / {cfg[B + 'viability.biomaterial_effect']}")
    print(f"  seeds      {sorted(cases[case])}")
    print(f"  overrides  {diff if diff else '(none)'}")
