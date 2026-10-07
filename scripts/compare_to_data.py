#!/usr/bin/env python3
"""Compare Output_Biomarkers.csv runs with the ByGroup sheet (days 3/6/9).

Usage: compare_to_data.py GH10 output/job_*/Output_Biomarkers.csv
Experimental values are per 0.3 mL sample; scaled to the world volume.
"""
import csv, statistics as st, sys

DAYS = {3: 144, 6: 288, 9: 432}
SAMPLE_ML = 0.3
# mean and SD of the three samples; GH10 day 3 PicoGreen is excluded
# (identical to GH2 day 6, suspected copy error).
DATA = {
 "GH2":  {"cells": {3: (200298, 1756), 6: (138371, 3465), 9: (150318, 2138)},
          "via":   {3: (90.1, 0.9), 6: (87.2, 0.6), 9: (87.2, 0.5)},
          "col":   {3: (76.53, 6.29), 6: (138.90, 31.35), 9: (105.64, 18.39)}},
 "GH5":  {"cells": {3: (100358, 7329), 6: (138742, 10342), 9: (150260, 9293)},
          "via":   {3: (89.5, 0.9), 6: (85.0, 2.4), 9: (86.7, 1.2)},
          "col":   {3: (126.65, 15.94), 6: (90.17, 5.62), 9: (69.90, 18.14)}},
 "GH10": {"cells": {6: (138090, 12988), 9: (160039, 2622)},
          "via":   {3: (90.0, 0.5), 6: (90.3, 2.4), 9: (86.8, 0.3)},
          "col":   {3: (89.91, 21.88), 6: (59.42, 20.94), 9: (85.76, 19.04)}},
}

def col(head, prefix):
    return next(i for i, h in enumerate(head) if h.startswith(prefix))

def main():
    group, files = sys.argv[1], sys.argv[2:]
    sims = {k: {d: [] for d in DAYS} for k in ("cells", "via", "col")}
    for f in files:
        rows = list(csv.reader(open(f)))
        h = rows[0]
        world_ml = None
        for d, t in DAYS.items():
            r = rows[t + 1]
            sims["cells"][d].append(float(r[col(h, "Total Cells")]))
            sims["via"][d].append(float(r[col(h, "Viability")]))
            sims["col"][d].append(float(r[col(h, "Collagen")]))
    world_ml = float(input("world edge in mm (e.g. 3): ")) ** 3 / 1000
    scale = world_ml / SAMPLE_ML
    E = 0.0
    print(f"{'output':6s} {'day':>3s} {'sim mean':>12s} {'exp mean':>12s} {'exp SD':>10s} {'z':>7s}")
    for key, exp_scale in (("cells", scale), ("via", 1.0), ("col", scale * 1e6)):  # ug -> pg
        for d in DAYS:
            if d not in DATA[group][key]:
                continue
            m_e, sd_e = DATA[group][key][d]
            m_e, sd_e = m_e * exp_scale, sd_e * exp_scale
            m_s = st.mean(sims[key][d])
            z = (m_s - m_e) / sd_e
            if key != "via":
                E += z * z
            print(f"{key:6s} {d:3d} {m_s:12.4g} {m_e:12.4g} {sd_e:10.3g} {z:7.2f}")
    print(f"\nSD-weighted error over cells + collagen: E = {E:.2f}  ({len(files)} runs)")

if __name__ == "__main__":
    main()
