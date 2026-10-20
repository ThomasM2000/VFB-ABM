#!/usr/bin/env python3
"""Compare model runs of one gel with the wet-lab data, as fold change from day 0.

Prints a Notion-ready Markdown table: one row per output x day with the model
mean ± SD over runs, the wet-lab mean ± SD and n, and z = (model - wet) / wet SD.

  Total cells : model "Total Cells (live + dead)" vs PicoGreen (cells per gel)
  Live fraction: model "Viability Rate (%)"       vs LIVE/DEAD (%)  - absolute, not a ratio
  Collagen    : model "Collagen (pg)"             vs Sircol (ug per gel), day-0 amount included

Fold change is each quantity divided by its own day-0 value. Model runs are
paired with their own day 0 (one ratio per run). Wet-lab gels are destroyed at
each time point, so day 0 and day t are different gels: FC = mean_t / mean_0 and
its SD is propagated, SD_FC = FC * sqrt((SD_t/mean_t)^2 + (SD_0/mean_0)^2).

With --world (or a stdout.log next to each CSV) the absolute values are also
shown, scaled to one 0.3 mL gel.

Usage:
  python3 scripts/compare_to_data.py GH10 output/job_*/Output_Biomarkers.csv --world 3
  python3 scripts/compare_to_data.py GH5 output/verify_123/baseline_base_s*/Output_Biomarkers.csv
"""
import argparse, csv, math, os, re, statistics as st, sys

DAYS = {0: 0, 3: 144, 6: 288, 9: 432}   # day -> tick (30 min per tick)
GEL_ML = 0.3                              # volume of one wet-lab gel

# Per-sample values from "ScaffoldABMExperiment_Summary by Sample.xlsx", sheet
# ByGroup. None = no sample. PicoGreen in cells per gel, viability as a
# fraction, Sircol in ug collagen per gel (acid-pepsin soluble, not normalised).
DATA = {
    "GH2": {
        "cells": {0: [156870.75538676282, 129967.33738398763, 135428.0453157124],
                  3: [199242.82216632584, 202328.65836520752, 199322.57430589298],
                  6: [142278.29137880268, 135489.83280950671, 137343.72981493876],
                  9: [152198.85871942536, 150815.58099492011, 147938.78794684916]},
        "via":   {0: [0.8730097087378641, 0.85763146643864896, 0.90272614622057001],
                  3: [0.91051660516605171, 0.89949263107030686, 0.89240065070880781],
                  6: [0.8707865168539326, 0.86727456940222902, 0.87936046511627908],
                  9: [0.87800888057957471, 0.86810661764705888, 0.87125748502994016]},
        "col":   {0: [40.354773510878566, 27.880870288907364, 27.7615027939603],
                  3: [69.33123290928539, 79.358102484841226, 80.909879919153425],
                  6: [102.93318273689208, 153.36594935203891, 160.40863155391742],
                  9: [107.46914754488161, 86.400784686719732, 123.04660563547725]},
    },
    "GH5": {
        "cells": {0: [73919.820105684892, 53630.930238132692, 112182.33806891956],
                  3: [95709.574470073436, 96575.960341206912, 108787.2922094637],
                  6: [130807.32066966716, 134900.81018144111, 150516.71461524541],
                  9: [148466.83965592831, 141988.13513041174, 160325.1390155894]},
        "via":   {0: [0.89815526518063027, 0.89563106796116509, 0.89636010002778554],
                  3: [0.88633093525179851, 0.89342153619992648, 0.90446260213702079],
                  6: [0.82834850455136544, 0.84651162790697676, 0.8751023751023751],
                  9: [0.85642034279291934, 0.87972620599739249, 0.86612903225806448]},
        "col":   {0: [12.959933420520724, 78.820948757579288, 32.446676970633668],
                  3: [112.06479610034471, 143.66734038758759, 124.21043871121141],
                  6: [96.636547378432923, 87.445250267506736, 86.430626560456503],
                  9: [49.635596243015051, 75.418975151587148, 84.640114136250105]},
    },
    "GH10": {
        "cells": {0: [198251.77254119329, 141202.04578205064, 154631.43507202354],
                  3: [142278.29137880268, 135489.83280950671, 137343.72981493876],
                  6: [135355.64234941598, 126707.29855943039, 152208.1132339144],
                  9: [158185.16863580575, 161892.9626466698, None]},
        "via":   {0: [0.89331366764995079, 0.91443452380952384, 0.90049621311047268],
                  3: [0.90441533236293059, 0.90111025045184612, 0.89493751323872062],
                  6: [0.88142292490118579, 0.899070385126162, 0.92934051144010765],
                  9: [0.86974894357444688, 0.86596289187460018, None]},
        "col":   {0: [41.399239091665649, 28.358340268695716, 38.444893591725091],
                  3: [64.884793722506203, 105.43990013078101, 99.411841635952854],
                  6: [81.357508025205007, 57.275115919629023, 39.638568541196022],
                  9: [68.316609202235142, 82.90928545951725, 106.06657947925325]},
    },
}

# GH2 day 6 and GH10 day 3 PicoGreen triplets are identical to every decimal;
# either could be the copy, so both are excluded (Nicole, Oct 14). Viability
# and collagen at those time points are kept.
EXCLUDED = {("GH2", "cells", 6), ("GH10", "cells", 3)}

OUTPUTS = [  # key, model column prefix, label, compare as fold change?
    ("cells", "Total Cells", "Total cells (PicoGreen)", True),
    ("via", "Viability Rate", "Live fraction % (LIVE/DEAD)", False),
    ("col", "Collagen", "Collagen (Sircol)", True),
]


def wet(gel, key, day):
    """mean, sample SD, n of the wet-lab samples (viability in %)."""
    vals = [v for v in DATA[gel][key][day] if v is not None]
    if key == "via":
        vals = [100 * v for v in vals]
    return st.mean(vals), (st.stdev(vals) if len(vals) > 1 else 0.0), len(vals)


def read_run(path):
    rows = list(csv.reader(open(path)))
    head = rows[0]
    out = {}
    for key, prefix, _, _ in OUTPUTS:
        i = next((k for k, h in enumerate(head) if h.startswith(prefix)), None)
        if i is None:
            sys.exit(f"column starting with '{prefix}' not found in {path}")
        if len(rows) <= max(DAYS.values()) + 1:
            sys.exit(f"{path} has only {len(rows) - 2} ticks; need {max(DAYS.values())}")
        out[key] = {d: float(rows[t + 1][i]) for d, t in DAYS.items()}
    return out


def world_ml(path, world_mm):
    if world_mm:
        return world_mm ** 3 / 1000.0
    log = os.path.join(os.path.dirname(path), "stdout.log")
    if os.path.exists(log):
        # "Construct volume" is printed with 3 decimals (0.000 mL for a 0.6 mm
        # world), so use the patch grid instead: 15 um patches.
        m = re.search(r"nx, ny, nz: (\d+), (\d+), (\d+)", open(log).read())
        if m:
            nx, ny, nz = (int(g) for g in m.groups())
            return nx * ny * nz * 0.015 ** 3 / 1000.0
    return None


def fmt(x):
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    return f"{x:.{max(0, 2 - e)}f}" if -2 <= e < 4 else f"{x:.3g}"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("gel", choices=sorted(DATA))
    p.add_argument("csv", nargs="+", help="Output_Biomarkers.csv of each run (one per seed)")
    p.add_argument("--world", type=float, help="world edge in mm, for absolute values")
    a = p.parse_args()

    runs = [read_run(f) for f in a.csv]
    vols = {world_ml(f, a.world) for f in a.csv}
    vol = vols.pop() if len(vols) == 1 else None
    if vol is None and (a.world or any(os.path.exists(os.path.join(os.path.dirname(f), "stdout.log")) for f in a.csv)):
        print("warning: runs have different or unknown world sizes; absolute values skipped\n", file=sys.stderr)
    to_gel = (GEL_ML / vol) if vol else None   # model total -> per 0.3 mL gel

    print(f"## {a.gel}: model vs wet lab ({len(runs)} runs)\n")
    print("| Output | Day | Model FC | Wet FC | Wet n | z | Model (per gel) | Wet (per gel) |")
    print("|---|---|---|---|---|---|---|---|")
    for key, _, label, as_fc in OUTPUTS:
        for d in DAYS:
            excluded = (a.gel, key, d) in EXCLUDED
            mw, sw, nw = wet(a.gel, key, d)
            model = [r[key][d] for r in runs]
            if as_fc:
                m_fc = [r[key][d] / r[key][0] for r in runs]
                mm, ms = st.mean(m_fc), (st.pstdev(m_fc) if len(m_fc) > 1 else 0.0)
                m0, s0, _ = wet(a.gel, key, 0)
                wfc = mw / m0
                wsd = wfc * math.sqrt((sw / mw) ** 2 + (s0 / m0) ** 2) if d else 0.0
                model_cell, wet_cell = f"{fmt(mm)} ± {fmt(ms)}", f"{fmt(wfc)} ± {fmt(wsd)}"
                z = (mm - wfc) / wsd if wsd > 0 else None
            else:
                mm, ms = st.mean(model), (st.pstdev(model) if len(model) > 1 else 0.0)
                model_cell, wet_cell = "–", "–"
                z = (mm - mw) / sw if sw > 0 else None

            if key == "via":
                abs_model, abs_wet = f"{fmt(mm)} ± {fmt(ms)} %", f"{fmt(mw)} ± {fmt(sw)} %"
            elif to_gel:
                scale = to_gel * (1e-6 if key == "col" else 1.0)          # pg -> ug for collagen
                unit = " ug" if key == "col" else ""
                ma = [v * scale for v in model]
                abs_model = f"{fmt(st.mean(ma))} ± {fmt(st.pstdev(ma) if len(ma) > 1 else 0.0)}{unit}"
                abs_wet = f"{fmt(mw)} ± {fmt(sw)}{unit}"
            else:
                abs_model, abs_wet = "–", f"{fmt(mw)} ± {fmt(sw)}"

            z_cell = "–" if z is None or excluded or (as_fc and d == 0) else f"{z:.2f}"
            n_cell = f"{nw} (excluded)" if excluded else str(nw)
            if excluded:
                wet_cell = f"~~{wet_cell}~~"
            print(f"| {label} | {d} | {model_cell} | {wet_cell} | {n_cell} | {z_cell} "
                  f"| {abs_model} | {abs_wet} |")
    print("\nFC = value / day-0 value. Wet FC SD is propagated from the day-0 and day-t SDs "
          "(independent gels); model FC is per run, then mean ± SD over runs. "
          "Live fraction is compared as an absolute %. z = (model − wet) / wet SD; |z| < 2 counts as agreement. "
          "Excluded: GH2 day 6 and GH10 day 3 PicoGreen (duplicated triplets). GH10 day 9: n = 2.")


if __name__ == "__main__":
    main()
