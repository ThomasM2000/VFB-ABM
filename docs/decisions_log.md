# Decisions log (VFB-ABM)

| Date | Decision | Reason / source |
|---|---|---|
| 2026-10-14 | k3 = -2.75, k4 = 92.8 from OLS fit of viability vs ln(t), pooled gels (n=26, R²=0.34) | Wet-lab LIVE/DEAD; treated as input fit, not validation output (Nicole) |
| 2026-10-14 | t < 1 day → t = 1 in vr(t) and Prolif_b, so vr = k4 | ln t negative below day 1 |
| 2026-10-14 | Day-0 dead fraction seeded from measured day-0 viability (87.8–90.3%) | Nicole |
| 2026-10-14 | Death mechanism: P_death = 1 − vr_target/vr_current when current > target | Live = total × vr(t) by construction; Total cells carries growth info |
| 2026-10-14 | Total cells = live + cumulative dead; compared with PicoGreen | PicoGreen measures all DNA; GH2 day-9 drop reported as limit of count definition |
| 2026-10-14 | Exclude PicoGreen triplets GH2 d6 and GH10 d3; keep viability and collagen there | Identical to every decimal, either could be the copy |
| 2026-10-14 | GH10 d9 n = 2; Q4 restored to formula (11,120); SD factor = 0.6³/1000/0.3 | Spreadsheet errors |
| 2026-10-14 | Collagen: start each gel at its day-0 Sircol value, calibrate on change from day 0; no degradation rule | Sircol = absolute µg per 0.3 mL gel; decline only in GH5, unexplained as in v13 |
| 2026-10-14 | ECM unit = pg; time terms use natural log | Fit wet data; ln gives day 1.67 crossing |
| 2026-10-14 | HAf ≡ 0 in Study 1; k39 inert | No fragmentation rule; equal/else collagen branches only |
| 2026-10-14 | IL-10 has no source; k21 inert | Figure 1C lists five synthesized cytokines |
| 2026-10-14 | ρXL output only; c7, c8 inert. No time term in E, ρXL, p | Source fits have none |
| 2026-10-14 | Keep activation/deactivation (k27–k31) | Nicole; overrides handover line to remove them |
| 2026-10-14 | Bounds fix A: k7–k10, k35, c13–c19 from Emily's Sensitivity Analysis.xlsx | Supp. Table 1 bounds shifted one row from c13 |
| 2026-10-14 | Bounds fix C: k33, k34, k44, k45, k52, k53 → [X/10, X×10] | Original bounds too wide for SA |
| 2026-10-14 | Bounds fix D: k13, k37, k48, k56 → [0, 10] | Estimate 0 lay outside [0.1, 10]; set by calibration |
| 2026-10-14 | Each gel seeded from its own day-0 PicoGreen mean; Table 1's 50,000 not used | Handover; GH5 starts at 57% of GH2 |
| PENDING | k27 = 25 (intermediate), k29 = 500 (high) | Names were inverted in v13 |
| PENDING | Bounds convention for unstated ranges: ×0.25/×2.25 or ×0.1/×10 | To settle at meeting |
| 2026-10-19 | Migration verification judged at tick 24 with On = 2× v (k2 = 1.73) | Mean displacement saturates at world scale in the 0.6 mm world before day 1; 10× run saturated |
| 2026-10-20 | k27 = 25 (intermediate), k29 = 500 (high); should_activate checks k29 before k27 | Names were inverted in v13; with the old check order the high threshold could never be reached |
| 2026-10-20 | k28/k30/k31 stored as per-tick probabilities in [0, 1] (0.005/0.002/0.002); k57/k58 as fractions (0.5/0.05), bounds divided by 100 | Nicole: probabilities as fractions. Values chosen to keep behaviour identical; HA term in the k57/k58 rules keeps its original percent scale (units not defined in v13) |
| 2026-10-20 | Seeding set by day0_total_cells_per_mL (day-0 PicoGreen mean / 0.3 mL) x initial_viability x world volume, replacing fibroblast_count | Absolute counts were for the 3 mm world and over-seeded the 0.6 mm world ~125x; the measured density works for any world size |
| 2026-10-20 | All fibroblasts synthesise cytokines and ECM at the Table 3 rates; activated cells multiply the amounts (not the event probabilities) by synthesis_multiplier, fixed at 1 in Study 1 | Measured collagen rises 2.4–3.1× by day 3, which the activation gate could not produce; Liao 2008 cells behind Col_b/Eln_b/HA_b were not TGF-stimulated (Nicole, Oct 14). With multiplier = 1, k27–k31 do not affect cell counts or collagen and are inert for the Study 1 SA |
| 2026-10-20 | Growth cap uses N = live + cumulative dead cells; K = 5e5 cells/mL × world volume | K read from day-9 PicoGreen; N = total so total cells plateau at K. Verified (0.6 mm, 5 seeds, jobs 24633895/97/98): GH5 PASS, GH2 PASS, GH10 CHECK as expected (starts above K). Day-9 fold change GH2 1.06 / GH5 1.49 / GH10 1.00 vs data 1.07 / 1.88 / 0.97. No growth before day 3 in any gel: Prolif_b (k7–k10) is the calibration target |
