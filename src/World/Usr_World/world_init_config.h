#ifndef IVDBM_WORLD_INIT_CONFIG_H
#define IVDBM_WORLD_INIT_CONFIG_H

/**
 * @file world_init_config.h
 * @brief JSON-backed initial configuration of the VFB-ABM.
 *
 * Mirrors Manuscript Table 1 and the HA-Gtn-PEGDA recipe of Section 2.2.3.
 * The config carries the *recipe* as weighed out at the bench; the Table 2
 * world variables (HA_ww, HA_wv, TP_wv, XL_ww) are derived in
 * BMWorld::userInput() so their unit convention lives in one place.
 */

#include <string>

/** HA-Gtn-PEGDA hydrogel recipe (Manuscript Section 2.2.3). */
struct WorldInitBiomaterialParams {
  /** HA : gelatin volumetric ratio; 2, 5 or 10 for GH2 / GH5 / GH10. */
  double ha_gtn_ratio = 10.0;
  /** Thiolated HA (CMHA-S) stock concentration, % w/v. */
  double ha_wv_percent = 1.0;
  /** Thiolated gelatin (Gtn-DTPH) stock concentration, % w/v. */
  double gtn_wv_percent = 1.0;
  /** PEGDA crosslinker final concentration, % w/v. */
  double pegda_wv_percent = 0.25;
  /** PEGDA stock concentration before mixing, % w/v (Shung 2017 §2.2.8 p.43). */
  double pegda_stock_wv_percent = 1.67;
  /** Free thiol content of CMHA-S, mmol/g (Vanderhooft 2009, measured). */
  double ha_thiol_mmol_per_g = 0.89;
  /** Free thiol content of Gtn-DTPH, mmol/g (Vanderhooft 2009, measured). */
  double gtn_thiol_mmol_per_g = 0.57;
  /** PEGDA molecular weight, g/mol (two acrylate C=C per molecule). */
  double pegda_mw_g_per_mol = 3400.0;
};

/**
 * Day-0 ECM present in the construct, per patch. Calibration initial
 * conditions measured per experimental condition (GH2/GH5/GH10): collagen
 * from the Sircol assay (Section 2.2.5).
 */
struct WorldInitInitialEcmParams {
  double collagen_per_patch = 0.0; // Col at t = 0, ug/patch
  double elastin_per_patch  = 0.0; // Eln at t = 0, ug/patch
  double ha_per_patch       = 0.0; // HA  at t = 0, ug/patch
};

struct WorldInitParams {
    /**
   * Day-0 cell density (live + dead), cells/mL: day-0 PicoGreen mean per
   * 0.3 mL gel / 0.3. Live agents seeded = density x initial_viability x
   * construct volume, so one config works for any world size.
   */
  double day0_total_cells_per_mL = 1.0e6;
  double initial_viability = 1.0;
  WorldInitBiomaterialParams biomaterial;
  WorldInitInitialEcmParams initial_ecm;
};

WorldInitParams load_world_init_config(const std::string &path);

#endif
