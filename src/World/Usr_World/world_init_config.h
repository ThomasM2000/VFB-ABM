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
  /** Thiol : double bond molar ratio (Table 2 TDB_MR); not in the manuscript. */
  double thiol_double_bond_molar_ratio = 0.1;
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
   * Fibroblast agents seeded at t = 0 (Manuscript Table 1: 50 000).
   * 0 falls back to the in vitro density of 1e6 cells/mL (Section 2.2.3).
   */
  int fibroblast_count = 0;
  double initial_viability = 1.0;
  WorldInitBiomaterialParams biomaterial;
  WorldInitInitialEcmParams initial_ecm;
};

WorldInitParams load_world_init_config(const std::string &path);

#endif
