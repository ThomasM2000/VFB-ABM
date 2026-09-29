#ifndef IVDBM_BIOLOGY_PARAMETERS_CONFIG_H
#define IVDBM_BIOLOGY_PARAMETERS_CONFIG_H

/**
 * @file biology_parameters_config.h
 * @brief JSON-backed cell rule and hydrogel calibration parameters.
 *
 * Loaded from the @c biology section of simulation_config.json.
 */

#include <string>

#include <nlohmann/json.hpp>


/* -------------------------------------------------------------------------- */
/*                       AGENT RULES  (Manuscript Table 3)                     */
/* -------------------------------------------------------------------------- */

/** Migration speed: v = -k1 ln(E) + k2   [um/min] */
struct MigrationParams {
  double elasticity_effect = 0; // k1
  double baseline_speed = 0;    // k2
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(MigrationParams, elasticity_effect,
                                   baseline_speed)

/** Viability rate: vr = k3 ln(t) + k4 */
struct ViabilityParams {
  double time_effect = 0;         // k3
  double biomaterial_effect = 0;  // k4
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ViabilityParams, time_effect,
                                   biomaterial_effect)

/**
 * Proliferation (evaluated every k6 hours; x = -1 if TGF > k5 else +1):
 *   Prolif_b = (k7 - k8*HAww) log(t) + k9*HAww - k10
 *   Prolif   = Prolif_b * ( k11 + (HAf + k12 log10(1+TNF+xTGF+FGF+HAf) + k13)/k14 )
 */
struct ProliferationParams {
  double tgf_threshold = 0;            // k5
  double hours_between = 0;            // k6
  double time_effect = 0;              // k7
  double ha_time_effect = 0;           // k8
  double ha_effect = 0;                // k9
  double biomaterial_baseline = 0;     // k10
  double baseline = 0;                 // k11
  double activating_factor_effect = 0; // k12
  double baseline_offset = 0;          // k13
  double chemical_effect = 0;          // k14
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ProliferationParams, tgf_threshold,
                                   hours_between, time_effect, ha_time_effect,
                                   ha_effect, biomaterial_baseline, baseline,
                                   activating_factor_effect, baseline_offset,
                                   chemical_effect)

/** TGF-beta synthesis: TGFinc = k15 + k16 (k17 + TNF + IL10) */
struct TgfSynthesisParams {
  double baseline = 0;          // k15
  double feedback_effect = 0;   // k16
  double feedback_baseline = 0; // k17
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(TgfSynthesisParams, baseline,
                                   feedback_effect, feedback_baseline)

/** FGF synthesis: FGFinc = k18 */
struct FgfSynthesisParams {
  double rate = 0; // k18
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(FgfSynthesisParams, rate)

/** TNF-alpha synthesis: TNFinc = k19 / (k20 + TGF + k21 IL10 + HA) */
struct TnfSynthesisParams {
  double baseline = 0;          // k19
  double inhibitory_effect = 0; // k20
  double il10_effect = 0;       // k21
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(TnfSynthesisParams, baseline,
                                   inhibitory_effect, il10_effect)

/** IL-6 synthesis: IL6inc = k22 (k23 + TNF) / (k24 + IL10) */
struct Il6SynthesisParams {
  double baseline = 0;          // k22
  double activating_effect = 0; // k23
  double inhibitory_effect = 0; // k24
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Il6SynthesisParams, baseline,
                                   activating_effect, inhibitory_effect)

/** IL-8 synthesis: IL8inc = k25 / (k26 + IL10 + HA) */
struct Il8SynthesisParams {
  double baseline = 0;          // k25
  double inhibitory_effect = 0; // k26
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Il8SynthesisParams, baseline,
                                   inhibitory_effect)

/**
 * Activation / deactivation:
 *   if TGF > k27 activate with probability k28,
 *   OR if TGF > k29 activate,
 *   OR activate with probability k30.
 *   Activated fibroblasts deactivate with probability k31.
 */
struct ActivationParams {
  double tgf_intermediate_threshold = 0; // k27
  double intermediate_probability = 0;   // k28
  double tgf_high_threshold = 0;         // k29
  double independent_probability = 0;    // k30
  double deactivation_probability = 0;   // k31
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ActivationParams, tgf_intermediate_threshold,
                                   intermediate_probability,
                                   tgf_high_threshold, independent_probability,
                                   deactivation_probability)

/**
 * Collagen synthesis (every k32 hours):
 *   Col_b = k33 - k34 mesh - k35 E
 *   Col   = Col_b k36 (log10(1+TGF+IL6) + k37) / (1 + FGF + k38 IL8)
 *   if HAf > HAn      -> rate = Col + k39
 *   else if HAf == HAn-> rate = Col/k40 + k41
 *   else              -> rate = k42 + Col/k43
 */
struct CollagenSynthesisParams {
  double hours_between = 0;            // k32
  double baseline = 0;                 // k33
  double mesh_effect = 0;              // k34
  double elasticity_effect = 0;        // k35
  double activating_strength = 0;      // k36
  double activating_offset = 0;        // k37
  double il8_effect = 0;               // k38
  double high_frag_ratio_effect = 0;   // k39
  double equal_frag_ratio_effect = 0;  // k40
  double equal_frag_offset = 0;        // k41
  double baseline_rate = 0;            // k42
  double chemical_effect = 0;          // k43
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(CollagenSynthesisParams, hours_between,
                                   baseline, mesh_effect, elasticity_effect,
                                   activating_strength, activating_offset,
                                   il8_effect, high_frag_ratio_effect,
                                   equal_frag_ratio_effect, equal_frag_offset,
                                   baseline_rate, chemical_effect)

/**
 * Elastin synthesis (also every k32 hours):
 *   Eln_b = k44 - k45 mesh - k46 E
 *   Eln   = Eln_b ( (k47 log10(1+TGF) + k48) / (k49 (1+FGF+TNF)) + k50 )
 */
struct ElastinSynthesisParams {
  double baseline = 0;            // k44
  double mesh_effect = 0;         // k45
  double elasticity_effect = 0;   // k46
  double activating_strength = 0; // k47
  double activating_offset = 0;   // k48
  double inhibitory_strength = 0; // k49
  double baseline_rate = 0;       // k50
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ElastinSynthesisParams, baseline,
                                   mesh_effect, elasticity_effect,
                                   activating_strength, activating_offset,
                                   inhibitory_strength, baseline_rate)

/**
 * Hyaluronic acid synthesis (every k51 hours):
 *   HA_b = k52 - k53 mesh + k54 E
 *   HA   = HA_b (k55 log10(1+TGF+TNF+FGF) + k56)
 *   with probability HA + k57       -> move, then deposit HA
 *   with probability k58 + HA/k59   -> deposit HA on the current patch
 */
struct HaSynthesisParams {
  double hours_between = 0;         // k51
  double baseline = 0;              // k52
  double mesh_effect = 0;           // k53
  double elasticity_effect = 0;     // k54
  double activating_strength = 0;   // k55
  double activating_offset = 0;     // k56
  double move_probability = 0;      // k57
  double same_patch_probability = 0;// k58
  double same_patch_effect = 0;     // k59
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(HaSynthesisParams, hours_between, baseline,
                                   mesh_effect, elasticity_effect,
                                   activating_strength, activating_offset,
                                   move_probability, same_patch_probability,
                                   same_patch_effect)

/** All fibroblast agent rules. */
struct FibroblastParams {
  MigrationParams migration;
  ViabilityParams viability;
  ProliferationParams proliferation;
  TgfSynthesisParams tgf_synthesis;
  FgfSynthesisParams fgf_synthesis;
  TnfSynthesisParams tnf_synthesis;
  Il6SynthesisParams il6_synthesis;
  Il8SynthesisParams il8_synthesis;
  ActivationParams activation;
  CollagenSynthesisParams collagen_synthesis;
  ElastinSynthesisParams elastin_synthesis;
  HaSynthesisParams ha_synthesis;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(FibroblastParams, migration, viability,
                                   proliferation, tgf_synthesis, fgf_synthesis,
                                   tnf_synthesis, il6_synthesis, il8_synthesis,
                                   activation, collagen_synthesis,
                                   elastin_synthesis, ha_synthesis)

/* -------------------------------------------------------------------------- */
/*                    BIOMATERIAL RULES  (Manuscript Table 4)                  */
/* -------------------------------------------------------------------------- */

/** E = c1 TPwv HAww + c2 HAww + c3 TPwv + c4 HAwv XLww + c5 XLww + c6  [Pa] */
struct ElasticModulusParams {
  double polymer_ha_interaction = 0;    // c1
  double ha_concentration = 0;          // c2
  double polymer_concentration = 0;     // c3
  double ha_crosslinker_interaction = 0;// c4
  double crosslinker_concentration = 0; // c5
  double baseline = 0;                  // c6
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ElasticModulusParams, polymer_ha_interaction,
                                   ha_concentration, polymer_concentration,
                                   ha_crosslinker_interaction,
                                   crosslinker_concentration, baseline)

/** rho_XL = c7 - c8 TDB_MR   [mmol/mL] */
struct CrosslinkDensityParams {
  double baseline = 0;                 // c7
  double thiol_double_bond_effect = 0; // c8
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(CrosslinkDensityParams, baseline,
                                   thiol_double_bond_effect)

/** Q = (c9 HAww + c10) ln(t_min) + c11 HAww + c12   [% w/w] */
struct SwellRatioParams {
  double ha_time_effect = 0; // c9
  double time_effect = 0;    // c10
  double ha_effect = 0;      // c11
  double baseline = 0;       // c12
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SwellRatioParams, ha_time_effect,
                                   time_effect, ha_effect, baseline)

/** w_l = (c13 HAww - c14) t_weeks + c15 HAww + c16   [%] */
struct MassLossParams {
  double ha_time_effect = 0; // c13
  double time_effect = 0;    // c14
  double ha_effect = 0;      // c15
  double baseline = 0;       // c16
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(MassLossParams, ha_time_effect, time_effect,
                                   ha_effect, baseline)

/** p = -c17 HAww^2 + c18 HAww + c19   [um] */
struct PoreSizeParams {
  double ha_quadratic_effect = 0; // c17
  double ha_linear_effect = 0;    // c18
  double baseline = 0;            // c19
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PoreSizeParams, ha_quadratic_effect,
                                   ha_linear_effect, baseline)

struct BiomaterialParams {
  ElasticModulusParams elastic_modulus;
  CrosslinkDensityParams crosslink_density;
  SwellRatioParams swell_ratio;
  MassLossParams mass_loss;
  PoreSizeParams pore_size;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(BiomaterialParams, elastic_modulus,
                                   crosslink_density, swell_ratio, mass_loss,
                                   pore_size)

struct BiologyParametersConfig {
  FibroblastParams fibroblast;
  BiomaterialParams biomaterial;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(BiologyParametersConfig, fibroblast, biomaterial)

BiologyParametersConfig load_biology_parameters_config(const std::string &path);

void apply_biology_parameters(const BiologyParametersConfig &cfg);

void log_biology_parameters(const BiologyParametersConfig &cfg,
                            const std::string &source_path);

void record_biology_parameters_in_run_params(const BiologyParametersConfig &cfg,
                                             const std::string &source_path);

#endif
