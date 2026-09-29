#include "biology_parameters_config.h"

#include "../World/Usr_World/biomaterialWorld.h"
#include "Usr_Agents/Cell.h"

#include <nlohmann/json.hpp>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {

using json = nlohmann::json;

void reject_unknown_keys(const json &actual, const json &allowed,
                         const std::string &context) {
  for (auto it = actual.begin(); it != actual.end(); ++it) {
    if (it.key() == "description")
      continue;
    if (!allowed.contains(it.key()))
      throw std::invalid_argument("biology parameters config: unknown key '" +
                                  it.key() + "' in " + context);
    if (it.value().is_object() && allowed.at(it.key()).is_object())
      reject_unknown_keys(it.value(), allowed.at(it.key()),
                          context + "." + it.key());
  }
}

template <typename T>
T parse_strict(const json &obj, const std::string &context) {
  reject_unknown_keys(obj, json(T{}), context);
  return obj.get<T>();
}

void log_json_leaves(const json &node, const std::string &prefix) {
  for (auto it = node.begin(); it != node.end(); ++it) {
    const std::string key = prefix.empty() ? it.key() : prefix + "." + it.key();
    if (it.value().is_object())
      log_json_leaves(it.value(), key);
    else
      std::cout << "[biology_parameters] " << key << " = " << it.value()
                << std::endl;
  }
}

void flatten_json_leaves(const json &node, const std::string &prefix,
                         json &out) {
  for (auto it = node.begin(); it != node.end(); ++it) {
    const std::string key = prefix.empty() ? it.key() : prefix + "." + it.key();
    if (it.value().is_object())
      flatten_json_leaves(it.value(), key, out);
    else
      out[key] = it.value();
  }
}

} // namespace

BiologyParametersConfig
load_biology_parameters_config(const std::string &path) {
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("Cannot open biology parameters config: " + path);

  json root;
  try {
    root = json::parse(in, /*callback=*/nullptr, /*allow_exceptions=*/true,
                       /*ignore_comments=*/true);
  } catch (const json::parse_error &e) {
    throw std::runtime_error(std::string("Invalid JSON in ") + path + ": " +
                             e.what());
  }

  if (!root.contains("biology"))
    throw std::runtime_error(
        "biology parameters config: missing 'biology' section in " + path);

  return parse_strict<BiologyParametersConfig>(root.at("biology"), "biology");
}

void apply_biology_parameters(const BiologyParametersConfig &cfg) {
  /* ---- Agent rules: Manuscript Table 3, parameters k1 .. k59 ---- */
  const auto &fb = cfg.fibroblast;

  // k1, k2
  Fibroblast::migration[Fibroblast::MIGRATION_ELASTICITY_EFFECT] =
      static_cast<float>(fb.migration.elasticity_effect);
  Fibroblast::migration[Fibroblast::MIGRATION_BASELINE_SPEED] =
      static_cast<float>(fb.migration.baseline_speed);

  // k3, k4
  Fibroblast::viability[Fibroblast::VIABILITY_TIME_EFFECT] =
      static_cast<float>(fb.viability.time_effect);
  Fibroblast::viability[Fibroblast::VIABILITY_BIOMATERIAL_EFFECT] =
      static_cast<float>(fb.viability.biomaterial_effect);

  // k5 .. k14
  Fibroblast::proliferation[Fibroblast::PROLIFERATION_TGF_THRESHOLD] =
      static_cast<float>(fb.proliferation.tgf_threshold);
  Fibroblast::proliferation[Fibroblast::PROLIFERATION_HOURS_BETWEEN] =
      static_cast<float>(fb.proliferation.hours_between);
  Fibroblast::proliferation[Fibroblast::PROLIFERATION_TIME_EFFECT] =
      static_cast<float>(fb.proliferation.time_effect);
  Fibroblast::proliferation[Fibroblast::PROLIFERATION_HA_TIME_EFFECT] =
      static_cast<float>(fb.proliferation.ha_time_effect);
  Fibroblast::proliferation[Fibroblast::PROLIFERATION_HA_EFFECT] =
      static_cast<float>(fb.proliferation.ha_effect);
  Fibroblast::proliferation[Fibroblast::PROLIFERATION_BIOMATERIAL_BASELINE] =
      static_cast<float>(fb.proliferation.biomaterial_baseline);
  Fibroblast::proliferation[Fibroblast::PROLIFERATION_BASELINE] =
      static_cast<float>(fb.proliferation.baseline);
  Fibroblast::proliferation[Fibroblast::PROLIFERATION_ACTIVATING_FACTOR_EFFECT] =
      static_cast<float>(fb.proliferation.activating_factor_effect);
  Fibroblast::proliferation[Fibroblast::PROLIFERATION_BASELINE_OFFSET] =
      static_cast<float>(fb.proliferation.baseline_offset);
  Fibroblast::proliferation[Fibroblast::PROLIFERATION_CHEMICAL_EFFECT] =
      static_cast<float>(fb.proliferation.chemical_effect);

  // k15 .. k17
  Fibroblast::tgfSynthesis[Fibroblast::TGF_BASELINE] =
      static_cast<float>(fb.tgf_synthesis.baseline);
  Fibroblast::tgfSynthesis[Fibroblast::TGF_FEEDBACK_EFFECT] =
      static_cast<float>(fb.tgf_synthesis.feedback_effect);
  Fibroblast::tgfSynthesis[Fibroblast::TGF_FEEDBACK_BASELINE] =
      static_cast<float>(fb.tgf_synthesis.feedback_baseline);

  // k18
  Fibroblast::fgfSynthesis[Fibroblast::FGF_RATE] =
      static_cast<float>(fb.fgf_synthesis.rate);

  // k19 .. k21
  Fibroblast::tnfSynthesis[Fibroblast::TNF_BASELINE] =
      static_cast<float>(fb.tnf_synthesis.baseline);
  Fibroblast::tnfSynthesis[Fibroblast::TNF_INHIBITORY_EFFECT] =
      static_cast<float>(fb.tnf_synthesis.inhibitory_effect);
  Fibroblast::tnfSynthesis[Fibroblast::TNF_IL10_EFFECT] =
      static_cast<float>(fb.tnf_synthesis.il10_effect);

  // k22 .. k24
  Fibroblast::il6Synthesis[Fibroblast::IL6_BASELINE] =
      static_cast<float>(fb.il6_synthesis.baseline);
  Fibroblast::il6Synthesis[Fibroblast::IL6_ACTIVATING_EFFECT] =
      static_cast<float>(fb.il6_synthesis.activating_effect);
  Fibroblast::il6Synthesis[Fibroblast::IL6_INHIBITORY_EFFECT] =
      static_cast<float>(fb.il6_synthesis.inhibitory_effect);

  // k25, k26
  Fibroblast::il8Synthesis[Fibroblast::IL8_BASELINE] =
      static_cast<float>(fb.il8_synthesis.baseline);
  Fibroblast::il8Synthesis[Fibroblast::IL8_INHIBITORY_EFFECT] =
      static_cast<float>(fb.il8_synthesis.inhibitory_effect);

  // k27 .. k31
  Fibroblast::activation_params[Fibroblast::ACTIVATION_TGF_INTERMEDIATE_THRESHOLD] =
      static_cast<float>(fb.activation.tgf_intermediate_threshold);
  Fibroblast::activation_params[Fibroblast::ACTIVATION_INTERMEDIATE_PROBABILITY] =
      static_cast<float>(fb.activation.intermediate_probability);
  Fibroblast::activation_params[Fibroblast::ACTIVATION_TGF_HIGH_THRESHOLD] =
      static_cast<float>(fb.activation.tgf_high_threshold);
  Fibroblast::activation_params[Fibroblast::ACTIVATION_INDEPENDENT_PROBABILITY] =
      static_cast<float>(fb.activation.independent_probability);
  Fibroblast::activation_params[Fibroblast::ACTIVATION_DEACTIVATION_PROBABILITY] =
      static_cast<float>(fb.activation.deactivation_probability);

  // k32 .. k43
  Fibroblast::collagenSynth[Fibroblast::COLLAGEN_HOURS_BETWEEN] =
      static_cast<float>(fb.collagen_synthesis.hours_between);
  Fibroblast::collagenSynth[Fibroblast::COLLAGEN_BASELINE] =
      static_cast<float>(fb.collagen_synthesis.baseline);
  Fibroblast::collagenSynth[Fibroblast::COLLAGEN_MESH_EFFECT] =
      static_cast<float>(fb.collagen_synthesis.mesh_effect);
  Fibroblast::collagenSynth[Fibroblast::COLLAGEN_ELASTICITY_EFFECT] =
      static_cast<float>(fb.collagen_synthesis.elasticity_effect);
  Fibroblast::collagenSynth[Fibroblast::COLLAGEN_ACTIVATING_STRENGTH] =
      static_cast<float>(fb.collagen_synthesis.activating_strength);
  Fibroblast::collagenSynth[Fibroblast::COLLAGEN_ACTIVATING_OFFSET] =
      static_cast<float>(fb.collagen_synthesis.activating_offset);
  Fibroblast::collagenSynth[Fibroblast::COLLAGEN_IL8_EFFECT] =
      static_cast<float>(fb.collagen_synthesis.il8_effect);
  Fibroblast::collagenSynth[Fibroblast::COLLAGEN_HIGH_FRAG_RATIO_EFFECT] =
      static_cast<float>(fb.collagen_synthesis.high_frag_ratio_effect);
  Fibroblast::collagenSynth[Fibroblast::COLLAGEN_EQUAL_FRAG_RATIO_EFFECT] =
      static_cast<float>(fb.collagen_synthesis.equal_frag_ratio_effect);
  Fibroblast::collagenSynth[Fibroblast::COLLAGEN_EQUAL_FRAG_OFFSET] =
      static_cast<float>(fb.collagen_synthesis.equal_frag_offset);
  Fibroblast::collagenSynth[Fibroblast::COLLAGEN_BASELINE_RATE] =
      static_cast<float>(fb.collagen_synthesis.baseline_rate);
  Fibroblast::collagenSynth[Fibroblast::COLLAGEN_CHEMICAL_EFFECT] =
      static_cast<float>(fb.collagen_synthesis.chemical_effect);

  // k44 .. k50
  Fibroblast::elastinSynth[Fibroblast::ELASTIN_BASELINE] =
      static_cast<float>(fb.elastin_synthesis.baseline);
  Fibroblast::elastinSynth[Fibroblast::ELASTIN_MESH_EFFECT] =
      static_cast<float>(fb.elastin_synthesis.mesh_effect);
  Fibroblast::elastinSynth[Fibroblast::ELASTIN_ELASTICITY_EFFECT] =
      static_cast<float>(fb.elastin_synthesis.elasticity_effect);
  Fibroblast::elastinSynth[Fibroblast::ELASTIN_ACTIVATING_STRENGTH] =
      static_cast<float>(fb.elastin_synthesis.activating_strength);
  Fibroblast::elastinSynth[Fibroblast::ELASTIN_ACTIVATING_OFFSET] =
      static_cast<float>(fb.elastin_synthesis.activating_offset);
  Fibroblast::elastinSynth[Fibroblast::ELASTIN_INHIBITORY_STRENGTH] =
      static_cast<float>(fb.elastin_synthesis.inhibitory_strength);
  Fibroblast::elastinSynth[Fibroblast::ELASTIN_BASELINE_RATE] =
      static_cast<float>(fb.elastin_synthesis.baseline_rate);

  // k51 .. k59
  Fibroblast::haSynth[Fibroblast::HA_HOURS_BETWEEN] =
      static_cast<float>(fb.ha_synthesis.hours_between);
  Fibroblast::haSynth[Fibroblast::HA_BASELINE] =
      static_cast<float>(fb.ha_synthesis.baseline);
  Fibroblast::haSynth[Fibroblast::HA_MESH_EFFECT] =
      static_cast<float>(fb.ha_synthesis.mesh_effect);
  Fibroblast::haSynth[Fibroblast::HA_ELASTICITY_EFFECT] =
      static_cast<float>(fb.ha_synthesis.elasticity_effect);
  Fibroblast::haSynth[Fibroblast::HA_ACTIVATING_STRENGTH] =
      static_cast<float>(fb.ha_synthesis.activating_strength);
  Fibroblast::haSynth[Fibroblast::HA_ACTIVATING_OFFSET] =
      static_cast<float>(fb.ha_synthesis.activating_offset);
  Fibroblast::haSynth[Fibroblast::HA_MOVE_PROBABILITY] =
      static_cast<float>(fb.ha_synthesis.move_probability);
  Fibroblast::haSynth[Fibroblast::HA_SAME_PATCH_PROBABILITY] =
      static_cast<float>(fb.ha_synthesis.same_patch_probability);
  Fibroblast::haSynth[Fibroblast::HA_SAME_PATCH_EFFECT] =
      static_cast<float>(fb.ha_synthesis.same_patch_effect);

  /* ---- Biomaterial rules: Manuscript Table 4, parameters c1 .. c19 ---- */
  const auto &bm = cfg.biomaterial;

  // c1 .. c6
  BMWorld::ElasticMod[BMWorld::ELASTIC_POLYMER_HA_INTERACTION] =
      static_cast<float>(bm.elastic_modulus.polymer_ha_interaction);
  BMWorld::ElasticMod[BMWorld::ELASTIC_HA_CONCENTRATION] =
      static_cast<float>(bm.elastic_modulus.ha_concentration);
  BMWorld::ElasticMod[BMWorld::ELASTIC_POLYMER_CONCENTRATION] =
      static_cast<float>(bm.elastic_modulus.polymer_concentration);
  BMWorld::ElasticMod[BMWorld::ELASTIC_HA_CROSSLINKER_INTERACTION] =
      static_cast<float>(bm.elastic_modulus.ha_crosslinker_interaction);
  BMWorld::ElasticMod[BMWorld::ELASTIC_CROSSLINKER_CONCENTRATION] =
      static_cast<float>(bm.elastic_modulus.crosslinker_concentration);
  BMWorld::ElasticMod[BMWorld::ELASTIC_BASELINE] =
      static_cast<float>(bm.elastic_modulus.baseline);

  // c7, c8
  BMWorld::XLDensity[BMWorld::XLDENSITY_BASELINE] =
      static_cast<float>(bm.crosslink_density.baseline);
  BMWorld::XLDensity[BMWorld::XLDENSITY_THIOL_DOUBLE_BOND_EFFECT] =
      static_cast<float>(bm.crosslink_density.thiol_double_bond_effect);

  // c9 .. c12
  BMWorld::SwellRatio[BMWorld::SWELL_HA_TIME_EFFECT] =
      static_cast<float>(bm.swell_ratio.ha_time_effect);
  BMWorld::SwellRatio[BMWorld::SWELL_TIME_EFFECT] =
      static_cast<float>(bm.swell_ratio.time_effect);
  BMWorld::SwellRatio[BMWorld::SWELL_HA_EFFECT] =
      static_cast<float>(bm.swell_ratio.ha_effect);
  BMWorld::SwellRatio[BMWorld::SWELL_BASELINE] =
      static_cast<float>(bm.swell_ratio.baseline);

  // c13 .. c16
  BMWorld::MassLoss[BMWorld::MASSLOSS_HA_TIME_EFFECT] =
      static_cast<float>(bm.mass_loss.ha_time_effect);
  BMWorld::MassLoss[BMWorld::MASSLOSS_TIME_EFFECT] =
      static_cast<float>(bm.mass_loss.time_effect);
  BMWorld::MassLoss[BMWorld::MASSLOSS_HA_EFFECT] =
      static_cast<float>(bm.mass_loss.ha_effect);
  BMWorld::MassLoss[BMWorld::MASSLOSS_BASELINE] =
      static_cast<float>(bm.mass_loss.baseline);

  // c17 .. c19
  BMWorld::PoreSize[BMWorld::PORE_HA_QUADRATIC_EFFECT] =
      static_cast<float>(bm.pore_size.ha_quadratic_effect);
  BMWorld::PoreSize[BMWorld::PORE_HA_LINEAR_EFFECT] =
      static_cast<float>(bm.pore_size.ha_linear_effect);
  BMWorld::PoreSize[BMWorld::PORE_BASELINE] =
      static_cast<float>(bm.pore_size.baseline);
}

void log_biology_parameters(const BiologyParametersConfig &cfg,
                            const std::string &source_path) {
  std::cout << "[biology_parameters] source: " << source_path << std::endl;
  log_json_leaves(json(cfg), "");
}

void record_biology_parameters_in_run_params(const BiologyParametersConfig &cfg,
                                             const std::string &source_path) {
  const char *path = std::getenv("IVDBM_RUN_PARAMS_JSON");
  if (path == nullptr || path[0] == '\0')
    return;

  json root;
  {
    std::ifstream in(path);
    if (in)
      in >> root;
  }

  json values = json::object();
  flatten_json_leaves(json(cfg), "", values);

  if (!root.contains("simulation"))
    root["simulation"] = json::object();
  root["simulation"]["simulation_config"] = source_path;
  root["biology_parameters"] = values;

  std::ofstream out(path);
  if (!out) {
    std::fprintf(stderr, "Warning: cannot update run params at %s\n", path);
    return;
  }
  out << root.dump(2) << std::endl;
}
