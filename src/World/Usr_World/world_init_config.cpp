#include "world_init_config.h"

#include <nlohmann/json.hpp>

#include <fstream>
#include <set>
#include <stdexcept>

namespace
{

using json = nlohmann::json;

void reject_unknown_keys(const json &obj, const std::set<std::string> &allowed,
                         const std::string &context)
{
  for (auto it = obj.begin(); it != obj.end(); ++it)
  {
    if (allowed.count(it.key()) == 0)
      throw std::invalid_argument("world_init config: unknown key '" +
                                  it.key() + "' in " + context);
  }
}

} // namespace

WorldInitParams load_world_init_config(const std::string &path)
{
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("Cannot open simulation config: " + path);

  json root;
  try
  {
    root = json::parse(in, /*callback=*/nullptr, /*allow_exceptions=*/true,
                       /*ignore_comments=*/true);
  }
  catch (const json::parse_error &e)
  {
    throw std::runtime_error(std::string("Invalid JSON in ") + path + ": " +
                             e.what());
  }

  if (!root.contains("world_init"))
    throw std::runtime_error(
        "world_init config: missing 'world_init' section in " + path);

  const json &section = root.at("world_init");
  reject_unknown_keys(
      section, {"description", "day0_total_cells_per_mL", "biomaterial", "initial_ecm", "initial_viability"},
      "world_init");

  WorldInitParams cfg;
  cfg.day0_total_cells_per_mL = section.at("day0_total_cells_per_mL").get<double>();
  if (cfg.day0_total_cells_per_mL <= 0.0)
    throw std::invalid_argument(
        "world_init config: day0_total_cells_per_mL must be > 0");

  
  if (section.contains("initial_viability")) {
    cfg.initial_viability = section.at("initial_viability").get<double>();
    if (cfg.initial_viability <= 0.0 || cfg.initial_viability > 1.0)
      throw std::invalid_argument(
          "world_init config: initial_viability must be in (0, 1]");
  }


  const json &bm = section.at("biomaterial");
  reject_unknown_keys(bm,
    {"description", "ha_gtn_ratio", "ha_wv_percent",
      "gtn_wv_percent", "pegda_wv_percent",
      "pegda_stock_wv_percent",
      "ha_thiol_mmol_per_g",
      "gtn_thiol_mmol_per_g",
      "pegda_mw_g_per_mol"},
    "world_init.biomaterial");
  cfg.biomaterial.ha_gtn_ratio     = bm.at("ha_gtn_ratio").get<double>();
  cfg.biomaterial.ha_wv_percent    = bm.at("ha_wv_percent").get<double>();
  cfg.biomaterial.gtn_wv_percent   = bm.at("gtn_wv_percent").get<double>();
  cfg.biomaterial.pegda_wv_percent = bm.at("pegda_wv_percent").get<double>();
  cfg.biomaterial.pegda_stock_wv_percent = bm.at("pegda_stock_wv_percent").get<double>();
  cfg.biomaterial.ha_thiol_mmol_per_g = bm.at("ha_thiol_mmol_per_g").get<double>();
  cfg.biomaterial.gtn_thiol_mmol_per_g = bm.at("gtn_thiol_mmol_per_g").get<double>();
  cfg.biomaterial.pegda_mw_g_per_mol = bm.at("pegda_mw_g_per_mol").get<double>();

  if (cfg.biomaterial.ha_thiol_mmol_per_g <= 0.0)
    throw std::invalid_argument(
        "world_init config: biomaterial.ha_thiol_mmol_per_g must be > 0");
  if (cfg.biomaterial.gtn_thiol_mmol_per_g <= 0.0)
    throw std::invalid_argument(
        "world_init config: biomaterial.gtn_thiol_mmol_per_g must be > 0");
  if (cfg.biomaterial.pegda_mw_g_per_mol <= 0.0)
    throw std::invalid_argument(
        "world_init config: biomaterial.pegda_mw_g_per_mol must be > 0");

  if (cfg.biomaterial.ha_gtn_ratio <= 0.0)
    throw std::invalid_argument(
        "world_init config: biomaterial.ha_gtn_ratio must be > 0");

  if (cfg.biomaterial.pegda_stock_wv_percent <= cfg.biomaterial.pegda_wv_percent)
    throw std::invalid_argument(
        "world_init config: biomaterial.pegda_stock_wv_percent must be > biomaterial.pegda_wv_percent");

  /* initial_ecm is optional; struct defaults apply when absent. */
  if (section.contains("initial_ecm")) {
    const json &ecm = section.at("initial_ecm");
    reject_unknown_keys(ecm,
                        {"description", "collagen_per_patch",
                          "elastin_per_patch", "ha_per_patch"},
                        "world_init.initial_ecm");
    if (ecm.contains("collagen_per_patch"))
      cfg.initial_ecm.collagen_per_patch = ecm.at("collagen_per_patch").get<double>();
    if (ecm.contains("elastin_per_patch"))
      cfg.initial_ecm.elastin_per_patch = ecm.at("elastin_per_patch").get<double>();
    if (ecm.contains("ha_per_patch"))
      cfg.initial_ecm.ha_per_patch = ecm.at("ha_per_patch").get<double>();
  }

  return cfg;
}
