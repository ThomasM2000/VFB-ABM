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
      section, {"description", "fibroblast_count", "biomaterial", "initial_ecm"},
      "world_init");

  WorldInitParams cfg;
  cfg.fibroblast_count = section.at("fibroblast_count").get<int>();
  if (cfg.fibroblast_count < 0)
    throw std::invalid_argument(
        "world_init config: fibroblast_count must be >= 0");

  const json &bm = section.at("biomaterial");
  reject_unknown_keys(bm,
                      {"description", "ha_gtn_ratio", "ha_wv_percent",
                        "gtn_wv_percent", "pegda_wv_percent",
                        "thiol_double_bond_molar_ratio"},
                      "world_init.biomaterial");
  cfg.biomaterial.ha_gtn_ratio     = bm.at("ha_gtn_ratio").get<double>();
  cfg.biomaterial.ha_wv_percent    = bm.at("ha_wv_percent").get<double>();
  cfg.biomaterial.gtn_wv_percent   = bm.at("gtn_wv_percent").get<double>();
  cfg.biomaterial.pegda_wv_percent = bm.at("pegda_wv_percent").get<double>();
  cfg.biomaterial.thiol_double_bond_molar_ratio =
      bm.at("thiol_double_bond_molar_ratio").get<double>();

  if (cfg.biomaterial.ha_gtn_ratio <= 0.0)
    throw std::invalid_argument(
        "world_init config: biomaterial.ha_gtn_ratio must be > 0");

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
