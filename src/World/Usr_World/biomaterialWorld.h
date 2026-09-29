/*
 * biomaterialWorld.h
 *
 * File Contents: Contains declarations for the BMWorld class.
 *
 * Author: Yvonna
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 */

#ifndef BMWORLD_H
#define BMWORLD_H

#include "../../Agent/Usr_Agents/Cell.h"
#include "../../Agent/biology_parameters_config.h"
#include "../../ArrayChain/ArrayChain.h"
#include "../../Chemistry/chemical_environment.h"
#include "../../ECM/ECM.h"
#include "../../common.h"
#include "../World.h"

#include <map>
#include <memory>
#include <new>
#include <stdlib.h>
#include <vector>

class Cell;
class Fibroblast;
class Collagen;
class Elastin;
class HA;
class fHA;
class Hyaluronan;
class ECM;

using namespace std;

/*
 BMWORLD (BIOMATERIAL WORLD) CLASS DESCRIPTION: BMWorld is a derived class of
 the parent class World.
 *                                                The BMWorld class manages the
 model world.
 *                                                It is used to initialize
 cells, ECM, patches, and chemicals; to destroy
 *                                                agent ArrayChains; to execute
 each timestep of the model; to sprout agents;
 *                                                to count patches; and to
 output data.
 */
class BMWorld : public World {
public:
  /*
   * Description:	BMWorld constructor.
   *
   * Return: void
   *
   * Parameters: width    -- Width (x dimension) of the world in millimeters
   *             length   -- Length (y dimension) of the world in millimeters
   *             height   -- Height (z dimension) of the world in millimeters
   *             plength  -- Length of each patch (grid point) in millimeters
   */
  BMWorld(double width = 3,     // mm
          double length = 3,    // mm
          double height = 3,    // mm
          double plength = 0.015 // mm (15 um)
  );

  /*
   * Description:	BMWorld destructor.
   *
   * Return: void
   * Parameters: void
   * NOTE: Function is called implicitly: chonds.~ArrayChain()
   */
  ~BMWorld();

  /*
   * Description:	Destructor function for the Cell ArrayChain
   *
   * Return: void
   *
   * Parameters: &agent  -- Reference to cell that will be destroyed
   */
  void destroyCell(Cell *&agent);

  /*
   * Description:	Assign a patch type to each patch within bounds of type
   *
   * Return: void
   * Parameters: void
   */
  void assignPatches(int type, int xmin, int xmax, int ymin, int ymax, int zmin,
                     int zmax);

  /*
   * Description:	Assign a patch type to each patch in the world in row
   * major index manner
   *
   * Return: void
   * Parameters: void
   */
  void initializePatches();

  /** Copy baseline_total_mass from simulation_config.json (chemistry section)
   * into baselineChem.
   */
  void sync_baseline_chem_from_config();

  /** Set baseline p* on patches via ChemicalEnvironment (after facade is
   * wired). */
  void initializeChemBaseline();

  /*
   * Description:	Initializes all cells to their correct patches
   *
   * Return: void
   * Parameters: void
   */
  void initializeCells();

  /*
   * Description:	Initializes collagen and hyaluronan to their
   * correct patches
   *
   * Return: void
   * Parameters: void
   */
  void initializeECM();

  /*
   * Description:	Initializes all damaged patches. Sprouts platelets and
   * fragments ECM on the damaged patches.
   *
   * Return: void
   * Parameters: void
   */
  // void initializeDamage();

  /*
   * Description:	Each call to go() simulates 30 minutes or 'real-world'
   * time of the biological model.
   *
   * Return: 0 on success
   * Parameters: void
   */
  int go();

  /*
   * Description:	Entry function for sprouting cells
   *              Selects the appropriate sprouting function to apply
   *
   * Return: void
   *
   * Parameters: num          -- Number of cells to sprout
   *             patchType    -- Type of patches of sprout on
   *             agentType    -- Type of agent to sprout
   *
   *             Physical boundaries of the sprouting area/volume:
   *             xmin         -- Left
   *             xmax         -- Right
   *             ymin         -- Top
   *             ymax         -- Bottom
   *             zmin         -- Near
   *             zmax         -- Far
   *
   *             bloodOrTiss  -- Pass in true if should be sprouted in blood
   * 							   Pass in false if
   * should be sprouted in tissue (default)
   */
  void sproutAgent(int num, int patchType, int agentType, int xmin, int xmax,
                   int ymin, int ymax, int zmin, int zmax);

  /*
   * Description:	Function for sprouting cells in a given area/volume
   *
   * Return: void
   *
   * Parameters: num          -- Number of cells to sprout
   *             patchType    -- Type of patches of sprout on
   *             agentType    -- Type of agent to sprout
   *
   *             Physical boundaries of the sprouting area/volume:
   *             xmin         -- Left
   *             xmax         -- Right
   *             ymin         -- Top
   *             ymax         -- Bottom
   *             zmin         --
   *             zmax         --
   *
   *             bloodOrTiss  -- Pass in true if should be sprouted in blood
   *                             Pass in false if should be sprouted in tissue
   * (default)
   */
  void sproutAgentInArea(int num, int patchType, int agentType, int xmin,
                         int xmax, int ymin, int ymax, int zmin, int zmax);

  /*
   * Description:	Function for sprouting cells in the whole world
   *
   * Return: void
   *
   * Parameters: num          -- Number of cells to sprout
   *             patchType    -- Type of patches of sprout on
   *             agentType    -- Type of agent to sprout
   *             bloodOrTiss  -- Pass in true if should be sprouted in blood
   *                             Pass in false if should be sprouted in tissue
   * (default)
   */
  void sproutAgentInWorld(int num, int patchType, int agentType);

  /*
   * Description:	Counts the number of patches of a given type
   *
   * Return: int  -- Number of patches of given patch type
   *
   * Parameters: int  -- Enumic value for the patch type to count
   */
  int countPatchType(int);

  /****************************************************************
   * HELPER SUBROUTINES                                           *
   ****************************************************************/

  /*
   * Description:	Converts a length in millimeters to the number of
   * patches
   *
   * Return: The number of patches
   *
   * Parameters: mm  -- Length in millimeters
   */
  int mmToPatch(double mm);

  /*
   * Description:	Converts hours and days into ticks
   *
   * Return: The number of ticks
   *
   * Parameters: hour  -- Number of hours
   *             day   -- Number of days
   */
  static int reportTick(int hour = 0, int day = 0);

  /*
   * Description:	Determines the number of minutes elapsed
   *
   * Return: The number of minutes elapsed
   *
   * Parameters: void
   */
  static double reportMinute();

  /*
   * Description:	Determines the number of hours elapsed
   *
   * Return: The number of hours elapsed
   *
   * Parameters: void
   */
  static double reportHour();

  /** Chemistry facade (null before constructor finishes chem init). */
  ChemicalEnvironment *chemical_environment();
  const ChemicalEnvironment *chemical_environment() const;

  /** Integrated cytokine masses (from environment when available). */
  float world_total_tnf() const;
  float world_total_tgf() const;
  float world_total_fgf() const;
  float world_total_il6() const;
  float world_total_il8() const;
  float world_total_il10() const;

  float chem_concentration(SpeciesId species, int patch_index) const;
  void chem_add_secretion(SpeciesId species, int patch_index,
                          float delta) const;
  float chem_concentration_channel(int concentration_channel,
                                   int patch_index) const;
  float chemotaxis_at(int patch_index) const;

  /*
   * Description:	Determines the number of days elapsed
   *
   * Return: The number of days elapsed
   *
   * Parameters: void
   */
  static double reportDay();

  /*
   * Description:	Determines the number of neighbors with given patch type
   *
   * Return: The number of neighbors with given patch type
   *
   * Parameters: ix         -- x coordinate of current patch
   *             iy         -- y coordinate of current patch
   *             iz         -- z coordinate of current patch
   *             patchType  -- patch type to count
   */
  int countNeighborPatchType(int ix, int iy, int iz, int patchType);

  /** Load scaffold seeding and alginate composition from simulation_config.json. */
  int userInput();

  // #ifdef MODEL_SCAFFOLD
  /*
   * Description: Calculate Hydrogel initial Elastic modulus, Crosslink density,
   * pore size
   *
   * Return:
   * Parameters:
   */
  void initializeBiomaterial();

  /** Q = (c9 HAww + c10) ln(t_min) + c11 HAww + c12    (Table 4). */
  void updateSwellingRatio();

  /** w_l = (c13 HAww - c14) t_weeks + c15 HAww + c16   (Table 4). */
  void updateMassLoss();

#ifdef PEPTIDE_BM
  /*
   * Description: Update effective stiffness E from peptide stress-relaxation
   * parameters (E_0, E_inf, t).
   */
  void updateE();
#endif

  /*
   * Description:Degrade Ca-Alg Patches and replace with immature tissue patch
   * (no new ECM produced)
   *
   * Return:
   *
   * Parameters:        numOfPatches  -- number of biomaterial patches to "degrade"
   * and be replaced with patch type tissue
   */
  void degradebiomaterial(int numOfPatches);

  /*
   * Description: Print out extra info for debugging purposes
   *
   * Return:
   * Parameters:
   */
  void debugInfo();
  // #endif

  /****************************************************************
   * OUTPUT SUBROUTINES & VISUALIZATION                           *
   ****************************************************************/

  ///*
  // * Description:	Outputs cell counts and cytokine levels from the current
  // tick to the file "Output/Output_Biomarkers.csv".
  // *              Used for testing.
  // *
  // * Return: void
  // * Parameters: void
  // */
  // void outputWorld_csv();

  /*
   * Description:	Outputs all patch assignments (patch type, agent type,
   * ECM type) to files in output directory.
   *
   * Return: void
   * Parameters: void
   */
  void patchassign_csv();

  /****************************************************************
   * WORLD VARIABLES  (Manuscript Table 2)                        *
   ****************************************************************/
  static unsigned seed;

  static float E;         // Elastic modulus, Pa
  static float poreWidth; // p, pore size, um
  static float meshSize;  // mesh = 15000 / p, um^-1
  static float Q;         // Swelling ratio, % w/w
  static float massLoss;  // w_l, mass loss of wet weight, %
  static float pXL;       // rho_XL, crosslinking density, mmol/mL

  static float HAww;  // HA concentration, % w/w of the polymer
  static float HAwv;  // HA concentration, % w/v
  static float TPwv;  // Total HA-Gtn polymer, % w/v
  static float XLww;  // Crosslinker (PEGDA) concentration, % w/w
  static float TDBMR; // Thiol : double bond molar ratio, mol/mol

  static float patchpermm;      // Patches per millimetre
  static int   initialPatches;  // Biomaterial patches present at t = 0
  static float totalVolumeML;   // Construct volume, mL

  static float liveCells, deadCells, deletedCells, prevCells;

  /****************************************************************
  * BIOMATERIAL-RULE PARAMETERS  (Manuscript Table 4, c1 .. c19) *
  ****************************************************************/
  /** E = c1 TPwv HAww + c2 HAww + c3 TPwv + c4 HAwv XLww + c5 XLww + c6 */
  enum ElasticModIdx {
    ELASTIC_POLYMER_HA_INTERACTION = 0,  // c1
    ELASTIC_HA_CONCENTRATION,            // c2
    ELASTIC_POLYMER_CONCENTRATION,       // c3
    ELASTIC_HA_CROSSLINKER_INTERACTION,  // c4
    ELASTIC_CROSSLINKER_CONCENTRATION,   // c5
    ELASTIC_BASELINE,                    // c6
    ELASTIC_MOD_COUNT
  };
  static_assert(sizeof(ElasticModulusParams) / sizeof(double) == ELASTIC_MOD_COUNT,
                "ElasticModulusParams field count must match BMWorld::ElasticModIdx");
  static float ElasticMod[ELASTIC_MOD_COUNT];

  /** rho_XL = c7 - c8 TDB_MR */
  enum XLDensityIdx {
    XLDENSITY_BASELINE = 0,               // c7
    XLDENSITY_THIOL_DOUBLE_BOND_EFFECT,   // c8
    XLDENSITY_COUNT
  };
  static_assert(sizeof(CrosslinkDensityParams) / sizeof(double) == XLDENSITY_COUNT,
                "CrosslinkDensityParams field count must match BMWorld::XLDensityIdx");
  static float XLDensity[XLDENSITY_COUNT];

  /** Q = (c9 HAww + c10) ln(t_m) + c11 HAww + c12 */
  enum SwellRatioIdx {
    SWELL_HA_TIME_EFFECT = 0, // c9
    SWELL_TIME_EFFECT,        // c10
    SWELL_HA_EFFECT,          // c11
    SWELL_BASELINE,           // c12
    SWELL_RATIO_COUNT
  };
  static_assert(sizeof(SwellRatioParams) / sizeof(double) == SWELL_RATIO_COUNT,
                "SwellRatioParams field count must match BMWorld::SwellRatioIdx");
  static float SwellRatio[SWELL_RATIO_COUNT];

  /** w_l = (c13 HAww - c14) t_w + c15 HAww + c16 */
  enum MassLossIdx {
    MASSLOSS_HA_TIME_EFFECT = 0, // c13
    MASSLOSS_TIME_EFFECT,        // c14
    MASSLOSS_HA_EFFECT,          // c15
    MASSLOSS_BASELINE,           // c16
    MASS_LOSS_COUNT
  };
  static_assert(sizeof(MassLossParams) / sizeof(double) == MASS_LOSS_COUNT,
                "MassLossParams field count must match BMWorld::MassLossIdx");
  static float MassLoss[MASS_LOSS_COUNT];

  /** p = -c17 HAww^2 + c18 HAww + c19 */
  enum PoreSizeIdx {
    PORE_HA_QUADRATIC_EFFECT = 0, // c17
    PORE_HA_LINEAR_EFFECT,        // c18
    PORE_BASELINE,                // c19
    PORE_SIZE_COUNT
  };
  static_assert(sizeof(PoreSizeParams) / sizeof(double) == PORE_SIZE_COUNT,
                "PoreSizeParams field count must match BMWorld::PoreSizeIdx");
  static float PoreSize[PORE_SIZE_COUNT];

  /****************************************************************
   * CONSTANT VARIABLES                                           *
   ****************************************************************/
  double patchlength; // The length of each patch

  // float pXL;       // Crosslink Density (mmol/mL)
  // float Q;         // Swelling Ratio (%)
  // float w_l;         // Mass Loss (%)
  // float p; // Pore Size (um)

#ifdef PEPTIDE_BM
  string peptide; // Type of peptide used for biomaterial conjugation
#endif

  int typesOfChem; // The number of different chemicals there are in the world
  vector<float> baselineChem; // Initial amount of each chemical in the world
  ECM *worldECM;              // Pointer to array of ECM
  vector<int>
      initialCells; // Initial amount of each cell type (agent) in the world
  
  /** Day-0 ECM from world_init.initial_ecm (ug per patch, Manuscript Table 2). */
  double initialCollagenPerPatch = 0.0;
  double initialElastinPerPatch  = 0.0;
  double initialHAPerPatch       = 0.0;
  
  ArrayChain<Cell *> cells; // ArrayChain to manage all cell data
  vector<Cell *>
      *localNewCells[MAX_NUM_THREADS]; // Vector of pointers to local lists of
                                       // cell pointers to add to global list
  vector<int> initHAcenters;       // Vector of patches which can be centers for
                                   // sprouting original hyaluronan
  unsigned seeds[MAX_NUM_THREADS]; // Seeds used to generate random numbers for
                                   // each thread

  vector<float> tgfLine; // Vector to store TGF values along an x-face line from
                         // boundary to center of ABM grid
  vector<float> fgfLine; // Vector to store FGF values along an x-face line from
                         // boundary to center of ABM grid
  vector<float> il6Line; // Vector to store IL6 values along an x-face line from
                         // boundary to center of ABM grid
  vector<float> il8Line; // Vector to store IL8 values along an x-face line from
                         // boundary to center of ABM grid
  vector<float> il10Line; // Vector to store IL10 values along an x-face line from
                         // boundary to center of ABM grid
  vector<float> tnfLine; // Vector to store TNF values along an x-face line from
                         // boundary to center of ABM grid

  // boundary to center of ABM grid

  int lineY;
  int lineZ;

  // float Alg_v, Alg_wv; // Volume (mL) and final concentration (% w/v) of Alg in
  //                      // Ca-Alg hydrogel
  float Ca_v, Ca_wv;   // Volume (mL) and final concentration (% w/v) of Ca 3400
  // float highMW_alg, lowMW_alg; // ratio components of high and low MW kDa in the
  //                              // alginate hydrogel

  /** Owns patch grids, registry, and per-tick diffusion. */
  std::unique_ptr<ChemicalEnvironment> chemical_environment_;

  /** Fallback tick length (minutes) if environment not initialized. */
  static constexpr double kTickIntervalMinutes = 30.0;

  double tick_interval_minutes() const;

protected:
  // --- output file hooks ---
  char *get_output_filename() override;
  void write_csv_header(std::ofstream &file) override;
  void write_data_row(std::ofstream &file,
                      std::map<std::string, int> &agent_counts,
                      std::map<std::string, float> &env_counts) override;

  // --- auxiliary output hooks ---
  void write_auxiliary_header() override;
  void write_auxiliary_outputs() override;

  // --- agent counting hooks ---
  std::vector<std::string> get_agent_type_names() override;
  void count_agent_types(std::map<std::string, int> &agent_counts) override;

  // --- agent population hooks ---
  int get_total_agent_count() override;

  // --- environment element counting hooks (e.g. ecm) ---
  std::vector<std::string> get_env_type_names() override;
  void count_env(std::map<std::string, float> &env_counts) override;

private:
  /****************************************************************
   * MAJOR SECTION SUBROUTINES - begin                            *
   ****************************************************************/

  /*
   * Description: (Stage 0)	Sprout agents on various patches.
   *
   * Return: void
   * Parameters: hours  -- Current hour in model execution
   */
  // void seedCells(float hours);

  /*
   * Description: (Stage 1) Diffuse all registry species over one tick.
   *
   * Clears d*, reads p*, writes diffusion increment into d* (cells add
   * secretion later). PDE numerics live in ChemicalEnvironment /
   * diffusion3d_core ? not in this file.
   *
   * Return: void
   * Parameters: void
   */
  void diffuseCytokines();

  /*
   * Description:	Helper function for cell execution. Execute all alive
   * chondrocytes.
   *
   * Return: void
   *
   * Parameters: void
   */
  void inline runCells();

  /*
   * Description:	(Stage 2)	Execute cell functions for all living
   * cells.
   *
   * Return: void
   *
   * Parameters: void
   */
  void executeCells();

  /*
   * Description:	(Stage 3)	Execute ECM functions.
   *
   * Return: void
   *
   * Parameters: void
   */
  void executeECMs();

  /*
   * Description:	(Stage 3)	Fragment ECM proteins if necessary.
   *
   * Return: void
   *
   * Parameters: void
   */
  void requestECMfragments();

  /*
   * Description: (Stage 4a) Commit per-tick chemical changes to patch
   * concentrations.
   *
   * For each diffusing species: p* += d* (diffusion increment + cell
   * secretion), then d* = 0. Called after executeCells() so d* holds both
   * contributions when merged.
   *
   * Return: void
   * Parameters: void
   */
  void updateChem();

  /*
   * Description: Same merge as updateChem(); CPU implementation used each tick.
   *
   * Return: void
   * Parameters: void
   */
  void updateChemCPU();

  /*
   * Description: Update O2 across boundary patches (coming in from external
   * environment). Gets called in updateChem().
   *
   * TO-DO: move this to within updateChem() using a generic
   * "apply external source" option or something
   *
   * Return: void
   * Parameters: void
   */
  // void updateO2();

  /*
   * Description:	Helper function for ECM updates. Execute updates for ALL
   * ECM managers.
   *
   * Return: void
   *
   * Parameters: void
   */
  void inline executeAllECMUpdates();

  /*
   * Description:	Helper function for ECM updates. Execute request resets
   * for ALL ECM managers.
   *
   * Return: void
   *
   * Parameters: void
   */
  void inline executeAllECMResetRequests();

  /*
   * Description:	(Stage 4c)	Update ECM managers to reflect next
   * tick's states
   *
   * Return: void
   *
   * Parameters: void
   */
  void updateECMManagers();

  /*
   * Description:	(Stage 4b)	Update cells to reflect next tick's
   * states
   *
   * Return: void
   *
   * Parameters: void
   */
  void updateCells();

  /*
   * Description:	Update cells to reflect next tick's states.
   *              This is called instead of updateCells() during setup.
   *
   * Return: void
   *
   * Parameters: void
   */
  void updateCellsInitial();

  // --- output function hooks ---
  // viability and differentiation are internal calculations
  // used only in write_data_row ? not exposed as hooks
  float calculate_viability();
  // float calculate_pct_differentiated(std::map<std::string, int> &agent_counts);

  float calculate_mean_displacement();
};

#ifdef PEPTIDE_BM
/* Experimental values for peptide biomaterial */
struct peptideCondition {
  float E_init;   // initial modulus, represents elastic portion of model
  float E_eq;     // 'equilibrium' modulus, represents viscous portion of model
  float t_stress; // stress relaxation time (t = t_half/ln(2))
};
extern peptideCondition MAL, CHAD, hA5G26,
    IKVAV; // names of the peptide structs
// struct peptideCondition MAL, CHAD, hA5G26, IKVAV; // names of the peptide
// structs
#endif // PEPTIDE_BM

#endif /* BMWorld_H */
