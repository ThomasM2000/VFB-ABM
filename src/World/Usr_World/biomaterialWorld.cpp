/*
 * biomaterialWorld.cpp
 *
 * File contents: Contains the BMWorld class.
 *
 * Author: Yvonna
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 */

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <tgmath.h>
#include <time.h>
#include <vector>
#define PI 3.14159
#include "../../Utilities/error_utils.h"
#include "../../Utilities/input_utils.h"
#include "../../Utilities/timer.h"
#include "../../enums.h"
#include "biomaterialWorld.h"
#include "world_init_config.h"
#include <omp.h>
#include <sstream>
#include <string>

using namespace std;

/* ------------------------------------------------------------------------------------
 */
/*                            STATIC VARIABLES INITIALIZATIONS */
/* ------------------------------------------------------------------------------------
 */
unsigned BMWorld::seed = 27000;
float BMWorld::patchpermm = 0;
float BMWorld::liveCells = 0;
float BMWorld::deadCells = 0;
float BMWorld::deletedCells = 0;
float BMWorld::prevCells = 0;

#ifdef MODEL_SCAFFOLD
int   BMWorld::initialPatches = 0;
float BMWorld::totalVolumeML = 0;

/* World variables, Manuscript Table 2 */
float BMWorld::E         = 0;   // Pa
float BMWorld::poreWidth = 0;   // um
float BMWorld::meshSize  = 0;   // um^-1
float BMWorld::Q         = 0;   // % w/w
float BMWorld::massLoss  = 0;   // %
float BMWorld::massLoss0 = 0;   // %
float BMWorld::pXL       = 0;   // mmol/mL

float BMWorld::HAww  = 0;
float BMWorld::HAwv  = 0;
float BMWorld::TPwv  = 0;
float BMWorld::XLww  = 0;
float BMWorld::TDBMR = 0;

/* Biomaterial-rule parameters, Manuscript Table 4 (c1..c19), filled from JSON */
float BMWorld::ElasticMod[BMWorld::ELASTIC_MOD_COUNT] = {};
float BMWorld::XLDensity [BMWorld::XLDENSITY_COUNT]   = {};
float BMWorld::SwellRatio[BMWorld::SWELL_RATIO_COUNT] = {};
float BMWorld::MassLoss  [BMWorld::MASS_LOSS_COUNT]   = {};
float BMWorld::PoreSize  [BMWorld::PORE_SIZE_COUNT]   = {};
#endif

#ifdef PEPTIDE_BM
/* Experimental values for peptide biomaterial */
peptideCondition MAL = {2.0031, 3.8757, 44.00};
peptideCondition CHAD = {3.0180, 5.0709, 37.87};
peptideCondition hA5G26 = {2.9628, 4.6913, 35.35};
peptideCondition IKVAV = {2.9067, 5.4795, 42.56};
#endif

BMWorld::BMWorld(double length, double width, double height, double plength) {
  // Random seeds from --seed so replicate runs differ (Section 2.5).
  srand(util::randomSeed);                    // rand() and random_shuffle()
  for (int i = 0; i < NUM_THREAD; i++)
    seeds[i] = util::randomSeed + 17 * i;     // per-thread rand_r() streams

  // Allocate memory for local lists of cell pointers to add:
  for (int i = 0; i < MAX_NUM_THREADS; i++)
    localNewCells[i] = new vector<Cell *>;

  /* --------------------------------------------------------------------------
   */
  /*                                 GRID SETUP */
  /* --------------------------------------------------------------------------
   */
  this->patchlength = plength;

  // Number of patches in x,y,z dimensions:
  int nx = width / patchlength;
  int ny = length / patchlength;
  int nz = height / patchlength;
  World::setupGrid(nx,     // number of grid points (patches) in x dimension
                   ny,     // number of grid points (patches) in y dimension
                   nz,     // number of grid points (patches) in z dimension
                   0.0,    // min coordinates in x
                   width,  // max corodinates in x
                   0.0,    // min coordinates in y
                   length, // max coordinates in y
                   0.0,    // min coordinates in z
                   height  // max coordinates in z
  );

  // Read input parameters (chem baseline, wound dimensions, initial cells) from
  // config file
    /* Construct volume: needed by userInput(), initializeECM() and
     * initializeCells(), so compute it once as soon as the grid is known. */
  BMWorld::totalVolumeML = (nx * ny * nz) * pow(this->patchlength, 3) * 1e-3;  // mm^3 -> mL

  int temp = BMWorld::userInput();
  cout << "length, width, height: " << length << " mm, " << width << " mm, "
       << height << " mm" << endl;
  cout << "Number of patches: nx, ny, nz: " << nx << ", " << ny << ", " << nz
       << " " << endl;

  // Allocate and initialize Patches/ECM
  if (util::ABMerror(!(this->worldPatch = new Patch[(nx) * (ny) * (nz)]),
                     "Patch mem alloc error!", __FILE__, __LINE__))
    exit(1);
  if (util::ABMerror(!(this->worldECM = new ECM[(nx) * (ny) * (nz)]),
                     "ECM mem alloc error!", __FILE__, __LINE__))
    exit(1);
  cout << "worldPatch size: " << (nx) * (ny) * (nz)
       << " (Number of world patches) " << endl;
  cout << "worldECM size: " << (nx) * (ny) * (nz) << " (Number of ECM patches) "
       << endl;

  /* Try initializing Patches and ECMs with the threads that will access them
   * later since the default allocation policy on Linux platforms is
   * first-touch. This is a best-effort implementation, since we cannot
   * guarantee size of data accessed per thread to be an integer multiple of
   * page size. */
  std::cout << std::fixed;
  std::cout << std::setprecision(3);
  // cout << "	allocating ECM Managers (also Patches) with best-effort first
  // touch" << endl;

  for (int iz = 0; iz < nz; iz++) {
#pragma omp parallel for
    for (int iy = 0; iy < ny; iy++) {
      for (int ix = 0; ix < nx; ix++) {
        int in = ix + iy * nx + iz * nx * ny;
        this->worldPatch[in] = Patch(ix, iy, iz, in);
        this->worldECM[in] = ECM(ix, iy, iz, in);
      }
    }
  }

  // Define Class static variables and pointers
  BMWorld::patchpermm = nx / width;
  Agent::nx = this->nx;
  Agent::ny = this->ny;
  Agent::nz = this->nz;
  Agent::agentWorldPtr = this;
  Agent::agentPatchPtr = this->worldPatch;
  Agent::agentECMPtr = this->worldECM;
  ECM::ECMWorldPtr = this;
  ECM::ECMPatchPtr = this->worldPatch;

  /* ----------------------- INITIALIZATION SUBROUTINES -----------------------
   */

  /* Define initial attributes of patches, damage, ECM, chem and cells based on
   * user defined values (in config file) and traits of native tissue */
  this->initializePatches();
  this->initializeECM();
  this->initializeCells();
#ifdef MODEL_SCAFFOLD
  this->initializeBiomaterial();
#endif
  /* Chemistry controller. */
  {
    const double swelling_Q =
        (this->Q > 0.0f) ? static_cast<double>(this->Q) : 1.0;
    const double stiffness_E =
        (this->E > 0.0f) ? static_cast<double>(this->E) : 1.0;
    chemical_environment_.reset(
        new ChemicalEnvironment(nx, ny, nz, this->patchlength));
    chemical_environment_->load_from_config(
        util::getSimulationConfigPath(), swelling_Q, stiffness_E);
    chemical_environment_->allocate_channels_from_config();
    this->sync_baseline_chem_from_config();
    this->initializeChemBaseline();
  }
  // this->initializeDamage();

  /* Calling update functions to synchronize read and write portion of the
   * attributes */
  // this->updateChemCPU();
  this->updateCellsInitial(); // Add cells to list before removal and updates
  this->updateECMManagers();
  this->updatePatches();
  // cout << "setupGrid complete" << endl;
  cout << "-------------------------------------------" << endl;
}

BMWorld::~BMWorld() {
  for (int i = 0; i < MAX_NUM_THREADS; i++)
    delete localNewCells[i];
  cerr << " removing dead cells" << endl;

  int cellsSize = cells.size();

#pragma omp parallel for
  for (int i = 0; i < cellsSize; i++) {
#ifdef _OMP
    int tid = omp_get_thread_num();
#else
    int tid = DEFAULT_TID;
#endif

    Cell *cell = cells.getDataAt(i);
    if (!cell)
      continue;
    cells.deleteData(i, tid);
    delete cell;
  }

  if (worldPatch != NULL)
    delete[] worldPatch;
  if (worldECM != NULL)
    delete[] worldECM;
  cout << "BMWorld has been successfully destructed." << endl;
}

void destroyCell(Cell *&agent) {
  if (agent) {
    free(agent);
    agent = NULL;
  }
}

void BMWorld::assignPatches(int type, int xmin, int xmax, int ymin, int ymax,
                            int zmin, int zmax) {
  // Assign patches within bounds of type 'type'
  for (int iz = 0; iz < nz; iz++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int ix = 0; ix < nx; ix++) {
        int in = ix + iy * nx + iz * nx * ny;
        switch (type) {
        case damage:
          this->worldPatch[in].type[read_t] = damage;
          this->worldPatch[in].type[write_t] = damage;
          this->worldPatch[in].color[read_t] = cdamage;
          this->worldPatch[in].color[write_t] = cdamage;
          this->worldPatch[in].dirty = true;
          break;
        case biomaterial:
          this->worldPatch[in].type[read_t] = biomaterial;
          this->worldPatch[in].type[write_t] = biomaterial;
          this->worldPatch[in].color[read_t] = cbiomaterial;
          this->worldPatch[in].color[write_t] = cbiomaterial;
          this->worldPatch[in].dirty = true;
          break;
        }
      }
    }
  }
}

void BMWorld::initializePatches() {
#ifdef MODEL_3D
  assignPatches(biomaterial, 0, nx, 0, ny, 0, nz);
#else
  assignPatches(biomaterial, 0, nx, 0, ny, 0, 0);
#endif
  tnfLine.resize(nx / 2 + 1, 0.0f);
  fgfLine.resize(nx / 2 + 1, 0.0f);
  il6Line.resize(nx / 2 + 1, 0.0f);
  il8Line.resize(nx / 2 + 1, 0.0f);
  il10Line.resize(nx / 2 + 1, 0.0f);
  tgfLine.resize(nx / 2 + 1, 0.0f);
  // o2Line.resize(nx / 2 + 1, 0.0f);
  lineY = ny / 2;
  lineZ = nz / 2;

  // Assign values to initial:
  BMWorld::initialPatches = this->countPatchType(biomaterial);
  // cout << "Finished building Ca-Alg Hydrogel" << endl;
}

#ifdef MODEL_SCAFFOLD
/*
 * Manuscript Table 4 - biomaterial rules evaluated once from the composition.
 */
 void BMWorld::initializeBiomaterial() {
  cout << "Computing biomaterial properties (Manuscript Table 4)..." << endl;

  /* rho_XL = c7 - c8 TDB_MR   [mmol/mL] */
  BMWorld::pXL = BMWorld::XLDensity[XLDENSITY_BASELINE]                      // c7
                 - BMWorld::XLDensity[XLDENSITY_THIOL_DOUBLE_BOND_EFFECT]    // c8
                       * BMWorld::TDBMR;

  /* E = c1 TPwv HAww + c2 HAww + c3 TPwv + c4 HAwv XLww + c5 XLww + c6  [Pa] */
  BMWorld::E = BMWorld::ElasticMod[ELASTIC_POLYMER_HA_INTERACTION]           // c1
                   * BMWorld::TPwv * BMWorld::HAww
               + BMWorld::ElasticMod[ELASTIC_HA_CONCENTRATION]               // c2
                     * BMWorld::HAww
               + BMWorld::ElasticMod[ELASTIC_POLYMER_CONCENTRATION]          // c3
                     * BMWorld::TPwv
               + BMWorld::ElasticMod[ELASTIC_HA_CROSSLINKER_INTERACTION]     // c4
                     * BMWorld::HAwv * BMWorld::XLww
               + BMWorld::ElasticMod[ELASTIC_CROSSLINKER_CONCENTRATION]      // c5
                     * BMWorld::XLww
               + BMWorld::ElasticMod[ELASTIC_BASELINE];                      // c6

  /* p = -c17 HAww^2 + c18 HAww + c19   [um] */
  BMWorld::poreWidth =
      -BMWorld::PoreSize[PORE_HA_QUADRATIC_EFFECT] * BMWorld::HAww * BMWorld::HAww // c17
      + BMWorld::PoreSize[PORE_HA_LINEAR_EFFECT] * BMWorld::HAww                   // c18
      + BMWorld::PoreSize[PORE_BASELINE];                                          // c19

  /* mesh = 15000 / p   [um^-1]   (Manuscript Table 2) */
  BMWorld::meshSize =
      (BMWorld::poreWidth > 0.f) ? 15000.f / BMWorld::poreWidth : 0.f;

  /* Q and w_l at t = 0. */
  this->updateSwellingRatio();
  this->updateMassLoss();

  cout << "  Crosslink density rho_XL (mmol/mL) = " << BMWorld::pXL << endl;
  cout << "  Elastic modulus   E      (Pa)      = " << BMWorld::E << endl;
  cout << "  Pore size         p      (um)      = " << BMWorld::poreWidth << endl;
  cout << "  Mesh size         mesh   (um^-1)   = " << BMWorld::meshSize << endl;
  cout << "  Swelling ratio    Q      (% w/w)   = " << BMWorld::Q << endl;
  cout << "  Mass loss         w_l    (%)       = " << BMWorld::massLoss << endl;
}
#endif // MODEL_SCAFFOLD

void BMWorld::sync_baseline_chem_from_config() {
  if (!chemical_environment_) {
    if (util::ABMerror(1, "ChemicalEnvironment not initialized", __FILE__,
                       __LINE__))
      exit(1);
  }

  const ChemicalEnvironmentConfig &cfg = chemical_environment_->configuration();
  this->typesOfChem = cfg.channel_count;

  this->baselineChem.assign(6, 0.f);
  this->baselineChem[TNF] =
      chemical_environment_->baseline_total_mass_for("TNF");
  this->baselineChem[TGF] =
      chemical_environment_->baseline_total_mass_for("TGF");
  this->baselineChem[FGF] =
      chemical_environment_->baseline_total_mass_for("FGF");
  this->baselineChem[IL6] =
      chemical_environment_->baseline_total_mass_for("IL6");
  this->baselineChem[IL8] =
      chemical_environment_->baseline_total_mass_for("IL8");
  this->baselineChem[IL10] =
      chemical_environment_->baseline_total_mass_for("IL10");
  
  cout << "Chemical environment config: " << cfg.model << endl;
  cout << "  tick_interval_minutes = " << cfg.tick_interval_minutes << endl;
  cout << "  channel_count = " << cfg.channel_count << endl;
  cout << "  diffusion algorithm: "
       << chemical_environment_->diffusion_algorithm_label() << endl;
  const DiffusionAlgorithm requested =
      chemical_environment_->diffusion_algorithm();
  if (std::string(chemical_environment_->diffusion_algorithm_label()) !=
      diffusion_algorithm_label(requested)) {
    cout << "  diffusion requested: " << diffusion_algorithm_label(requested)
         << endl;
  }
}

void BMWorld::initializeChemBaseline() {
  if (!chemical_environment_) {
    if (util::ABMerror(1, "ChemicalEnvironment not initialized", __FILE__,
                       __LINE__))
      exit(1);
  }

  chemical_environment_->clear_delta_channels();

  if (this->baselineChem.size() != 6) {
    if (util::ABMerror(1, "Error initializing chemicals!!", __FILE__, __LINE__))
      exit(1);
    return;
  }

  const int countbiomaterial = this->countPatchType(biomaterial);
  const float tnf0 = this->baselineChem[TNF] / countbiomaterial;
  const float tgf0 = this->baselineChem[TGF] / countbiomaterial;
  const float fgf0 = this->baselineChem[FGF] / countbiomaterial;
  const float il60 = this->baselineChem[IL6] / countbiomaterial;
  const float il80 = this->baselineChem[IL8] / countbiomaterial;
  const float il100 = this->baselineChem[IL10] / countbiomaterial;

  const int countBoundary =
      nx * ny * nz -
      (nx - 2) * (ny - 2) *
          (nz - 2); // boundary patches = all patches - interior patches
  const double volumeBoundary =
      (BMWorld::totalVolumeML / (nx * ny * nz)) / 1000 *
      countBoundary; // volume of all boundary patches in L
  // const double molO2 = this->initialO2 * pow(10, 9) *
  //                      volumeBoundary; // total fmol of O2 needed to distribute
  //                                      // across all boundary patches

  // debug
  cout << " Number of boundary patches = " << countBoundary << endl;
  cout << " Total volume of boundary patches = " << volumeBoundary << endl;
  // cout << " Total fmol of O2 for boundary patches = " << molO2 << endl;

  // // TO-DO: move setting of O2 baseline to chem config
  // this->baselineChem[o2] = molO2; // manually set baseline O2

  // const float o20 = this->baselineChem[o2] / countBoundary;

  for (int iz = 0; iz < this->nz; iz++) {
#pragma omp parallel for
    for (int iy = 0; iy < this->ny; iy++) {
      for (int ix = 0; ix < this->nx; ix++) {
        const int in = ix + iy * nx + iz * nx * ny;
        if (this->worldPatch[in].type[read_t] == biomaterial) {
          chemical_environment_->set_concentration(in, TNF, tnf0);
          chemical_environment_->set_concentration(in, TGF, tgf0);
          chemical_environment_->set_concentration(in, FGF, fgf0);
          chemical_environment_->set_concentration(in, IL6, il60);
          chemical_environment_->set_concentration(in, IL8, il80);
          chemical_environment_->set_concentration(in, IL10, il100);
        } else {
          chemical_environment_->set_concentration(in, TNF, 0.f);
          chemical_environment_->set_concentration(in, TGF, 0.f);
          chemical_environment_->set_concentration(in, FGF, 0.f);
          chemical_environment_->set_concentration(in, IL6, 0.f);
          chemical_environment_->set_concentration(in, IL8, 0.f);
          chemical_environment_->set_concentration(in, IL10, 0.f);
        }

        // // set O2 separately only on boundary patches
        // bool isBoundary = (ix == 0 || ix == nx - 1 || iy == 0 || iy == ny - 1 ||
        //                    iz == 0 || iz == nz - 1);
        // if (isBoundary) {
        //   chemical_environment_->set_concentration(in, o2, o20);
        // } else {
        //   chemical_environment_->set_concentration(in, o2, 0.f);
        // }
      }
    }
  }

  chemical_environment_->update_chemotaxis_from_species(TGF);
  chemical_environment_->recompute_world_totals();

  cout << "		Initial cytokine concentrations: totalTNF = "
       << this->world_total_tnf() << ", totalTGF = " << this->world_total_tgf()
       << ", totalFGF = " << this->world_total_fgf()
       << ", totalIL6 = " << this->world_total_il6()
       << ", totalIL8 = " << this->world_total_il8()
       << ", totalIL10 = " << this->world_total_il10()
       << endl;
}

void BMWorld::initializeCells() {
  cells = ArrayChain<Cell *>(DEFAULT_DATA_SMALL, 4, NULL, NULL);

  /* Manuscript Table 1: 50 000 fibroblast agents. A count of 0 falls back to
   * the in vitro seeding density of 1e6 cells/mL (Section 2.2.3). */
  int seedCount = this->initialCells[0];
  if (seedCount <= 0)
    seedCount = static_cast<int>(1.0e6 * BMWorld::totalVolumeML);

  const int available = this->countPatchType(biomaterial);
  if (seedCount > available) seedCount = available;
  this->initialCells[0] = seedCount;

  cout << "Construct volume: " << BMWorld::totalVolumeML << " mL" << endl;
  cout << "Seeding " << seedCount << " fibroblasts ("
       << (seedCount / BMWorld::totalVolumeML) << " cells/mL)" << endl;

  sproutAgent(seedCount, biomaterial, fibroblast, 0, nx, 0, ny, 0, nz);
  BMWorld::prevCells = seedCount;
}


void BMWorld::initializeECM() {
  /* Manuscript Table 2: Col, Eln, HA and HA_f are per-patch state variables.
   * Day-0 values are a measured calibration initial condition (Section 2.2.5),
   * not something the model generates. */
  const int totalPatches = nx * ny * nz;

  const float col0 = static_cast<float>(this->initialCollagenPerPatch);
  const float eln0 = static_cast<float>(this->initialElastinPerPatch);
  const float ha0  = static_cast<float>(this->initialHAPerPatch);

  for (int in = 0; in < totalPatches; in++) {
    this->worldECM[in].ncollagen[read_t]  = col0;
    this->worldECM[in].ncollagen[write_t] = col0;
    this->worldECM[in].nelastin[read_t]   = eln0;
    this->worldECM[in].nelastin[write_t]  = eln0;
    this->worldECM[in].HA[read_t]         = ha0;
    this->worldECM[in].HA[write_t]        = ha0;
    this->worldECM[in].fHA[read_t]        = 0.f;
    this->worldECM[in].fHA[write_t]       = 0.f;
    this->worldECM[in].isEmpty();
  }

  cout << "Initial ECM per patch (pg): collagen " << col0 << ", elastin " << eln0
       << ", HA " << ha0 << "  -> construct totals (pg): "
       << static_cast<double>(col0) * totalPatches << ", "
       << static_cast<double>(eln0) * totalPatches << ", "
       << static_cast<double>(ha0) * totalPatches << endl;
}

int BMWorld::go() {
  cout << "-------------------------------------------" << endl;

  Patch *tempPatchPtr;
  Agent *tempAgentPtr;
  double hours = this->reportHour();
  double days = this->reportDay();

// Profiling options defined in common.h
#ifdef PROFILE_MAJOR_STEPS
  struct timeval start, end;
  long elapsed_time; // in milliseconds
#endif

  // Increment Clock in ticks (1 tick = 30 min)
  BMWorld::clock++;
  cout << "tick: " << clock << " , hour: " << hours << " , day: " << days
       << endl;

#ifdef PROFILE_MAJOR_STEPS
  /* TIME_STAGE is a macro for timing a command/function and printing the timing
   * info (See Utilities/time.h) */

  /* --------------------------- CHEMICAL DIFFUSION ---------------------------
   */
  TIME_STAGE(this->diffuseCytokines(), "Chemical diffusion", "1");

  /* ------------------------------ CELL SEEDING ------------------------------
   */
  // TIME_STAGE(this->seedCells(hours), "Cell seeding", "0");

  /* ------------------------------ CELL FUNCTION -----------------------------
   */
  TIME_STAGE(this->executeCells(), "Cells function", "2");

  /* ------------------------------ ECM FUNCTION ------------------------------
   */
  TIME_STAGE(this->executeECMs(), "ECM function", "3(a)");
  TIME_STAGE(this->requestECMfragments(), "ECM fragmentation",
             "3(b)"); // Request fragment<ECM>

  /* ----------------------- ATTRIBUTES SYNCHRONIZATION -----------------------
   */
  cerr << " begin update... " << endl;
  TIME_STAGE(this->updateCells(), "Update cells", "4(b)");
  TIME_STAGE(this->updateECMManagers(), "Update ECM Managers", "4(c)");
  TIME_STAGE(this->updatePatches(), "Update Patches", "4(d)");
  TIME_STAGE(this->updateChem(), "Update chem", "4(a)");
#else // PROFILE_MAJOR_STEPS

  // For testing purposes:
  // cout << "  TNF: " << this->WorldChem.totalTNF << ", TGF: " <<
  // this->WorldChem.totalTGF << ", IL1beta: " << this->WorldChem.totalIL1beta
  // << endl;

  /* --------------------------- CHEMICAL DIFFUSION ---------------------------
   */
  this->diffuseCytokines();

  /* ------------------------------ CELL SEEDING ------------------------------
   */
  // this->seedCells(hours);

  /* ------------------------------ CELL FUNCTION -----------------------------
   */
  this->executeCells();

  /* ------------------------------ ECM FUNCTION ------------------------------
   */
  this->executeECMs();
  this->requestECMfragments();

#ifdef MODEL_SCAFFOLD
  /* ------------------------- UPDATE biomaterial Properties ------------------------
   */
  this->updateSwellingRatio();
  this->updateMassLoss();
#ifdef PEPTIDE_BM
  this->updateE();
#endif // PEPTIDE_BM
#endif // MODEL_SCAFFOLD

  /* ----------------------- ATTRIBUTES SYNCHRONIZATION -----------------------
   */
  cerr << " begin update... " << endl;
  this->updateCells();
  this->updateECMManagers();
  this->updatePatches();
  this->updateChem();

#endif
  return 0;
}

/* -------------------------------------------------------------------------- */
/*                      MAJOR SECTION SUBROUTINES - begin                     */
/* -------------------------------------------------------------------------- */
/*
 * Chemical diffusion tick - delegated to ChemicalEnvironment (see go()
 * ordering). PDE numerics run inside ChemicalEnvironment (Diffusion3D).
 */

ChemicalEnvironment *BMWorld::chemical_environment() {
  return chemical_environment_.get();
}

const ChemicalEnvironment *BMWorld::chemical_environment() const {
  return chemical_environment_.get();
}

float BMWorld::world_total_tnf() const {
  if (chemical_environment_)
    return chemical_environment_->total_tnf();
  return 0.f;
}

float BMWorld::world_total_tgf() const {
  if (chemical_environment_)
    return chemical_environment_->total_tgf();
  return 0.f;
}

float BMWorld::world_total_fgf() const {
  return chemical_environment_ ? chemical_environment_->total_fgf() : 0.f;
}
float BMWorld::world_total_il6() const {
  return chemical_environment_ ? chemical_environment_->total_il6() : 0.f;
}
float BMWorld::world_total_il8() const {
  return chemical_environment_ ? chemical_environment_->total_il8() : 0.f;
}
float BMWorld::world_total_il10() const {
  return chemical_environment_ ? chemical_environment_->total_il10() : 0.f;
}

double BMWorld::tick_interval_minutes() const {
  if (chemical_environment_)
    return chemical_environment_->tick_interval_minutes();
  return kTickIntervalMinutes;
}

float BMWorld::chem_concentration(SpeciesId species, int patch_index) const {
  if (!chemical_environment_)
    return 0.f;
  return chemical_environment_->concentration_at(patch_index, species);
}

void BMWorld::chem_add_secretion(SpeciesId species, int patch_index,
                                 float delta) const {
  if (chemical_environment_)
    chemical_environment_->accumulate_secretion(patch_index, species, delta);
}

float BMWorld::chem_concentration_channel(int concentration_channel,
                                          int patch_index) const {
  if (!chemical_environment_)
    return 0.f;
  return chemical_environment_->concentration_at_channel(patch_index,
                                                         concentration_channel);
}

float BMWorld::chemotaxis_at(int patch_index) const {
  if (!chemical_environment_)
    return 0.f;
  return chemical_environment_->chemotaxis_at(patch_index);
}

void BMWorld::diffuseCytokines() {
  if (chemical_environment_)
    chemical_environment_->run_diffusion_phase(tick_interval_minutes());
}

void BMWorld::runCells() {
  int cellsSize =
      cells.size(); /* This is only an upper bound on cell list size. It is NOT
                       an actual count of cells (some entries are NULL) */
#pragma omp parallel for
  for (int i = 0; i < cellsSize; i++) {
    Cell *cell = cells.getDataAt(i);
    if (!cell)
      continue;
    cell->cellFunction();
  }
}

void BMWorld::executeCells() {
#ifdef PROFILE_CELL_FUNC
  TIME_STAGE(this->runCells(), "Cell Function: Chondrocytes", "	");
#else
  cout << " execute cells " << endl;
  this->runCells();
#endif
}

void BMWorld::executeECMs() {
  cerr << " ECM function  " << endl;
  const int numPatches = nx * ny * nz;
#pragma omp parallel for
  for (int in = 0; in < numPatches; in++) {
    if (worldECM[in].empty[read_t] == false)
      this->worldECM[in].ECMFunction();
  }
}

void BMWorld::requestECMfragments() {
//   if (BMWorld::highTNFdamage == true) {
//     cout << " high TNF damage " << endl;
//     BMWorld::highTNFdamage = false;
//     for (int in = 0; in < (nx - 1) + (ny - 1) * nx + (nz - 1) * nx * ny; in++) {
// #ifndef CALIBRATION
//       if (this->chem_concentration(TNF, in) > 10.f) {
// #else
//       if (this->chem_concentration(TNF, in) > 10.f) {
// #endif
//         cout << " Degrade ECM " << endl;
//         this->worldECM[in].fragmentNCollagen();
//       }
//     }
//   }
}

// void BMWorld::updateO2() {
//   // TO-DO: move this O2 boundary update to somewhere within the chem
//   // environment, not directly in BMWorld
//   const int countBoundary =
//       nx * ny * nz -
//       (nx - 2) * (ny - 2) *
//           (nz - 2); // boundary patches = all patches - interior patches
//   const double volumeBoundary =
//       (BMWorld::totalVolumeML / (nx * ny * nz)) / 1000 *
//       countBoundary; // volume of all boundary patches in L
//   const double molO2 = this->initialO2 * pow(10, 9) *
//                        volumeBoundary; // total fmol of O2 needed to distribute
//                                        // across all boundary patches
//   const double incrO2 =
//       this->incrementO2 *
//       (molO2 / countBoundary); // oxygen increment for each patch

//   // update boundary O2 across Z faces
//   for (int zi : {0, nz - 1}) {
// #pragma omp parallel for
//     for (int yi = 0; yi < ny; yi++) {
// #pragma omp simd
//       for (int xi = 0; xi < nx; xi++) {
//         int in = xi + yi * nx + zi * nx * ny;
//         chemical_environment_->accumulate_secretion(in, o2, incrO2);
//       }
//     }
//   }

//   // update boundary O2 across Y faces
//   for (int yi : {0, ny - 1}) {
// #pragma omp parallel for
//     for (int zi = 1; zi < nz - 1; zi++) {
//       // zi starts at 1 and ends at nz-2 to exclude patches
//       // already counted by the Z face loops above
// #pragma omp simd
//       for (int xi = 0; xi < nx; xi++) {
//         int in = xi + yi * nx + zi * nx * ny;
//         chemical_environment_->accumulate_secretion(in, o2, incrO2);
//       }
//     }
//   }

//   // update boundary O2 across X faces
//   for (int xi : {0, nx - 1}) {
// #pragma omp parallel for
//     for (int zi = 1; zi < nz - 1; zi++) {
//       for (int yi = 1; yi < ny - 1; yi++) {
//         int in = xi + yi * nx + zi * nx * ny;
//         chemical_environment_->accumulate_secretion(in, o2, incrO2);
//       }
//     }
//   }
// }

void BMWorld::updateChemCPU() {
  if (!chemical_environment_)
    return;
  chemical_environment_->merge_and_reset_secretion();

  // temp for TGF output for diffusion debugging
  for (int xi = 0; xi < nx / 2; xi++) {
    int in = xi + lineY * nx + lineZ * nx * ny;
    tgfLine[xi] = chemical_environment_->concentration_at(in, TGF);
  }

  // // temp for O2 output for diffusion debugging
  // for (int xi = 0; xi < nx / 2; xi++) {
  //   int in = xi + lineY * nx + lineZ * nx * ny;
  //   o2Line[xi] = chemical_environment_->concentration_at(in, o2);
  // }
}

void BMWorld::updateChem() { updateChemCPU(); }

void BMWorld::executeAllECMUpdates() {
  for (int iz = 0; iz < nz; iz++) {
#pragma omp parallel for
    for (int iy = 0; iy < ny; iy++) {
      for (int ix = 0; ix < nx; ix++) {
        int in = ix + iy * nx + iz * nx * ny;
        this->worldECM[in].updateECM();
      }
    }
  }
}

void BMWorld::executeAllECMResetRequests() {
  for (int iz = 0; iz < nz; iz++) {
#pragma omp parallel for
    for (int iy = 0; iy < ny; iy++) {
      for (int ix = 0; ix < nx; ix++) {
        int in = ix + iy * nx + iz * nx * ny;
        this->worldECM[in].resetrequests();
      }
    }
  }
}

void BMWorld::updateECMManagers() {
#ifdef PROFILE_ECM_UPDATE
  TIME_STAGE(this->executeAllECMUpdates(), "	updateECM()", "	");
  TIME_STAGE(this->executeAllECMResetRequests(), "	resetrequests()",
             "	");
#else
  this->executeAllECMUpdates();
  this->executeAllECMResetRequests();
#endif
}

/*
 * Steps:
 * 1. Perform updates
 * 2. Remove all dead cells
 * 3. If OMP, add cells from thread-local lists to corresponding global lists
 */
void BMWorld::updateCells() {
  cout << " previous tick's cells " << prevCells << endl;
  cerr << "	removing dead cells" << endl;
  int cellsSize = cells.size();
  liveCells = 0;
  int dcells = 0;
  deletedCells = 0;
#pragma omp parallel for
  for (int i = 0; i < cellsSize; i++) {
#ifdef _OMP
    int tid = omp_get_thread_num();
#else
    int tid = DEFAULT_TID;
#endif
    // Get pointer of cell i from the array chain
    Cell *cell = cells.getDataAt(i);
    if (!cell)
      continue; // cell was deleted
    cell->updateAgent();

    if (cell->isRealDead() == true) {
      dcells++;
    }

    if (cell->isAlive() == true) {
      liveCells++;
      // Update cell stage counts
      /* Added by MM to check types of cell stages and add to respective
       * counters: */
      if (typeid(*cell) == typeid(Fibroblast)) {
        Fibroblast::numOfFibroblast++;
      } 

      // Remove dead cells
    } else if (cell->isAlive() == false) {
      // Get residing patch index and update its occupancy
      // int in = cell->getIndex();
      // this->worldPatch[in].clearOccupied();
      // this->worldPatch[in].occupiedby[write_t] = nothing;
      // this->worldPatch[in].dirty = true;
      /* Added by MM to check types of cell stages and subtract from respective
       * counters: */
      if (typeid(*cell) == typeid(Fibroblast)) {
        Fibroblast::numOfFibroblast--;
      } 
      cells.deleteData(i, tid);
      delete cell;
      deletedCells++;
    }
  }
  // if (BMWorld::clock == 0) {
  //	prevCells = this->initialCells[0];
  // }
  // else if (BMWorld::clock > 0) {
  //	prevCells = liveCells;
  // }

  deadCells += dcells;
  // if (prevCells - cells.actualSize() >= 0) {
  //   deadCells += prevCells - cells.actualSize();
  // } else if (prevCells - cells.actualSize() < 0) {
  //   deadCells += 0;
  // }
  // deadCells += prevCells - cells.actualSize();
  Cell::numOfCells = cells.actualSize();
  // cout << " number of dead cells in this tick " << prevCells -
  // cells.actualSize() << endl;
  cout << " number of dead cells in this tick (prev - cells actual size) "
       << prevCells - cells.actualSize() << endl;
  cout << " number of dead cells in this tick (dcells) " << dcells << endl;
  cout << " number of cells now (actualSize) = " << cells.actualSize() << endl;
  cout << " number of live cells = " << liveCells << endl;
  cout << " total number of dead cells " << deadCells << endl;
  cout << " cell viability " << std::setprecision(3)
       << (liveCells / (liveCells + deadCells)) * 100 << "%" << endl;

// Add new cells
#ifdef _OMP
  /* In OMP, cells were only added to each thread's local list when
   * sproutAgent() was called. Thus, this step is needed to add those cells onto
   * the global lists. */
  cerr << "	updateCell() _OMP" << endl;
  // TODO: parallelize
  // int numThreads = omp_get_num_threads();
  int numThreads = std::max(atoi(std::getenv("OMP_NUM_THREADS")), 1);
  // cout << "		numThreads = " << numThreads << endl;
  for (int tid = 0; tid < numThreads; tid++) {
    // Cells
    vector<Cell *> *fvec_ptr = localNewCells[tid];
    for (vector<Cell *>::iterator cell_it = fvec_ptr->begin();
         cell_it != fvec_ptr->end(); cell_it++) {
      Cell *newCell = *cell_it;
      if (!cells.addData(newCell, tid)) {
        cerr << "Error: Could not add cell" << endl;
        exit(-1);
      }
    }
    fvec_ptr->clear();
  }
#endif
  if (liveCells > nx * ny * nz) {
    cerr << "ERROR: " << liveCells << " live cells exceed " << nx * ny * nz
         << " patches (agents stacked on one patch)" << endl;
  }
  prevCells = cells.actualSize();
}

/****************************************************************
 * MAJOR SECTION SUBROUTINES - end                              *
 ****************************************************************/
// NOTE: only use this function to sprout new_coll, new_agg in initialization.

void BMWorld::sproutAgent(int num, int patchType, int agentType, int xmin,
                          int xmax, int ymin, int ymax, int zmin, int zmax) {
#ifdef OPT_CELL_SEEDING
  if (xmin != 0 || xmax != nx || ymin != 0 || ymax != ny || zmin != 0 ||
      zmax != nz)
    sproutAgentInArea(num, patchType, agentType, xmin, xmax, ymin, ymax, zmin,
                      zmax);
  else if (patchType == biomaterial)
    sproutAgentInArea(num, patchType, agentType, xmin, xmax, ymin, ymax, zmin,
                      zmax);
  else
    sproutAgentInWorld(num, patchType, agentType);
#else
  // Target a specific area of the world
  sproutAgentInArea(num, patchType, agentType, xmin, xmax, ymin, ymax, zmin,
                    zmax);
#endif
}

void BMWorld::sproutAgentInArea(int num, int patchType, int agentType, int xmin,
                                int xmax, int ymin, int ymax, int zmin,
                                int zmax) {
  int count = 0;
  vector<int> patchlist;
  int *reservoir = new int[num];
  for (int i = 0; i < num; i++)
    reservoir[i] = -1;
  Patch *tempPatchPtr;
  Agent *tempAgentPtr;
  int in, agentIndex, max;
  for (int izz = zmin; izz < zmax + 1; izz++) {
    for (int iyy = ymin; iyy < ymax + 1; iyy++) {
      for (int ixx = xmin; ixx < xmax + 1; ixx++) {
        in = ixx + iyy * nx + izz * nx * ny;

        // Try another patch if this one is out of bounds or the wrong type or
        // occupied
        if (ixx < 0 || ixx >= nx || iyy < 0 || iyy >= ny || izz < 0 ||
            izz >= nz)
          continue;
        if (BMWorld::worldPatch[in].type[read_t] != patchType)
          continue;
        if (this->worldPatch[in].isOccupied() == false)
          patchlist.push_back(in);
      }
    }
  }
  /* Sample without replacement so no two agents share a patch. */
  std::random_shuffle(patchlist.begin(), patchlist.end());
  if (static_cast<int>(patchlist.size()) < num) {
    cout << " sprout agent warning: only " << patchlist.size()
         << " free patches for " << num << " agents" << endl;
    num = static_cast<int>(patchlist.size());
  }
  for (int i = 0; i < num; i++)
    reservoir[i] = patchlist[i];


  // Sprout agent on each patch in reservoir
  for (int i = 0; i < num; i++) {
    int in = reservoir[i];
    if (in < 0 || in > (nx - 1) + (ny - 1) * nx + (nz - 1) * nx * ny)
      continue;
    switch (agentType) {
    case fibroblast: {
      // cout << "new fibroblast added" << endl; //added for debugging
      tempPatchPtr = &(this->worldPatch[in]);
      // std::shared_ptr<Cell> cell = std::make_shared<Fibroblast>(tempPatchPtr);
      Fibroblast *newFibroblast = new Fibroblast(tempPatchPtr);
#ifdef _OMP
      int tid = omp_get_thread_num();
      this->localNewCells[tid]->push_back(newFibroblast);
#else
      if (!this->cells.addData(newFibroblast, DEFAULT_TID)) {
        cerr << "Error: Could not add fibroblast cell in sproutAgent()" << endl;
        exit(-1);
      }
#endif
      ///* Added by MM to check types of cell stages and add to respective
      /// counters: */
      // if (typeid(*this) == typeid(Fibroblast)) {
      //	Fibroblast::numOfFibroblast++;
      // }
      // Cell::numOfCells++;

      this->worldPatch[in].setOccupied();
      this->worldPatch[in].occupiedby[write_t] = fibroblast;
      this->worldPatch[in].dirty = true;
      // cout << "patch index " << in << endl; //added for debugging

      break;
    }
    
    // case orig_coll: {
    //   this->worldECM[in].ocollagen[write_t] =
    //       this->worldECM[in].ocollagen[read_t] + 1;
    //   this->worldECM[in].dirty = true;
    //   this->worldECM[in].isEmpty();
    //   break;
    // }
    case new_coll: {
      this->worldPatch[in].initcollagen = true;
      this->worldECM[in].ncollagen[write_t] =
          this->worldECM[in].ncollagen[read_t] + 1;
      this->worldECM[in].dirty = true;
      this->worldECM[in].isEmpty();
      break;
    }
    }
  }
  delete[] reservoir;
  return;
}

#ifdef OPT_CELL_SEEDING
/*
 * Optimized by:
 *  - If sprout in tissue:
 *     (*) Randomly choosing a target patch:
 *          - If is tissue and occupied, sprout
 * 			    - Else repeat (*)
 *  - Else (sprout in blood):
 *     (**)	Look at the list of capillary patches initialized in the setup
 * stage
 *          - Create a list of unoccupied capillary patches
 *          - Pick randomly and sprout
 *          - Repeat until 'num' cells are sprouted
 */
void BMWorld::sproutAgentInWorld(int num, int patchType,
                                 int agentType) { // bool bloodORtiss
  int count = 0;
  vector<int> patchlist;
  int *reservoir = new int[num];
  for (int i = 0; i < num; i++)
    reservoir[i] = -1;
  Patch *tempPatchPtr;
  Agent *tempAgentPtr;
  int in, agentIndex, max;
  int totalNumPatches =
      this->(n - 1) x * this->(ny - 1) * this->(nz - 1); // Is this accurate?
  int numfound;
  int counter;
  int threshold;
  vector<int> unoccupiedCaps;
  switch (patchType) {
    for (int i = 0; i < num; i++) {
      if (reservoir[i] < 0 ||
          reservoir[i] > (nx - 1) + (ny - 1) * nx + (nz - 1) * nx * ny)
        continue;
      switch (agentType) {
      case fibroblast: {
        tempPatchPtr = &(this->worldPatch[reservoir[i]]);
        Fibroblast *newFibroblast = new Fibroblast(tempPatchPtr);
#ifdef _OMP
        int tid = omp_get_thread_num();
        this->localNewCells[tid]->push_back(newFibroblast);
#else
        if (!this->cells.addData(newFibroblast, DEFAULT_TID)) {
          cerr << "Error: Could not add fibroblast cell in sproutAgentInWorld()"
               << endl;
          exit(-1);
        }
#endif
        // this->worldPatch[reservoir[i]].occupied = true;
        this->worldPatch[reservoir[i]].setOccupied();
        this->worldPatch[reservoir[i]].occupiedby = fibroblast;
        break;
      }
      
      }
    }
    delete[] reservoir;
    return;
  }
#endif // OPT_CELL_SEEDING

  int BMWorld::countPatchType(int whichType) {
    if (whichType == biomaterial) {
      Patch::numOfEachTypes[whichType] = 0;
      for (int iz = 0; iz < this->nz; iz++) {
        int currCount = 0;
#pragma omp parallel for reduction(+ : currCount)

        for (int iy = 0; iy < this->ny; iy++) {
          for (int ix = 0; ix < this->nx; ix++) {
            int in = ix + iy * nx + iz * nx * ny;
            if (this->worldPatch[in].type[read_t] == whichType) {
              currCount++;
              // Patch::numOfEachTypes [whichType]++;
            }
          }
        }
        Patch::numOfEachTypes[whichType] += currCount;
      }
    } else if (whichType == damage) {
      Patch::numOfEachTypes[whichType] = 0;
      for (int iz = 0; iz < this->nz; iz++) {
        int currCount = 0;
#pragma omp parallel for reduction(+ : currCount)
        for (int iy = 0; iy < this->ny; iy++) {
          for (int ix = 0; ix < this->nx; ix++) {
            int in = ix + iy * nx + iz * nx * ny;
            currCount += this->worldPatch[in].damage[read_t];
          }
        }
        Patch::numOfEachTypes[whichType] += currCount;
      }
    } else
      cout << "type must be 0, 1, 2 , 3 or 4!"
           << endl; // cout << "type must be 0, 1, 2 , 3, 4 or 5!" << endl;
    return Patch::numOfEachTypes[whichType];
  }

  int BMWorld::mmToPatch(double mm) { return mm * (this->patchpermm); }

  int BMWorld::reportTick(int hour, int day) { return (hour * 2 + day * 48); }

  double BMWorld::reportMinute() { return (BMWorld::clock) * 30; }

  double BMWorld::reportHour() { return (BMWorld::clock) / 2; }

  double BMWorld::reportDay() { return (BMWorld::clock) / 48; }

  int BMWorld::countNeighborPatchType(int ix, int iy, int iz, int patchType) {
    int neighborcount = 0;
    for (int dx = -1; dx < 2; dx++) {
      for (int dy = -1; dy < 2; dy++) {
        for (int dz = -1; dz < 2; dz++) {
          int neighborindex = (ix + dx) + (iy + dy) * nx + (iz + dz) * ny * nx;
          if (ix + dx < 0 || ix + dx >= nx || iy + dy < 0 || iy + dy >= ny ||
              iz + dz < 0 || iz + dz >= nz)
            continue;
          if (dx == 0 && dy == 0 && dz == 0)
            continue;
          if (Agent::agentPatchPtr[neighborindex].type[read_t] == patchType)
            neighborcount++;
        }
      }
    }
    return neighborcount;
  }

#ifdef MODEL_SCAFFOLD
  /* Table 4: Q = (c9 HAww + c10) ln(t_m) + c11 HAww + c12   [% w/w] */
  void BMWorld::updateSwellingRatio() {
    double tmin = reportMinute();
    if (tmin < 1.0) tmin = 1.0;               // ln(t) domain guard

    BMWorld::Q = (BMWorld::SwellRatio[SWELL_HA_TIME_EFFECT] * BMWorld::HAww   // c9
                  + BMWorld::SwellRatio[SWELL_TIME_EFFECT])                   // c10
                     * static_cast<float>(log(tmin))
                 + BMWorld::SwellRatio[SWELL_HA_EFFECT] * BMWorld::HAww       // c11
                 + BMWorld::SwellRatio[SWELL_BASELINE];                       // c12

    /* Swelling reduces effective diffusivity: D_eff = base_D * Q. */
    if (BMWorld::Q > 0.0f && chemical_environment_)
      chemical_environment_->set_swelling_ratio(static_cast<double>(BMWorld::Q));
    if (BMWorld::E > 0.0f && chemical_environment_)
      chemical_environment_->set_stiffness(static_cast<double>(BMWorld::E));
  }

  /* Table 4: w_l = (c13 HAww - c14) t_w + c15 HAww + c16   [%] */
  void BMWorld::updateMassLoss() {
    const float tweek = static_cast<float>(reportDay()) / 7.f;

    float w_t = (BMWorld::MassLoss[MASSLOSS_HA_TIME_EFFECT] * BMWorld::HAww    // c13
                 - BMWorld::MassLoss[MASSLOSS_TIME_EFFECT])                    // c14
                    * tweek
                + BMWorld::MassLoss[MASSLOSS_HA_EFFECT] * BMWorld::HAww        // c15
                + BMWorld::MassLoss[MASSLOSS_BASELINE];                        // c16

    if (w_t < 0.f) w_t = 0.f;
        /* The fit has a non-zero intercept (w_l(0) = c15 HAww + c16 = 65.3 %). That
     * is a property of the regression, not mass the construct has already lost.
     * Without this guard the t = 0 call degrades 65 % of the patches before the
     * first tick; every cell on a degraded patch then fails the
     * `type != biomaterial` test in proliferate() and ecm_synthesis() for the
     * rest of the run. */
    if (BMWorld::clock == 0) {
      BMWorld::massLoss = w_t;
      BMWorld::massLoss0 = w_t;
      return;
    }

    /* Degrade against the cumulative target so the per-tick float->int
     * truncation does not accumulate. massLoss0 is w_l at t = 0 (see guard). */
     const int target = static_cast<int>(std::lround(
      0.01 * (w_t - BMWorld::massLoss0) * BMWorld::initialPatches));
     const int already = BMWorld::initialPatches - this->countPatchType(biomaterial);
     if (target > already)
       this->degradebiomaterial(target - already);
     BMWorld::massLoss = w_t;

  }
// #ifdef PEPTIDE_BM
//   void BMWorld::updateE() {
//     BMWorld::E =
//         BMWorld::E_inf + (BMWorld::E_0 - BMWorld::E_inf) *
//                              exp(-(BMWorld::clock * 30 * 60) /
//                                  BMWorld::t); // converts tick to seconds
//     if (BMWorld::E > 0.0f && chemical_environment_)
//       chemical_environment_->set_stiffness(static_cast<double>(BMWorld::E));
//   }
// #endif // PEPTIDE_BM
#endif // MODEL_SCAFFOLD

#ifdef MODEL_SCAFFOLD
  void BMWorld::degradebiomaterial(int numOfPatches) {
    int xmin = 0;
    int xmax = nx;
    int ymin = 0;
    int ymax = ny;
    int zmin = 0;
    int zmax = nz;
    int patchType = biomaterial;
    vector<int> patchlist;
    // int *reservoir = new int[numOfPatches];
    // for (int i = 0; i < numOfPatches; i++)
    //   reservoir[i] = -1;
    int in, agentIndex, max;
    int count = 0;

    // Make list of possible biomaterial Patches to degrade
    for (int iz = zmin; iz < zmax; iz++) {
      for (int iy = ymin; iy < ymax; iy++) {
        for (int ix = xmin; ix < xmax; ix++) {
          in = ix + iy * nx + iz * nx * ny;

          // Try another patch if this one is out of bounds or the wrong type or
          // occupied
          if (ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz)
            continue;
          if (BMWorld::worldPatch[in].type[read_t] != patchType)
            continue;
          patchlist.push_back(in);
        }
      }
    }

    // // Choose random patches from patch list
    // for (int i = 0; i < numOfPatches; i++) {
    //   if (patchlist.size() == 0) { // No available patches
    //     // cout << " biomaterial degrade error, no available patch within bounds! " <<
    //     // endl;
    //     delete[] reservoir;
    //     return;
    //   }
    //   int randnumber = rand() % patchlist.size();
    //   reservoir[i] = patchlist[randnumber]; // Prepare 'num' random patches
    // }

    // // Degrade 'numOfPatches" number of biomaterial patches in reservoir list
    // for (int i = 0; i < numOfPatches; i++) {
    //   int in = reservoir[i];
    //   if (in < 0 || in > (nx - 1) + (ny - 1) * nx + (nz - 1) * nx * ny)
    //     continue;
    //   BMWorld::worldPatch[in].type[write_t] = nothing;
    //   BMWorld::worldPatch[in].color[write_t] = cnothing;
    //   BMWorld::worldPatch[in].dirty = true;
    //   count++;
    // }
    // delete[] reservoir;
    if (numOfPatches <= 0 || patchlist.empty()) return;
    std::random_shuffle(patchlist.begin(), patchlist.end());
    const int n = std::min<int>(numOfPatches, static_cast<int>(patchlist.size()));
    for (int i = 0; i < n; i++) {
      const int in = patchlist[i];
      BMWorld::worldPatch[in].type[write_t]  = nothing;
      BMWorld::worldPatch[in].color[write_t] = cnothing;
      BMWorld::worldPatch[in].dirty = true;
    }

  }
#endif // MODEL_SCAFFOLD

  void BMWorld::debugInfo() {
    int alive = 0;
    int dead = 0;
    int fibroblastSize = 0;

    int cellsSize = cells.size();
    for (int i = 0; i < cellsSize; i++) {
      Cell *cell = cells.getDataAt(i);
      if (!cell)
        continue;
      if (cell->isAlive() == false) {
        dead++;
      } else if (cell->isAlive() == true) {
        alive++;
      }
      // if (cell->activate[read_t] == false) f++;
      else if (typeid(*cell) == typeid(Fibroblast)) {
        fibroblastSize++;
      } 
      // else af++;
    }

    int numbiomaterial = 0;
    numbiomaterial = countPatchType(biomaterial);
    cout << " total patches: " << numbiomaterial << endl;
    cout << " alive cells: " << alive << endl;
    cout << " dead cells: " << dead << endl;
    cout << " total cells: " << cells.actualSize() << endl;
    cout << " fibroblast cells: " << fibroblastSize << endl;
  }

  /*
   * Steps:
   *   1. If OMP, add cells from thread-local lists to corresponding global
   * lists
   *   2. Perform updates
   */
  void BMWorld::updateCellsInitial() {
// Cell lists should be empty
// Add new cells
#ifdef _OMP
    cerr << "	updateCell() _OMP" << endl;
    // TODO: parallelize

    int numThreads = omp_get_num_threads();
    for (int tid = 0; tid < numThreads; tid++) {
      /* ------------------------------ Cells ------------------------------ */
      vector<Cell *> *fvec_ptr = localNewCells[tid];
      for (vector<Cell *>::iterator cell_it = fvec_ptr->begin();
           cell_it != fvec_ptr->end(); cell_it++) {
        Cell *newCell = *cell_it;
        if (!cells.addData(newCell, tid)) {
          cerr << "Error: Could not add cell" << endl;
          exit(-1);
        }
      }
      fvec_ptr->clear();
    }
#endif

    /* ------------------------------ Cells ------------------------------ */
    // No need for deletion since these are new cells
    int cellsSize = cells.size();
#pragma omp parallel for
    for (int i = 0; i < cellsSize; i++) {
#ifdef _OMP
      int tid = omp_get_thread_num();
#else
    int tid = DEFAULT_TID;
#endif
      Cell *cell = cells.getDataAt(i);
      if (!cell)
        continue;
      cell->updateAgent();
      /* Added by MM to check types of cell stages and add to respective
       * counters: */
      if (typeid(*cell) == typeid(Fibroblast)) {
        Fibroblast::numOfFibroblast++;
      } 
    }
    Cell::numOfCells = cells.actualSize();
    BMWorld::liveCells = static_cast<float>(cells.actualSize());
    /* Day-0 dead cells implied by the measured viability (LIVE/DEAD, 2.2.7). */
    const double v0 = this->initialViability;
    BMWorld::deadCells = (v0 > 0.0 && v0 < 1.0)
        ? static_cast<float>(BMWorld::liveCells * (1.0 - v0) / v0) : 0.f;
  }

  int BMWorld::userInput() {
    const std::string config_path = util::getSimulationConfigPath();
    const WorldInitParams cfg = load_world_init_config(config_path);

    this->initialViability = cfg.initial_viability;
    this->initialCells.assign(1, cfg.fibroblast_count);
    this->initialCollagenPerPatch = cfg.initial_ecm.collagen_per_patch;
    this->initialElastinPerPatch  = cfg.initial_ecm.elastin_per_patch;
    this->initialHAPerPatch       = cfg.initial_ecm.ha_per_patch;

    /* Section 2.2.3: equal-concentration CMHA-S and Gtn-DTPH stocks are mixed
     * at a volumetric ratio r:1, then PEGDA is added to a final % w/v.
     * The Table 2 world variables follow from that recipe; keeping the
     * conversion here means the config never carries a derived quantity. */
    const double r       = cfg.biomaterial.ha_gtn_ratio;
    const double ha_wv   = cfg.biomaterial.ha_wv_percent  * (r / (r + 1.0));
    const double gtn_wv  = cfg.biomaterial.gtn_wv_percent * (1.0 / (r + 1.0));
    const double tp_wv   = ha_wv + gtn_wv;
    const double pegda   = cfg.biomaterial.pegda_wv_percent;

    BMWorld::HAwv = static_cast<float>(ha_wv);
    BMWorld::TPwv = static_cast<float>(tp_wv);
    /* HA_ww: HA as a mass fraction of the HA-Gtn polymer. */
    BMWorld::HAww = static_cast<float>((tp_wv > 0.0) ? ha_wv / tp_wv : 0.0);
    /* XL_ww: PEGDA as a mass fraction of polymer + crosslinker.
     * NOTE: this denominator is the open question for Vanderhooft 2009 [63];
     * if that fit used PEGDA/polymer instead, drop `+ pegda` below. */
    BMWorld::XLww = static_cast<float>(
        (tp_wv + pegda > 0.0) ? pegda / (tp_wv + pegda) : 0.0);
    BMWorld::TDBMR =
        static_cast<float>(cfg.biomaterial.thiol_double_bond_molar_ratio);

    cout << "-------------------------------------------" << endl;
    cout << "Biomaterial composition (Manuscript Table 2 world variables)" << endl;
    cout << "  Initial fibroblasts             = " << this->initialCells[0] << endl;
    cout << "  HA : Gtn volumetric ratio       = " << r << " : 1" << endl;
    cout << "  HA concentration  HAwv  (% w/v) = " << BMWorld::HAwv << endl;
    cout << "  Gtn concentration       (% w/v) = " << gtn_wv << endl;
    cout << "  Total polymer     TPwv  (% w/v) = " << BMWorld::TPwv << endl;
    cout << "  HA fraction       HAww  (w/w)   = " << BMWorld::HAww << endl;
    cout << "  PEGDA crosslinker XLww  (w/w)   = " << BMWorld::XLww << endl;
    cout << "  Thiol:double bond TDB_MR        = " << BMWorld::TDBMR << endl;
    cout << "-------------------------------------------" << endl;
    return 0;
  }

  char *BMWorld::get_output_filename() {
    return util::outputFileName;
    // return "output/Output_Biomarkers.csv";
  }

  vector<string> BMWorld::get_agent_type_names() {
    return {"Fibroblast", "Activated Fibroblast"};
  }

  void BMWorld::count_agent_types(map<string, int> &agent_counts) {
    int cellsSize = cells.size();
    for (int i = 0; i < cellsSize; i++) {
      Cell *cell = cells.getDataAt(i);
      if (!cell || !cell->isAlive()) continue;
      if (cell->isActivated()) agent_counts["Activated Fibroblast"]++;
      else                     agent_counts["Fibroblast"]++;
    }
  }

  int BMWorld::get_total_agent_count() { return cells.actualSize(); }

  vector<string> BMWorld::get_env_type_names() {
    return {"ncollagen", "nelastin", "HA", "fHA"};
  }

  void BMWorld::count_env(map<string, float> & env_counts) {
    /* Accumulate in double: summing ~8e6 per-patch floats into a float drifts
     * by a few percent once the total reaches ~1e6 (rounding of each add). */
    double col = 0.0, eln = 0.0, ha = 0.0, fha = 0.0;
    for (int in = 0; in < nx * ny * nz; in++) {
      col += this->worldECM[in].ncollagen[read_t];
      eln += this->worldECM[in].nelastin[read_t];
      ha  += this->worldECM[in].HA[read_t];
      fha += this->worldECM[in].fHA[read_t];
    }
    env_counts["ncollagen"] = static_cast<float>(col);
    env_counts["nelastin"]  = static_cast<float>(eln);
    env_counts["HA"]        = static_cast<float>(ha);
    env_counts["fHA"]       = static_cast<float>(fha);
  }


  void BMWorld::write_data_row(std::ofstream & file,
    std::map<std::string, int> & agent_counts,
    std::map<std::string, float> & env_counts) {
      /* Cytokines - Manuscript Table 2 */
      file << this->world_total_tnf()  << "," << this->world_total_tgf()  << ","
      << this->world_total_fgf()  << "," << this->world_total_il6()  << ","
      << this->world_total_il8()  << "," << this->world_total_il10() << ",";

      /* ECM - Manuscript Table 2 */
      file << fixed << setprecision(5)
      << env_counts["ncollagen"] << "," << env_counts["nelastin"] << ","
      << env_counts["HA"]        << "," << env_counts["fHA"]      << ",";

      /* Cells */
      /* Total = live + cumulative dead: every cell the gel contains, as counted
       * by PicoGreen (2.2.6) and as the LIVE/DEAD denominator (2.2.7). Dead
       * agents are deleted from `cells`, so actualSize() alone equals Live. */
      file << (liveCells + deadCells) << "," << liveCells << "," << deadCells << ","

      << agent_counts["Fibroblast"] << ","
      << agent_counts["Activated Fibroblast"] << ",";

      /* Scaffold - Manuscript Table 4 */
      file << BMWorld::E        << "," << BMWorld::Q         << ","
      << BMWorld::massLoss << "," << BMWorld::pXL       << ","
      << BMWorld::poreWidth<< "," << BMWorld::meshSize  << ",";

      /* Scaffold composition inputs */
      file << BMWorld::HAww << "," << BMWorld::HAwv << ","
      << BMWorld::TPwv << "," << BMWorld::XLww << ",";

      /* Cell behaviour */
      file << calculate_viability() << "," << calculate_mean_displacement() << endl;
  }


  // private helpers
  float BMWorld::calculate_viability() {
    const float total = liveCells + deadCells;
    if (total <= 0.f) return 100.0f;      // 播种 tick，还没有统计数据
    return (liveCells / total) * 100.f;
  }


  // float BMWorld::calculate_pct_differentiated(map<string, int> & agent_counts) {
  //   return (static_cast<float>(agent_counts["NP"]) / get_total_agent_count()) *
  //          100;
  // }

  float BMWorld::calculate_mean_displacement() {
    double sum = 0.0;
    int n = 0;
    const int cellsSize = cells.size();
    for (int i = 0; i < cellsSize; i++) {
      Cell *cell = cells.getDataAt(i);
      if (!cell || !cell->isAlive()) continue;
      sum += cell->displacementFromSeed();
      n++;
    }
    return (n > 0) ? static_cast<float>(sum / n) : 0.f;
  }


  void BMWorld::write_csv_header(ofstream & file) {
    file << "clock (30 min)" << ","          // skeleton
         << "Day" << ","                      // skeleton
         /* Cytokines - Manuscript Table 2 */
         << "Total TNF (pg)" << ","
         << "Total TGF (pg)" << ","
         << "Total FGF (pg)" << ","
         << "Total IL6 (pg)" << ","
         << "Total IL8 (pg)" << ","
         << "Total IL10 (pg)" << ","
         /* ECM - Manuscript Table 2 */
         << "Collagen (pg)" << ","
         << "Elastin (pg)" << ","
         << "HA (pg)" << ","
         << "fHA (pg)" << ","
         /* Cells */
         << "Total Cells (live + dead)" << ","
         << "Live Cells" << ","
         << "Dead Cells (cumulative)" << ","
         << "Unactivated Fibroblast" << ","
         << "Activated Fibroblast" << ","
         /* Scaffold - Manuscript Table 4 */
         << "Elastic Modulus (Pa)" << ","
         << "Swelling Ratio (%)" << ","
         << "Mass Loss (%)" << ","
         << "Crosslink Density (mmol/mL)" << ","
         << "Pore Size (um)" << ","
         << "Mesh Size (1/um)" << ","
         /* Scaffold composition inputs */
         << "HA_ww (%)" << ","
         << "HA_wv (%)" << ","
         << "TP_wv (%)" << ","
         << "XL_ww (%)" << ","
         /* Cell behaviour */
         << "Viability Rate (%)" << ","
         << "Mean Displacement (patches)" << endl;
  }
  
  /* Secondary output: TGF concentration along an x-face line through the
   * centre of the grid, for diffusion diagnostics. */
   void BMWorld::write_auxiliary_header() {
    char tgf_path[512];
    util::makeOutputPath(tgf_path, sizeof(tgf_path), "tgf_line.csv");
    remove(tgf_path);
    ofstream tgf_file(tgf_path, ios::app);
    for (int xi = 0; xi <= nx / 2; xi++)
      tgf_file << "x=" << xi << (xi < nx / 2 ? "," : "\n");
    tgf_file.close();
  }

  void BMWorld::write_auxiliary_outputs() {
    char tgf_path[512];
    util::makeOutputPath(tgf_path, sizeof(tgf_path), "tgf_line.csv");
    ofstream tgf_file(tgf_path, ios::app);
    tgf_file << fixed << setprecision(10);
    for (int xi = 0; xi <= nx / 2; xi++)
      tgf_file << tgfLine[xi] << (xi < nx / 2 ? "," : "\n");
    tgf_file.close();

    // char o2_path[512];
    // util::makeOutputPath(o2_path, sizeof(o2_path), "o2_line.csv");
    // ofstream o2_file(o2_path, ios::app);
    // o2_file << fixed << setprecision(10);
    // for (int xi = 0; xi <= nx / 2; xi++)
    //   o2_file << o2Line[xi] << (xi < nx / 2 ? "," : "\n");
    // o2_file.close();
  }

  void BMWorld::patchassign_csv() {
    int in = 0;
    int Number = 0;
    for (int iz = 0; iz < nz; iz++) {
      char patchassign[512];
      char cells[512];
      char cells_w[512];
      char initcollagen[512];
      char initHA[512];
      char damagezone[512];
      char initialdamage[512];
      char extension[10] = ".csv";
      char tempNumber[32] = "";
      sprintf(tempNumber, "_t%3.0f_z%d", this->clock, Number);
      snprintf(patchassign, sizeof(patchassign), "%s/patchassign%s%s",
               util::getOutputDir(), tempNumber, extension);
      snprintf(cells, sizeof(cells), "%s/cells_read%s%s", util::getOutputDir(),
               tempNumber, extension);
      snprintf(cells_w, sizeof(cells_w), "%s/cells_write%s%s",
               util::getOutputDir(), tempNumber, extension);
      snprintf(initHA, sizeof(initHA), "%s/initHA%s%s", util::getOutputDir(),
               tempNumber, extension);
      snprintf(initcollagen, sizeof(initcollagen), "%s/initcollagen%s%s",
               util::getOutputDir(), tempNumber, extension);
      snprintf(damagezone, sizeof(damagezone), "%s/damagezone%s%s",
               util::getOutputDir(), tempNumber, extension);
      snprintf(initialdamage, sizeof(initialdamage), "%s/initialdamage%s%s",
               util::getOutputDir(), tempNumber, extension);

      // Patch Assign
      ofstream output_file(patchassign, ios::app);
      for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
          in = ix + iy * nx + iz * nx * ny;
          if (worldPatch[in].type[read_t] == damage ||
              worldPatch[in].damage[read_t] != 0) {
            output_file << "x";
            continue;
          }
          if (worldPatch[in].type[read_t] == nothing) {
            output_file << "-";
          }
          if (worldPatch[in].type[read_t] == unidentifiable) {
            output_file << "?";
          }
        }
        output_file << endl;
      }
      output_file.close();

      // initHA
      ofstream output_file1(initHA, ios::app);
      for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
          in = ix + iy * nx + iz * nx * ny;
          if (worldPatch[in].initHA == true) {
            output_file1 << "u";
            continue;
          }
          if (worldPatch[in].type[read_t] == damage ||
              worldPatch[in].damage[read_t] != 0) {
            output_file1 << "x";
            continue;
          }
          if (worldPatch[in].type[read_t] == nothing) {
            output_file1 << "-";
          }
          if (worldPatch[in].type[read_t] == unidentifiable) {
            output_file1 << "?";
          }
        }
        output_file1 << endl;
      }
      output_file1.close();

      // Damage Zone
      ofstream output_file2(damagezone, ios::app);
      for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
          in = ix + iy * nx + iz * nx * ny;

          if (worldPatch[in].inDamzone == true) {
            output_file2 << "z";
            continue;
          }
          if (worldPatch[in].type[read_t] == damage ||
              worldPatch[in].damage[read_t] != 0) {
            output_file2 << "x";
            continue;
          }
          if (worldPatch[in].type[read_t] == nothing) {
            output_file2 << "-";
          }
          if (worldPatch[in].type[read_t] == unidentifiable) {
            output_file2 << "?";
          }
        }
        output_file2 << endl;
      }
      output_file2.close();

      // Initial Damage
      ofstream output_file3(initialdamage, ios::app);
      for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
          in = ix + iy * nx + iz * nx * ny;
          if (worldPatch[in].type[read_t] == damage ||
              worldPatch[in].damage[read_t] != 0) {
            output_file3 << "x";
            continue;
          }
          if (worldPatch[in].type[read_t] == nothing) {
            output_file3 << "-";
          }
          if (worldPatch[in].type[read_t] == unidentifiable) {
            output_file3 << "?";
          }
        }
        output_file3 << endl;
      }
      output_file3.close();

      // Cells
      ofstream output_file6(cells, ios::app);
      for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
          in = ix + iy * nx + iz * nx * ny;
          if (worldPatch[in].isOccupied()) {
            if (worldPatch[in].occupiedby[read_t] == fibroblast) {
              output_file6 << "f";
              continue;
            }
          }
          if (worldPatch[in].type[read_t] == nothing) {
            output_file6 << "-";
          }
          if (worldPatch[in].type[read_t] == unidentifiable) {
            output_file6 << "?";
          }
        }
        output_file6 << endl;
      }
      output_file6.close();

      ofstream output_file7(cells_w, ios::app);
      for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
          in = ix + iy * nx + iz * nx * ny;
          if (worldPatch[in].isOccupiedWrite()) {
            if (worldPatch[in].occupiedby[read_t] == fibroblast) {
              output_file6 << "f";
              continue;
            } 
          }
          if (worldPatch[in].type[write_t] == nothing) {
            output_file7 << "-";
          }
          if (worldPatch[in].type[write_t] == unidentifiable) {
            output_file7 << "?";
          }
        }
        output_file7 << endl;
      }
      output_file7.close();

      Number++;
    }
  }