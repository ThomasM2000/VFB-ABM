/* 
 * File: Cell.cpp
 *
 * File Contents: Contains the Cell class.
 *
 * Author: Yvonna
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 *               Meghana Munipalle
 */

#include "Cell.h"
#include "../../enums.h"
#include <iostream>
#include <algorithm>      
#include <cmath>  

using namespace std;
int Cell::numOfCells = 0;

int Fibroblast::numOfFibroblast  = 0;
int Fibroblast::numOfAFibroblast = 0;

float Fibroblast::migration      [Fibroblast::MIGRATION_COUNT]       = {};
float Fibroblast::viability      [Fibroblast::VIABILITY_COUNT]       = {};
float Fibroblast::proliferation  [Fibroblast::PROLIFERATION_COUNT]   = {};
float Fibroblast::tgfSynthesis   [Fibroblast::TGF_SYNTHESIS_COUNT]   = {};
float Fibroblast::fgfSynthesis   [Fibroblast::FGF_SYNTHESIS_COUNT]   = {};
float Fibroblast::tnfSynthesis   [Fibroblast::TNF_SYNTHESIS_COUNT]   = {};
float Fibroblast::il6Synthesis   [Fibroblast::IL6_SYNTHESIS_COUNT]   = {};
float Fibroblast::il8Synthesis   [Fibroblast::IL8_SYNTHESIS_COUNT]   = {};
float Fibroblast::activation_params[Fibroblast::ACTIVATION_COUNT]    = {};
float Fibroblast::collagenSynth  [Fibroblast::COLLAGEN_SYNTH_COUNT]  = {};
float Fibroblast::elastinSynth   [Fibroblast::ELASTIN_SYNTH_COUNT]   = {};
float Fibroblast::haSynth        [Fibroblast::HA_SYNTH_COUNT]        = {};

static inline float safe_denominator(float d) { return (d > 1e-6f) ? d : 1e-6f; }
static inline float safe_log10_arg(float a)   { return (a > 1.f) ? a : 1.f; }
static inline double sim_days_at_least_one() {
  double t = BMWorld::reportDay();
  return (t >= 1.0) ? t : 1.0;  // ln(t) is undefined at t = 0
}
/* Number of `interval_h`-hour boundaries crossed during the current tick.
 * Replaces fmod(hour, k) == 0, which only works when k is a multiple of the
 * 0.5 h tick and never fires for values such as k6 = 23.7 during SA/NM. */
static inline int intervals_elapsed(double interval_h) {
	if (interval_h <= 0.0) return 0;
	const double tick_h = Agent::agentWorldPtr->tick_interval_minutes() / 60.0;
	const double now  = BMWorld::reportHour();
	const double prev = now - tick_h;
	if (prev < 0.0) return 0;
	const double eps = 1e-9;
	return static_cast<int>(std::floor(now  / interval_h + eps)
						  - std::floor(prev / interval_h + eps));
}
  

//DEFAULT CONSTRUCTORS
Cell::Cell() {
	cout << "default cell alloc" << endl;
}

Fibroblast::Fibroblast() {
	cout << "default fibroblast alloc" << endl;
}



//CONSTRUCTORS FOR PATCH POINTERS
Cell::Cell(Patch* patchPtr) {
	this->ix[write_t] = patchPtr->indice[0];
	this->iy[write_t] = patchPtr->indice[1];
	this->iz[write_t] = patchPtr->indice[2];
	this->index[write_t] = patchPtr->index;
	this->ix0 = patchPtr->indice[0];
	this->iy0 = patchPtr->indice[1];
	this->iz0 = patchPtr->indice[2];
	this->alive[write_t] = true;
	this->realDeath[write_t] = false;
	this->doublings[write_t] = 0;


	this->activate[write_t] = false;
	//this->color[write_t] = ccell;
	this->size[write_t] = 2;
	//this->type[write_t] = cell;
	this->ix[read_t] = patchPtr->indice[0];
	this->iy[read_t] = patchPtr->indice[1];
	this->iz[read_t] = patchPtr->indice[2];
	this->index[read_t] = patchPtr->index;

	/* In OMP version, we wait to add cells at the end of the tick, whereas in serial version, cells are always added right away.
     * Thus, cells, when added in OMP, should be alive right away. */
	#ifdef _OMP
		this->alive[read_t] = true;
		this->life[read_t] = this->life[write_t];
	#else
		this->alive[read_t] = false;
		this->life[read_t] = 0;
	#endif
	this->life[write_t] = 0;
	
	this->activate[read_t] = false;
	//this->color[read_t] = ccell;
	this->size[read_t] = 2;
	//this->type[read_t] = cell;
}

Fibroblast::Fibroblast(Patch* patchPtr) : Cell(patchPtr) {
	this->alive[write_t] = true;
	this->realDeath[write_t] = false;
	this->doublings[write_t] = 0;
	this->color[write_t] = cfibroblast;
	this->type[write_t] = fibroblast;

	this->doublings[read_t] = 0;
	this->color[read_t] = cfibroblast;
	this->type[read_t] = fibroblast;
}



//CONSTRUCTORS WITH COORDINATES
Cell::Cell(int x, int y, int z) {
	this->ix[write_t] = x;
	this->iy[write_t] = y;
	this->iz[write_t] = z;
	this->index[write_t] = x + y*nx + z*nx*ny;
	this->ix0 = x;
	this->iy0 = y;
	this->iz0 = z;
	this->alive[write_t] = true;
	this->realDeath[write_t] = false;
	this->doublings[write_t] = 0;


	this->life[write_t] = 0;

	this->activate[write_t] = false;
	this->color[write_t] = ccell;
	this->size[write_t] = 2;
	this->type[write_t] = cell;

  /* In OMP version, we wait to add cells at the end of the tick, whereas in serial version, cells are always added right away.
   * Thus, cells, when added in OMP, should be alive right away. */
	#ifdef _OMP
		this->ix[read_t] = x;
		this->iy[read_t] = y;
		this->iz[read_t] = z;
		this->index[read_t] = x + y*nx + z*nx*ny;
		this->alive[read_t] = true;
		this->life[read_t] = this->life[write_t];
		this->activate[read_t] = false;
		this->color[read_t] = ccell;
		this->size[read_t] = 2;
		this->type[read_t] = cell;
		this->doublings[write_t] = 0;
		this->realDeath[read_t] = false;
	#else
		this->ix[read_t] = x;
		this->iy[read_t] = y;
		this->iz[read_t] = z;
		this->index[read_t] = x + y*nx + z*nx*ny;
		this->alive[read_t] = false;
		this->life[read_t] = 0;
		this->activate[read_t] = false;
		this->color[read_t] = ccell;
		this->size[read_t] = 2;
		this->type[read_t] = cell;
		this->doublings[write_t] = 0;
		this->realDeath[read_t] = false;
	#endif
}

Fibroblast::Fibroblast(int x, int y, int z) : Cell(x, y, z) {
	this->type[read_t]       = fibroblast;
	this->type[write_t]      = fibroblast;
	this->color[read_t]      = cfibroblast;
	this->color[write_t]     = cfibroblast;
	this->doublings[read_t]  = 0;
	this->doublings[write_t] = 0;
	this->realDeath[read_t]  = false;
	this->realDeath[write_t] = false;
}


//DESTRUCTORS
Cell::~Cell() {}

Fibroblast::~Fibroblast() {}

/* -------------------------------------------------------------------------- */
/*                            ECM DEPOSITION HELPERS                          */
/* -------------------------------------------------------------------------- */
namespace {
	/** Pick a random in-bounds Moore neighbour of the agent, or -1. */
	int random_neighbour_index(int x, int y, int z, unsigned *seed) {
	  vector<int> neighbours;
	  for (int i = 0; i < 27; i++) {
		int dx = Agent::dX[i], dy = Agent::dY[i], dz = Agent::dZ[i];
		if (x + dx < 0 || x + dx >= Agent::nx) continue;
		if (y + dy < 0 || y + dy >= Agent::ny) continue;
		if (z + dz < 0 || z + dz >= Agent::nz) continue;
		neighbours.push_back((x + dx) + (y + dy) * Agent::nx +
							 (z + dz) * Agent::nx * Agent::ny);
	  }
	  if (neighbours.empty()) return -1;
	  return neighbours[rand_r(seed) % neighbours.size()];
	}
}  // namespace

void Cell::depositCollagen(float amount) {
	if (amount <= 0) return;
	int tid = 0;
#ifdef _OMP
	tid = omp_get_thread_num();
#endif
	int in = random_neighbour_index(this->ix[read_t], this->iy[read_t],
									this->iz[read_t],
									&(Agent::agentWorldPtr->seeds[tid]));
	if (in < 0) return;
	Agent::agentECMPtr[in].ncollagen[write_t] += amount;
#ifdef OPT_ECM
	Agent::agentECMPtr[in].set_dirty();
#endif
}

void Cell::depositElastin(float amount) {
	if (amount <= 0) return;
	int tid = 0;
#ifdef _OMP
	tid = omp_get_thread_num();
#endif
	int in = random_neighbour_index(this->ix[read_t], this->iy[read_t],
									this->iz[read_t],
									&(Agent::agentWorldPtr->seeds[tid]));
	if (in < 0) return;
	Agent::agentECMPtr[in].nelastin[write_t]  += amount; 
#ifdef OPT_ECM
	Agent::agentECMPtr[in].set_dirty();
#endif
}

void Cell::depositHA(float amount, int here) {
	if (amount <= 0) return;
	int in;
	if (here) {
	in = isModified(this->index) ? this->index[write_t] : this->index[read_t];
	} else {
	int tid = 0;
#ifdef _OMP
	tid = omp_get_thread_num();
#endif
	in = random_neighbour_index(this->ix[read_t], this->iy[read_t],
								this->iz[read_t],
								&(Agent::agentWorldPtr->seeds[tid]));
	}
	if (in < 0) return;
	Agent::agentECMPtr[in].HA[write_t]        += amount;
#ifdef OPT_ECM
	Agent::agentECMPtr[in].set_dirty();
#endif
}

//CELL FUNCTIONS
/* Table 3, rules 9 and 10. */
void Cell::activation() {
	int in = this->index[read_t];
	float patchTGF = this->patchChemConcentration(TGF, in);

	if (this->activate[read_t] == false) {
		if (this->should_activate(patchTGF)) {
			this->activate[write_t] = true;
			this->type[write_t] = afibroblast;
			this->color[write_t] = cafibroblast;
			Agent::agentPatchPtr[in].occupiedby[write_t] = afibroblast;
			Agent::agentPatchPtr[in].dirty = true;
			Fibroblast::numOfFibroblast--;
			Fibroblast::numOfAFibroblast++;
		}
	} else {
		if (this->should_deactivate()) {
			this->activate[write_t] = false;
			this->type[write_t] = fibroblast;
			this->color[write_t] = cfibroblast;
			Agent::agentPatchPtr[in].occupiedby[write_t] = fibroblast;
			Agent::agentPatchPtr[in].dirty = true;
			Fibroblast::numOfAFibroblast--;
			Fibroblast::numOfFibroblast++;
		}
	}
}

void Cell::cytokine_synthesis() {
	this->create_cytokines();
}

/* Table 3, rules 11-13. Each hook owns its own synthesis interval. */
void Cell::ecm_synthesis() {
	int in = this->index[read_t];
	if (Agent::agentPatchPtr[in].type[read_t] != biomaterial) return;
  
	this->create_collagen();
	this->create_elastin();
	this->create_ha();
}


void Cell::cellFunction() {
	if (this->alive[read_t] == false) return;

	/* 1. Activation / deactivation  (Table 3, rules 9-10) */
	this->activation();

	/* 2. Migration  (Table 3, rule 1) */
	this->cellSniff();

	/* 3. Death  (Table 3, rule 2).
	 * Runs before proliferation: apoptose() compares against the
	 * liveCells/deadCells aggregate computed at the end of the previous tick,
	 * so it must act on that same population. With proliferation first the
	 * aggregate is stale by one day's births and the viability ratio drifts
	 * above vr(t) by a growing margin. A cell that dies this tick also no
	 * longer divides this tick. */
	this->apoptose();
	if (this->alive[write_t] == false) return;

	/* 4. Proliferation  (Table 3, rule 3) */
	this->proliferate();

	/* 5-6. Only activated fibroblasts secrete (Figure 1C) */
	if (this->activate[write_t]) {
		this->cytokine_synthesis();
		this->ecm_synthesis();
	}

	if (this->life[read_t] >= 0)
		this->life[write_t] = this->life[read_t] + 1;
}


void Cell::cellSniff() {
	int speed = static_cast<int>(this->get_migration_speed());
	if (speed < 1) return;                 // slower than one patch this tick
  
	if (this->moveTowardChemotaxis(speed)) return;
  
	for (int step = 0; step < speed; step++)
	  	this->wiggle();
}

void Cell::die() {
	int in = this->index[read_t];
	Agent::agentPatchPtr[in].clearOccupied();
	Agent::agentPatchPtr[in].occupiedby[write_t] = nothing;
	this->alive[write_t] = false;
	this->life[write_t] = -1;
}

/*
 * Table 3, rule 2. Chen & Thibeault 2010 [54] Section 3.4 / Fig. 5 measure
 * vr(t) as "the percentage of live cells in the total cell population"
 * (63.7 % at day 3, 67 % at day 7) -- a snapshot ratio, not a daily hazard.
 * Dead cells stay in the gel and stay in the denominator, which is what
 * BMWorld::liveCells / (liveCells + deadCells) tracks. So treat vr(t) as a
 * target live fraction and kill only the excess:
 *     p_die = 1 - target / current_fraction
 */
void Cell::apoptose() {
	/* No daily gate: vr(t) is defined continuously and the targeting formula
	 * drives the live fraction to it at whatever interval it is evaluated.
	 * Correcting once a day lets a full day of births accumulate first, which
	 * biases the sampled viability 2-4 points above vr(t). */
	const double t = BMWorld::reportDay();
	// if (t < 1.0) return;                  // ln(t) undefined below one day

	const double total = static_cast<double>(BMWorld::liveCells)
					   + static_cast<double>(BMWorld::deadCells);
	if (total <= 0.0) return;
	const double current = static_cast<double>(BMWorld::liveCells) / total;
	if (current <= 0.0) return;

	const double target = this->get_viability_rate(t) / 100.0;
	if (current <= target) return;

	const double p_die = 1.0 - target / current;
	if (Agent::rollDice(static_cast<float>(p_die * 100.0))) {
		this->realDeath[write_t] = true;
		this->die();
	}
}



void Cell::copyAndInitialize(Agent* original, int dx, int dy, int dz) {
	int in = this->index[read_t];

	// Initializes location of new Cell relative to original agent:
	this->ix[write_t] = original->getX() + dx;
	this->iy[write_t] = original->getY() + dy;
	this->iz[write_t] = original->getZ() + dz;
	this->index[write_t] = this->ix[write_t] + this->iy[write_t]*Agent::nx + this->iz[write_t]*Agent::nx*Agent::ny;
  	
	// Initializes new Cell:
	this->alive[read_t] = true;
	this->life[read_t] = 0; // 0 corresponds to ticks
	this->activate[read_t] = false;
	this->color[read_t]= ccell;
	this->size[read_t] = 2;
	this->type[read_t] = cell;
	this->alive[write_t] = true;
	this->life[write_t] = this->life[read_t];
	this->activate[write_t] = false;
	this->color[write_t]= ccell;
	this->size[write_t] = 2;
	this->type[write_t] = cell;

	Cell::numOfCells++;

  	// Assigns new Cell to this patch if it is unoccupied:
	if (Agent::agentPatchPtr[in].isOccupied() == false) {
		Agent::agentPatchPtr[in].setOccupied();
		Agent::agentPatchPtr[in].occupiedby[write_t] = this->type[read_t];
	} else {
		cout << "error in hatching and initialization!!!" << dx << " " << dy << endl;
	}
}

int Cell::get_max_doublings() { return 100; }   // ihVFF, passage 6-10

void Cell::proliferate() {
	int in = this->index[read_t];
	if (Agent::agentPatchPtr[in].type[read_t] != biomaterial) return;
  
	// const float hours_between = Fibroblast::proliferation[Fibroblast::PROLIFERATION_HOURS_BETWEEN]; // k6
	// if (hours_between <= 0) return;
	// if (fmod(BMWorld::reportHour(), hours_between) != 0) return;
	// if (BMWorld::reportHour() == 0) return;  // no division on the seeding tick
	if (intervals_elapsed(Fibroblast::proliferation[Fibroblast::PROLIFERATION_HOURS_BETWEEN]) == 0) return; // k6

	if (this->doublings[read_t] >= this->get_max_doublings()) return;

	float prob = this->get_prolif_prob();
	if (prob <= 0) return;
  
	if (Agent::rollDice(prob) && this->hatchnewcell(1, this->type[read_t]) > 0)
	  this->doublings[write_t] = this->doublings[read_t] + 1;

}

float Cell::get_migration_speed()          { return 0; }
float Cell::get_viability_rate(double)     { return 100; }
float Cell::get_prolif_prob()              { return 0; }
bool  Cell::should_activate(float)         { return false; }
bool  Cell::should_deactivate()            { return false; }
void  Cell::create_cytokines()             {}
void  Cell::create_collagen()              {}
void  Cell::create_elastin()               {}
void  Cell::create_ha()                    {}


int Cell::hatchnewcell(int number, int agentType, int here) {
	int newcells = 0;
	const int x = this->ix[read_t];
	const int y = this->iy[read_t];
	const int z = this->iz[read_t];
	const int nx = Agent::nx, ny = Agent::ny, nz = Agent::nz;

#ifdef MODEL_3D
	random_shuffle(&Agent::neighbor[0], &Agent::neighbor[27]);
	const int nNeighbors = 27;
#else
	random_shuffle(&Agent::neighbor[0], &Agent::neighbor[8]);
	const int nNeighbors = 8;
#endif

	for (int i = 0; i < nNeighbors && newcells < number; i++) {
		int lx, ly, lz, in;

		if (here == 0) {
			const int dx = Agent::dX[Agent::neighbor[i]];
			const int dy = Agent::dY[Agent::neighbor[i]];
			const int dz = Agent::dZ[Agent::neighbor[i]];
			lx = x + dx; ly = y + dy; lz = z + dz;
			if (lx < 0 || lx >= nx || ly < 0 || ly >= ny || lz < 0 || lz >= nz) continue;
			in = lx + ly * nx + lz * nx * ny;
			/* Only divide into an unoccupied biomaterial patch. */
			if (Agent::agentPatchPtr[in].type[read_t] != biomaterial) continue;
			if (Agent::agentPatchPtr[in].isOccupied()) continue;
		} else {
			lx = x; ly = y; lz = z;
			in = this->getIndex();
		}

		Cell *newcell = nullptr;
		switch (agentType) {
		  case fibroblast:
		  case afibroblast:
			newcell = new Fibroblast(lx, ly, lz);
			break;
		  default:
			continue;
		}
		if (!newcell) continue;

		Agent::agentPatchPtr[in].setOccupied();
		Agent::agentPatchPtr[in].occupiedby[write_t] = fibroblast;
		Agent::agentPatchPtr[in].dirty = true;

#ifdef _OMP
		Agent::agentWorldPtr->localNewCells[omp_get_thread_num()]->push_back(newcell);
#else
		Agent::agentWorldPtr->cells.addData(newcell, DEFAULT_TID);
#endif
		newcells++;
	}
	return newcells;
}

/* -------------------------------------------------------------------------- */
/*                                    FIBROBLAST                                    */
/* -------------------------------------------------------------------------- */


// Table 3, rule 1:  v = -k1 ln(E) + k2   [um/min] -> patches/tick.
float Fibroblast::get_migration_speed() {
	float E = BMWorld::E;                       // elastic modulus in Pa
	if (E <= 1.f) E = 1.f;                      // ln(E) domain guard

	float v_um_per_min = -Fibroblast::migration[MIGRATION_ELASTICITY_EFFECT] * log(E)
                       + Fibroblast::migration[MIGRATION_BASELINE_SPEED];
	if (v_um_per_min < 0.f) v_um_per_min = 0.f;

	const double tick_min = Agent::agentWorldPtr->tick_interval_minutes();
	const double patch_um = Agent::agentWorldPtr->patchlength * 1000.0; // mm -> um
	const float patches = static_cast<float>(v_um_per_min * tick_min / patch_um);

  /* Sub-patch speeds are realised stochastically so that, on average, the cell
   * covers `patches` patches per tick. The value is per-agent, so it is not
   * cached in a shared static (which would race under OpenMP). */
	return Agent::rollDice((patches - floor(patches)) * 100.f) ? ceil(patches)
                                                             : floor(patches);
}

/* Table 3, rule 2:  vr = k3 ln(t) + k4  [% live cells in the total population]. */
float Fibroblast::get_viability_rate(double t_days) {
	if (t_days < 1.0) t_days = 1.0;          // ln(t) undefined below one day
	return Fibroblast::viability[VIABILITY_TIME_EFFECT] * static_cast<float>(log(t_days))
		 + Fibroblast::viability[VIABILITY_BIOMATERIAL_EFFECT];
}


/*
 * Table 3, rule 3:
 *   x        = -1 if TGF > k5 else +1
 *   Every k6 hours,
 *   Prolif_b = (k7 - k8 HAww) log(t) + k9 HAww - k10
 *   Prolif   = Prolif_b ( k11 + (HAf + k12 log10(1+TNF+xTGF+FGF+HAf) + k13)/k14 )
 */
float Fibroblast::get_prolif_prob() {
	const int in = this->index[read_t];

	const float patchTGF = this->patchChemConcentration(TGF, in);
	const float patchTNF = this->patchChemConcentration(TNF, in);
	const float patchFGF = this->patchChemConcentration(FGF, in);
	const float HAf = Agent::agentECMPtr[in].fHA[read_t];

	const float x = (patchTGF > Fibroblast::proliferation[PROLIFERATION_TGF_THRESHOLD])
						? -1.f : 1.f;                                   // k5

	const double t = sim_days_at_least_one();
	const float HAww = BMWorld::HAww;

	const float prolif_b =
		(Fibroblast::proliferation[PROLIFERATION_TIME_EFFECT]                       // k7
			- Fibroblast::proliferation[PROLIFERATION_HA_TIME_EFFECT] * HAww)          // k8
			* static_cast<float>(log(t))
		+ Fibroblast::proliferation[PROLIFERATION_HA_EFFECT] * HAww                 // k9
		- Fibroblast::proliferation[PROLIFERATION_BIOMATERIAL_BASELINE];            // k10

	const float chem_arg =
		safe_log10_arg(1.f + patchTNF + x * patchTGF + patchFGF + HAf);

	const float chem_term =
		(HAf
			+ Fibroblast::proliferation[PROLIFERATION_ACTIVATING_FACTOR_EFFECT]        // k12
				* log10(chem_arg)
			+ Fibroblast::proliferation[PROLIFERATION_BASELINE_OFFSET])                // k13
		/ safe_denominator(Fibroblast::proliferation[PROLIFERATION_CHEMICAL_EFFECT]);// k14

	return prolif_b * (Fibroblast::proliferation[PROLIFERATION_BASELINE] + chem_term); // k11
}

/*
 * Table 3, rule 9: activation.
 *   if TGF > k27  -> activate with probability k28
 *   if TGF > k29  -> activate
 *   otherwise     -> activate with probability k30
 * Probabilities are expressed in percent, as rollDice() expects.
 */
bool Fibroblast::should_activate(float patchTGF) {
	if (patchTGF > Fibroblast::activation_params[ACTIVATION_TGF_INTERMEDIATE_THRESHOLD]) // k27
	return Agent::rollDice(
		Fibroblast::activation_params[ACTIVATION_INTERMEDIATE_PROBABILITY]);           // k28

	if (patchTGF > Fibroblast::activation_params[ACTIVATION_TGF_HIGH_THRESHOLD])         // k29
	return true;

	return Agent::rollDice(
		Fibroblast::activation_params[ACTIVATION_INDEPENDENT_PROBABILITY]);              // k30
}

/* Table 3, rule 10: deactivate with probability k31. */
bool Fibroblast::should_deactivate() {
	return Agent::rollDice(
		Fibroblast::activation_params[ACTIVATION_DEACTIVATION_PROBABILITY]);             // k31
}

/*
 * Table 3, rules 4-8: cytokine synthesis, in pg/tick, added to the per-tick
 * secretion channel of the agent's own patch.
 */
void Fibroblast::create_cytokines() {
	const int in = this->index[read_t];
  
	const float lTNF  = this->patchChemConcentration(TNF, in);
	const float lTGF  = this->patchChemConcentration(TGF, in);
	const float lFGF  = this->patchChemConcentration(FGF, in);
	const float lIL10 = this->patchChemConcentration(IL10, in);
	const float lHA   = Agent::agentECMPtr[in].HA[read_t];
	(void)lFGF;
  
	/* TGFinc = k15 + k16 (k17 + TNF + IL10) */
	const float tgfinc =
		Fibroblast::tgfSynthesis[TGF_BASELINE]
		+ Fibroblast::tgfSynthesis[TGF_FEEDBACK_EFFECT]
			  * (Fibroblast::tgfSynthesis[TGF_FEEDBACK_BASELINE] + lTNF + lIL10);
  
	/* FGFinc = k18 */
	const float fgfinc = Fibroblast::fgfSynthesis[FGF_RATE];
  
	/* TNFinc = k19 / (k20 + TGF + k21 IL10 + HA) */
	const float tnfinc =
		Fibroblast::tnfSynthesis[TNF_BASELINE]
		/ safe_denominator(Fibroblast::tnfSynthesis[TNF_INHIBITORY_EFFECT] + lTGF
						   + Fibroblast::tnfSynthesis[TNF_IL10_EFFECT] * lIL10 + lHA);
  
	/* IL6inc = k22 (k23 + TNF) / (k24 + IL10) */
	const float il6inc =
		Fibroblast::il6Synthesis[IL6_BASELINE]
		* (Fibroblast::il6Synthesis[IL6_ACTIVATING_EFFECT] + lTNF)
		/ safe_denominator(Fibroblast::il6Synthesis[IL6_INHIBITORY_EFFECT] + lIL10);
  
	/* IL8inc = k25 / (k26 + IL10 + HA) */
	const float il8inc =
		Fibroblast::il8Synthesis[IL8_BASELINE]
		/ safe_denominator(Fibroblast::il8Synthesis[IL8_INHIBITORY_EFFECT] + lIL10 + lHA);
  
	this->addPatchChemSecretion(TGF, in, tgfinc);
	this->addPatchChemSecretion(FGF, in, fgfinc);
	this->addPatchChemSecretion(TNF, in, tnfinc);
	this->addPatchChemSecretion(IL6, in, il6inc);
	this->addPatchChemSecretion(IL8, in, il8inc);
	/* IL-10 has no synthesis rule in Table 3; it enters only as a baseline. */
}

/*
 * Table 3, rule 11: collagen synthesis, every k32 hours.
 *   Col_b = k33 - k34 mesh - k35 E
 *   Col   = Col_b k36 (log10(1+TGF+IL6) + k37) / (1 + FGF + k38 IL8)
 *   HAf > HAn  -> Col + k39
 *   HAf == HAn -> Col/k40 + k41
 *   else       -> k42 + Col/k43
 */
void Fibroblast::create_collagen() {
	const int n_events = intervals_elapsed(Fibroblast::collagenSynth[COLLAGEN_HOURS_BETWEEN]); // k32
	if (n_events == 0) return;
  
	const int in = this->index[read_t];
	const float lTGF = this->patchChemConcentration(TGF, in);
	const float lIL6 = this->patchChemConcentration(IL6, in);
	const float lFGF = this->patchChemConcentration(FGF, in);
	const float lIL8 = this->patchChemConcentration(IL8, in);
  
	const float col_b = Fibroblast::collagenSynth[COLLAGEN_BASELINE]                 // k33
						- Fibroblast::collagenSynth[COLLAGEN_MESH_EFFECT] * BMWorld::meshSize  // k34
						- Fibroblast::collagenSynth[COLLAGEN_ELASTICITY_EFFECT] * BMWorld::E;  // k35
  
	const float col =
		col_b * Fibroblast::collagenSynth[COLLAGEN_ACTIVATING_STRENGTH]              // k36
		* (log10(safe_log10_arg(1.f + lTGF + lIL6))
		   + Fibroblast::collagenSynth[COLLAGEN_ACTIVATING_OFFSET])                  // k37
		/ safe_denominator(1.f + lFGF
						   + Fibroblast::collagenSynth[COLLAGEN_IL8_EFFECT] * lIL8); // k38
  
	const float HAf = Agent::agentECMPtr[in].fHA[read_t];
	const float HAn = Agent::agentECMPtr[in].HA[read_t];
  
	float rate;
	if (HAf > HAn) {
	  rate = col + Fibroblast::collagenSynth[COLLAGEN_HIGH_FRAG_RATIO_EFFECT];       // k39
	} else if (HAf == HAn) {
	  rate = col / safe_denominator(Fibroblast::collagenSynth[COLLAGEN_EQUAL_FRAG_RATIO_EFFECT]) // k40
			 + Fibroblast::collagenSynth[COLLAGEN_EQUAL_FRAG_OFFSET];                // k41
	} else {
	  rate = Fibroblast::collagenSynth[COLLAGEN_BASELINE_RATE]                       // k42
			 + col / safe_denominator(Fibroblast::collagenSynth[COLLAGEN_CHEMICAL_EFFECT]); // k43
	}
  
	this->depositCollagen(rate * n_events);
}

/*
 * Table 3, rule 12: elastin synthesis, on the same k32-hour clock.
 *   Eln_b = k44 - k45 mesh - k46 E
 *   Eln   = Eln_b ( (k47 log10(1+TGF) + k48) / (k49 (1+FGF+TNF)) + k50 )
 */
void Fibroblast::create_elastin() {
	const int n_events = intervals_elapsed(Fibroblast::collagenSynth[COLLAGEN_HOURS_BETWEEN]); // k32
	if (n_events == 0) return;
  
	const int in = this->index[read_t];
	const float lTGF = this->patchChemConcentration(TGF, in);
	const float lFGF = this->patchChemConcentration(FGF, in);
	const float lTNF = this->patchChemConcentration(TNF, in);
  
	const float eln_b = Fibroblast::elastinSynth[ELASTIN_BASELINE]                    // k44
						- Fibroblast::elastinSynth[ELASTIN_MESH_EFFECT] * BMWorld::meshSize // k45
						- Fibroblast::elastinSynth[ELASTIN_ELASTICITY_EFFECT] * BMWorld::E; // k46
  
	const float eln =
		eln_b * ((Fibroblast::elastinSynth[ELASTIN_ACTIVATING_STRENGTH]               // k47
					  * log10(safe_log10_arg(1.f + lTGF))
				  + Fibroblast::elastinSynth[ELASTIN_ACTIVATING_OFFSET])              // k48
					 / safe_denominator(Fibroblast::elastinSynth[ELASTIN_INHIBITORY_STRENGTH] // k49
										* (1.f + lFGF + lTNF))
				 + Fibroblast::elastinSynth[ELASTIN_BASELINE_RATE]);                  // k50
  
	this->depositElastin(eln * n_events);
}

/*
 * Table 3, rule 13: hyaluronic acid synthesis, every k51 hours.
 *   HA_b = k52 - k53 mesh + k54 E
 *   HA   = HA_b (k55 log10(1+TGF+TNF+FGF) + k56)
 *   with probability HA + k57      -> move, then deposit HA
 *   with probability k58 + HA/k59  -> deposit HA on the current patch
 */
void Fibroblast::create_ha() {
	const int n_events = intervals_elapsed(Fibroblast::haSynth[HA_HOURS_BETWEEN]); // k51
	if (n_events == 0) return;
  
	const int in = this->index[read_t];
	const float lTGF = this->patchChemConcentration(TGF, in);
	const float lTNF = this->patchChemConcentration(TNF, in);
	const float lFGF = this->patchChemConcentration(FGF, in);
	const float lHA  = Agent::agentECMPtr[in].HA[read_t];
  
	const float ha_b = Fibroblast::haSynth[HA_BASELINE]                          // k52
					   - Fibroblast::haSynth[HA_MESH_EFFECT] * BMWorld::meshSize // k53
					   + Fibroblast::haSynth[HA_ELASTICITY_EFFECT] * BMWorld::E; // k54
  
	const float ha = ha_b * (Fibroblast::haSynth[HA_ACTIVATING_STRENGTH]         // k55
								 * log10(safe_log10_arg(1.f + lTGF + lTNF + lFGF))
							 + Fibroblast::haSynth[HA_ACTIVATING_OFFSET]);       // k56
  
	if (ha <= 0.f) return;
  
	for (int e = 0; e < n_events; e++) {
		/* With probability HA + k57: move, then produce HA. */
		if (Agent::rollDice(lHA + Fibroblast::haSynth[HA_MOVE_PROBABILITY])) {       // k57
		this->wiggle();
		this->depositHA(ha, 1);
		}
	
		/* With probability k58 + HA/k59: produce HA in the same patch. */
		if (Agent::rollDice(Fibroblast::haSynth[HA_SAME_PATCH_PROBABILITY]           // k58
							+ lHA / safe_denominator(Fibroblast::haSynth[HA_SAME_PATCH_EFFECT]))) { // k59
		this->depositHA(ha, 1);
		}
	}
}

