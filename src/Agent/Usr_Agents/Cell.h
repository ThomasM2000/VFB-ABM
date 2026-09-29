/*
 * File: Cell.h
 *
 * File Contents: Declarations for the Cell base class and the Fibroblast agent
 *                of the VFB-ABM.
 *
 * The manuscript models a single agent type - the vocal fold fibroblast - with
 * two phenotypes, unactivated and activated (Figure 1C). Activation is a state
 * of the agent (Agent::activate[]), not a separate class: only activated
 * fibroblasts synthesize cytokines and extracellular matrix, while both
 * phenotypes migrate, proliferate and die.
 */

 #ifndef CELL_H
 #define	CELL_H
 
 #include "../Agent.h"
 #include "../../Patch/Patch.h"
 #include "../../World/Usr_World/biomaterialWorld.h"
 #include "../biology_parameters_config.h"
 
 class ECM;
 
 #include <stdlib.h>
 #include <vector>
 #include <cmath>
 #include <omp.h>
 
 using namespace std;
 
 /*
  * CELL CLASS DESCRIPTION: Cell is a derived class of Agent. It provides the
  *                         per-tick template methods (migration, proliferation,
  *                         activation, cytokine synthesis, ECM synthesis, death)
  *                         and leaves the numerical rules to virtual hooks that
  *                         concrete agents override.
  */
 class Cell: public Agent {
   public:
     Cell();
     Cell(Patch* patchPtr);
     Cell(int x, int y, int z);
     ~Cell();
 
     /* Per-tick agent behaviour, ordered as in Figure 1B/1C. */
     void cellFunction();
 
     /** Migration: move along the chemotactic gradient at the speed given by
      *  the migration-speed hook (Table 3, rule 1). */
     virtual void cellSniff() final;
 
     /** Proliferation: Table 3, rule 3. */
     virtual void proliferate() final;
 
     /** Activation / deactivation: Table 3, rules 9 and 10. */
     virtual void activation() final;
 
     /** Cytokine synthesis: Table 3, rules 4-8 (activated fibroblasts only). */
     virtual void cytokine_synthesis() final;
 
     /** ECM synthesis: Table 3, rules 11-13 (activated fibroblasts only). */
     virtual void ecm_synthesis() final;
 
     /** Cell death driven by the viability rate (Table 3, rule 2). */
     virtual void apoptose() final;
 
     /** Performs cell death and frees the patch. */
     void die();
 
     void copyAndInitialize(Agent* original, int dx, int dy, int dz = 0);
 
     /** Deposit ECM protein of the given type on a random neighbouring patch. */
     void depositCollagen(float amount);
     void depositElastin(float amount);
     /** Deposit HA either here (here = 1) or on a neighbour (here = 0). */
     void depositHA(float amount, int here);
 
     /** Hatches `number` new cells of `agentType` on unoccupied neighbours. */
     virtual void hatchnewcell(int number, int agentType, int here = 0);
 
     /* ---------------------------------------------------------------------- */
     /*                            STATIC VARIABLES                            */
     /* ---------------------------------------------------------------------- */
     static int numOfCells;  // Number of living cells
 
   protected:
     /** Senescence cap: cells stop dividing after this many doublings. */
     virtual int get_max_doublings();
     /* ------------------------- Agent-rule hooks ------------------------- */
     /** Table 3 rule 1, v = -k1 ln(E) + k2, returned in patches/tick. */
     virtual float get_migration_speed();
     /** Table 3 rule 2, vr = k3 ln(t) + k4, a survival percentage. */
     virtual float get_viability_rate();
     /** Table 3 rule 3, returned as a percentage for rollDice(). */
     virtual float get_prolif_prob();
     /** Table 3 rules 9/10; return true to change activation state. */
     virtual bool should_activate(float patchTGF);
     virtual bool should_deactivate();
     /** Table 3 rules 4-8; writes into the per-tick secretion channels. */
     virtual void create_cytokines();
     /** Table 3 rules 11-13. */
     virtual void create_collagen();
     virtual void create_elastin();
     virtual void create_ha();

 };
 
 /*
  * FIBROBLAST CLASS DESCRIPTION: the vocal fold fibroblast of the VFB-ABM.
  *                               All numerical agent rules of Table 3 and their
  *                               parameters k1..k59 live here.
  */
 class Fibroblast: public Cell {
   public:
     Fibroblast();
     Fibroblast(Patch* patchPtr);
     Fibroblast(int x, int y, int z);
     ~Fibroblast();
 
   /* ---------------------------------------------------------------------- */
   /*                            STATIC VARIABLES                            */
   /* ---------------------------------------------------------------------- */
   static int numOfFibroblast;   // Living unactivated fibroblasts
   static int numOfAFibroblast;  // Living activated fibroblasts
 
   /* --------------- Agent-rule parameters (Manuscript Table 3) ------------ */
 
   /** k1, k2 - migration speed. */
   enum MigrationIdx {
     MIGRATION_ELASTICITY_EFFECT = 0,  // k1
     MIGRATION_BASELINE_SPEED,         // k2
     MIGRATION_COUNT
   };
   static_assert(sizeof(MigrationParams) / sizeof(double) == MIGRATION_COUNT,
                 "MigrationParams field count must match Fibroblast::MigrationIdx");
   static float migration[MIGRATION_COUNT];
 
   /** k3, k4 - viability rate. */
   enum ViabilityIdx {
     VIABILITY_TIME_EFFECT = 0,     // k3
     VIABILITY_BIOMATERIAL_EFFECT,  // k4
     VIABILITY_COUNT
   };
   static_assert(sizeof(ViabilityParams) / sizeof(double) == VIABILITY_COUNT,
                 "ViabilityParams field count must match Fibroblast::ViabilityIdx");
   static float viability[VIABILITY_COUNT];
 
   /** k5 .. k14 - proliferation. */
   enum ProliferationIdx {
     PROLIFERATION_TGF_THRESHOLD = 0,       // k5
     PROLIFERATION_HOURS_BETWEEN,           // k6
     PROLIFERATION_TIME_EFFECT,             // k7
     PROLIFERATION_HA_TIME_EFFECT,          // k8
     PROLIFERATION_HA_EFFECT,               // k9
     PROLIFERATION_BIOMATERIAL_BASELINE,    // k10
     PROLIFERATION_BASELINE,                // k11
     PROLIFERATION_ACTIVATING_FACTOR_EFFECT,// k12
     PROLIFERATION_BASELINE_OFFSET,         // k13
     PROLIFERATION_CHEMICAL_EFFECT,         // k14
     PROLIFERATION_COUNT
   };
   static_assert(sizeof(ProliferationParams) / sizeof(double) == PROLIFERATION_COUNT,
                 "ProliferationParams field count must match Fibroblast::ProliferationIdx");
   static float proliferation[PROLIFERATION_COUNT];
 
   /** k15 .. k17 - TGF-beta synthesis. */
   enum TgfSynthesisIdx {
     TGF_BASELINE = 0,       // k15
     TGF_FEEDBACK_EFFECT,    // k16
     TGF_FEEDBACK_BASELINE,  // k17
     TGF_SYNTHESIS_COUNT
   };
   static_assert(sizeof(TgfSynthesisParams) / sizeof(double) == TGF_SYNTHESIS_COUNT,
                 "TgfSynthesisParams field count must match Fibroblast::TgfSynthesisIdx");
   static float tgfSynthesis[TGF_SYNTHESIS_COUNT];
 
   /** k18 - FGF synthesis. */
   enum FgfSynthesisIdx { FGF_RATE = 0, FGF_SYNTHESIS_COUNT };  // k18
   static_assert(sizeof(FgfSynthesisParams) / sizeof(double) == FGF_SYNTHESIS_COUNT,
                 "FgfSynthesisParams field count must match Fibroblast::FgfSynthesisIdx");
   static float fgfSynthesis[FGF_SYNTHESIS_COUNT];
 
   /** k19 .. k21 - TNF-alpha synthesis. */
   enum TnfSynthesisIdx {
     TNF_BASELINE = 0,      // k19
     TNF_INHIBITORY_EFFECT, // k20
     TNF_IL10_EFFECT,       // k21
     TNF_SYNTHESIS_COUNT
   };
   static_assert(sizeof(TnfSynthesisParams) / sizeof(double) == TNF_SYNTHESIS_COUNT,
                 "TnfSynthesisParams field count must match Fibroblast::TnfSynthesisIdx");
   static float tnfSynthesis[TNF_SYNTHESIS_COUNT];
 
   /** k22 .. k24 - IL-6 synthesis. */
   enum Il6SynthesisIdx {
     IL6_BASELINE = 0,       // k22
     IL6_ACTIVATING_EFFECT,  // k23
     IL6_INHIBITORY_EFFECT,  // k24
     IL6_SYNTHESIS_COUNT
   };
   static_assert(sizeof(Il6SynthesisParams) / sizeof(double) == IL6_SYNTHESIS_COUNT,
                 "Il6SynthesisParams field count must match Fibroblast::Il6SynthesisIdx");
   static float il6Synthesis[IL6_SYNTHESIS_COUNT];
 
   /** k25, k26 - IL-8 synthesis. */
   enum Il8SynthesisIdx {
     IL8_BASELINE = 0,      // k25
     IL8_INHIBITORY_EFFECT, // k26
     IL8_SYNTHESIS_COUNT
   };
   static_assert(sizeof(Il8SynthesisParams) / sizeof(double) == IL8_SYNTHESIS_COUNT,
                 "Il8SynthesisParams field count must match Fibroblast::Il8SynthesisIdx");
   static float il8Synthesis[IL8_SYNTHESIS_COUNT];
 
   /** k27 .. k31 - activation and deactivation. */
   enum ActivationIdx {
     ACTIVATION_TGF_INTERMEDIATE_THRESHOLD = 0, // k27
     ACTIVATION_INTERMEDIATE_PROBABILITY,       // k28
     ACTIVATION_TGF_HIGH_THRESHOLD,             // k29
     ACTIVATION_INDEPENDENT_PROBABILITY,        // k30
     ACTIVATION_DEACTIVATION_PROBABILITY,       // k31
     ACTIVATION_COUNT
   };
   static_assert(sizeof(ActivationParams) / sizeof(double) == ACTIVATION_COUNT,
                 "ActivationParams field count must match Fibroblast::ActivationIdx");
   static float activation_params[ACTIVATION_COUNT];
 
   /** k32 .. k43 - collagen synthesis. */
   enum CollagenSynthIdx {
     COLLAGEN_HOURS_BETWEEN = 0,      // k32
     COLLAGEN_BASELINE,               // k33
     COLLAGEN_MESH_EFFECT,            // k34
     COLLAGEN_ELASTICITY_EFFECT,      // k35
     COLLAGEN_ACTIVATING_STRENGTH,    // k36
     COLLAGEN_ACTIVATING_OFFSET,      // k37
     COLLAGEN_IL8_EFFECT,             // k38
     COLLAGEN_HIGH_FRAG_RATIO_EFFECT, // k39
     COLLAGEN_EQUAL_FRAG_RATIO_EFFECT,// k40
     COLLAGEN_EQUAL_FRAG_OFFSET,      // k41
     COLLAGEN_BASELINE_RATE,          // k42
     COLLAGEN_CHEMICAL_EFFECT,        // k43
     COLLAGEN_SYNTH_COUNT
   };
   static_assert(sizeof(CollagenSynthesisParams) / sizeof(double) == COLLAGEN_SYNTH_COUNT,
                 "CollagenSynthesisParams field count must match Fibroblast::CollagenSynthIdx");
   static float collagenSynth[COLLAGEN_SYNTH_COUNT];
 
   /** k44 .. k50 - elastin synthesis (shares the k32 interval with collagen). */
   enum ElastinSynthIdx {
     ELASTIN_BASELINE = 0,         // k44
     ELASTIN_MESH_EFFECT,          // k45
     ELASTIN_ELASTICITY_EFFECT,    // k46
     ELASTIN_ACTIVATING_STRENGTH,  // k47
     ELASTIN_ACTIVATING_OFFSET,    // k48
     ELASTIN_INHIBITORY_STRENGTH,  // k49
     ELASTIN_BASELINE_RATE,        // k50
     ELASTIN_SYNTH_COUNT
   };
   static_assert(sizeof(ElastinSynthesisParams) / sizeof(double) == ELASTIN_SYNTH_COUNT,
                 "ElastinSynthesisParams field count must match Fibroblast::ElastinSynthIdx");
   static float elastinSynth[ELASTIN_SYNTH_COUNT];
 
   /** k51 .. k59 - hyaluronic acid synthesis. */
   enum HaSynthIdx {
     HA_HOURS_BETWEEN = 0,       // k51
     HA_BASELINE,                // k52
     HA_MESH_EFFECT,             // k53
     HA_ELASTICITY_EFFECT,       // k54
     HA_ACTIVATING_STRENGTH,     // k55
     HA_ACTIVATING_OFFSET,       // k56
     HA_MOVE_PROBABILITY,        // k57
     HA_SAME_PATCH_PROBABILITY,  // k58
     HA_SAME_PATCH_EFFECT,       // k59
     HA_SYNTH_COUNT
   };
   static_assert(sizeof(HaSynthesisParams) / sizeof(double) == HA_SYNTH_COUNT,
                 "HaSynthesisParams field count must match Fibroblast::HaSynthIdx");
   static float haSynth[HA_SYNTH_COUNT];
 
   protected:
     float get_migration_speed() override;
     float get_viability_rate() override;
     float get_prolif_prob() override;
     bool should_activate(float patchTGF) override;
     bool should_deactivate() override;
     void create_cytokines() override;
     void create_collagen() override;
     void create_elastin() override;
     void create_ha() override;
 };
 
 #endif
 