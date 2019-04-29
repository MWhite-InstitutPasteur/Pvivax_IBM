/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///                                                                       ///
///  Individual-based Plasmodium vivax transmission model.                ///
///                                                                       ///
///  Dr Michael White                                                     ///
///  Institut Pasteur                                                     ///
///  michael.white@pasteur.fr                                             ///
///                                                                       ///
///  This file and accompanying .cpp contain:                             ///
///                                                                       ///
///     A structure called Population stores all individuals.             ///
///                                                                       ///
///  3. EQUILIBRIUM SETUP                                                 ///
///     This set of functions calculates the equilibrium set up of the    ///
///     population. It is only called once while the population is        ///
///     being initialised.                                                ///
///                                                                       ///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Include guard
#ifndef PVIVAX_MODEL_POPULATION
#define PVIVAX_MODEL_POPULATION

#include "Individual.hpp"


/////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
// 0.3. Define a structure for a population of individuals                             //
//                                                                                     //
//      Note that this object only stores information at a fixed point in time         //
//      and as such is memoryless                                                      //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

class Population
{
public:
    //////////////////////////////////////////////////////////////////////////
    //  Class member functions
    //////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////
    // Initialise a population of individuals at equilibrium
    void equi_pop_setup(Params& theta);

    ////////////////////////////////////////////////////
    // Update human individuals
    void human_step(Params& theta);
    
    ////////////////////////////////////////////////////
    // Summarise population outputs
    void summary();


private:
    void gauher(Params& theta);

public:
    //////////////////////////////////////////////////////////////////////////
    //  Data
    //////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////
    // 0.3.1. Human population

    int N_pop;                      // Population size - we have balanced demography at the moment so this will be effectively constant 

    vector<Individual> people;      // A vector of individuals

    vector<vector<double>> pi_n;    // Proportion of bites on humans that person n receives
    vector<vector<double>> lam_n;   // Biting rate of a single mosquito on person n


    ///////////////////////////////////////////////
    // 0.3.2. Mosquito population
    //
    //        Depends dynamically on vector control

    double yM[N_spec][N_M_comp];      // mosquito state

    double SUM_pi_w[N_spec];
    double SUM_pi_z[N_spec];


    double delta_1_VC[N_spec];        // duration spent foraging
    double delta_VC[N_spec];          // duration of gonotrophic cycle = duration between blood meals

    double Z_VC[N_spec];              // average probability of mosquito repeating during a single attempt
    double W_VC[N_spec];              // average probability of successfully feeding on a human during a single attempt
    double Q_VC[N_spec];              // human blood index
    double p_1_VC[N_spec];            // probability of surviving foraging
    double mu_M_VC[N_spec];           // mosquito death rate
    double aa_VC[N_spec];             // mosquito biting rate
    double exp_muM_tauM_VC[N_spec];   // probability of surviving sporogony
    double beta_VC[N_spec];           // egg oviposition rate


    ////////////////////////////////////////////////////////////////
    // 0.3.3. Objects for storing summary output of the population

    int yH[N_H_comp];   // Human compartmental states

    int prev_all[11];   // Contains {N_pop, PvPR_PCR, PvPR_LM, Pv_clin, PvHR, PvHR_batches, new_PCR, new_LM, new_D, new_ACT, new_PQ} 
    int prev_U5[11];    // Contains {N_pop, PvPR_PCR, PvPR_LM, Pv_clin, PvHR, PvHR_batches, new_PCR, new_LM, new_D, new_ACT, new_PQ} 
    int prev_U10[11];   // Contains {N_pop, PvPR_PCR, PvPR_LM, Pv_clin, PvHR, PvHR_batches, new_PCR, new_LM, new_D, new_ACT, new_PQ} 

    double EIR_t;       // EIR
    int LLIN_cov_t;     // LLIN coverage
    int IRS_cov_t;      // IRS coverage
    int ACT_treat_t;    // Coverage with front-line treatment (ACT) 
    int PQ_treat_t;     // Coverage with front-line treatment (primaquine or tafenoquine)
    int pregnant_t;     // Coverage with front-line treatment (primaquine or tafenoquine)

    int PQ_overtreat_t;
    int PQ_overtreat_9m_t;

    double A_par_mean_t;     // Average anti-parasite immunity 
    double A_clin_mean_t;    // Average anti-clinical immunity


    ///////////////////////////////////////////////////////////////
    // 0.3.4. Equilibrium settings of population
    //        These are only needed for initialising the simulation

    ////////////////////////////////////////
    // Age and heterogeneity demographics 

    double age_bounds[N_age + 1];

    double age_demog[N_age];
    double age_bite[N_age];
    double age_mids[N_age];

    double x_het_bounds[N_het + 1];

    double x_het[N_het];
    double w_het[N_het];

    double x_age_het[N_age][N_het];
    double w_age_het[N_age][N_het];


    ////////////////////////////////////////
    // Ageing rates 

    double r_age[N_age];

    double P_age_bite;     // Proportional change for age-dependent biting


    ////////////////////////////////////////
    // Indicator for age group of 20 year old woman

    int index_age_20;
};

#endif
