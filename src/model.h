/*
 * model.h
 *
 *  Created on: 6 Apr 2020
 *      Author: gc1610
 */

#ifndef SRC_MODEL_H_
#define SRC_MODEL_H_

#include <vector>
#include <Rcpp.h>

using namespace std;

#define t_step 1            // Time step for update in humans
#define mosq_steps 20       // Number of mosquito steps per human step

#define N_H_comp 6          // Number of human compartments (indexed by p)
#define N_M_comp 6          // Number of mossquito compartments (indexed by p)

#define N_age 58            // Number of age categories for calculation of equilibrium set up (indexed by i)
#define N_het 21            // Number of heterogeneity categories for calculation of equilibrium set up (indexed by j)
#define K_max 30            // Maximum umber of hypnozoites (indexed by k)
#define N_int 6             // Number of interventions

#define N_spec_max 3
#define N_spec 3

#define log2 0.69314718056

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
//          //                                      //
//   ####   //   ###  ##### ######   ##  ## #####   //
//  ##  ##  //  ##    ##      ##     ##  ## ##  ##  //
//  ##  ##  //   ###  ####    ##     ##  ## #####   //
//  ##  ##  //     ## ##      ##     ##  ## ##      //
//   ####   //   ###  #####   ##      ####  ##      //
//          //                                      //
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////
//                                                               //
// 0.1. Define structure for parameters                          //
//                                                               //
///////////////////////////////////////////////////////////////////

struct params
{
    /////////////////////////////////////
    // Equilibrium EIR (no seasonality)

    double EIR_equil;              // EIR at equilibrium


    /////////////////////////////////
    // Human demography

    double mu_H;                   // human death rate
    double age_mean;               // mean age of human population
    double age_max;                // maximum age in human population
    double het_max;                // maximum heterogeneity

    double G6PD_prev;              // Prevalence of G6PD deficiency
    double CYP2D6_prev;            // Prevalence of CYP2D6 phenotype


    ////////////////////////////////////////////////////////
    // Heterogeneity in exposure and age-dependent biting

    double age_0;                  // age-dependent biting parameter
    double rho_age;                // age-dependent biting parameter

    double sig_het;                // heterogeneity in exposure - standard deviation on a log scale


    /////////////////////////////////
    // Transmission probabilities

    double bb;                     // mosquito -> human transmission probability
    double c_PCR;                  // human -> mosquito transmission probability - PCR detectable
    double c_LM;                   // human -> mosquito transmission probability - LM detectable
    double c_D;                    // human -> mosquito transmission probability - clinical disease
    double c_T;                    // human -> mosquito transmission probability - treatment

    double d_latent;               // duration of latency in the liver
    int H_track;                   // Number of steps for tracking lagged lam_H (due to duration of latency)


    /////////////////////////////////
    // Human recovery parameters

    double r_LM;                   // rate of recovery from LM detectable infection
    double r_D;                    // rate of recovery from symptomatic disease
    double r_T;                    // rate of progression through treatment
    double r_P;                    // rate of loss of prophylaxis

    double d_PCR_min;              // minimum duration of PCR-detectable infection - full immunity
    double d_PCR_max;              // maximum duration of PCR-detectable infection - no immunity
    double A_PCR_50pc;             // scale parameter for effect of anti-parasite immunity on PCR-detectable infection
    double K_PCR;                  // shape parameter for effect of anti-parasite immunity on PCR-detectable infection

    double r_PCR;                  // recovery rate from PCR detectable infection - depends on level of immunity


    /////////////////////////////////
    // Baseline treatment parameters

    double BS_treat_cov_base;      // proportion of episodes of symptomatic disease treated (baseline)
    double BS_treat_eff_base;      // proportion of episodes of symptomatic disease treated (baseline)
    double BS_treat_BSproph_base;  // proportion of episodes of symptomatic disease treated (baseline)

    double treat_cov;              // proportion of episodes of symptomatic disease treated (changing)
    double treat_eff;              // blood-stage treatment efficacy


    ////////////////////////////////
    // Anti-parasite immunity parameters

    double u_par;                  // scale paramter for acquisition of blood-stage immunity
    double r_par;                  // rate of decay of blood-stage immunity

    double phi_LM_max;             // probability of LM-detectable infection with no immunity
    double phi_LM_min;             // probability of LM-detectable infection with full immunity
    double A_LM_50pc;              // blood-stage immunity scale parameter
    double K_LM;                   // blood-stage immunity shape parameter

    double phi_LM;                 // probability that PCR detectable infection becomes detectable by light microscopy (LM) - DYNAMICALLY UPDATED


    ////////////////////////////////
    // Clinical immunity parameters

    double u_clin;                 // scale paramter for acquisition of blood-stage immunity
    double r_clin;                 // rate of decay of clinical immunity

    double phi_D_max;              // probability of clinical episode with no immunity
    double phi_D_min;              // probability of clinical episode with full immunity
    double A_D_50pc;               // clinical immunity scale parameter
    double K_D;                    // clinical immunity shape parameter

    double phi_D;                  // probability that LM detectable infection progresses to symptomatic disease - DYNAMICALLY UPDATED


    /////////////////////////////////////////////
    // maternal im munity

    double P_mat;                  // New-born immunity relative to mother's
    double d_mat;                  // Inverse of decay rate of maternal immunity


    /////////////////////////////////
    // Relapse paramters

    double ff;                     // relapse rate
    double gamma_L;                // liver clearance rate


    /////////////////////////////////
    // Entomological paramters

    double Prop_mosq[N_spec_max];      // Proportion of An. farauti

    double mm_0[N_spec];               // number of mosquitoes per human (An. farauti)
    double aa[N_spec];                 // mosquito biting rate (in the absence of vector control)
    double mu_M[N_spec];               // mosquito death rate
    double tau_M[N_spec];              // duration of sporogony

    double lam_M[N_spec];              // Force of infection on mosquites - updated dynamically

    int M_track;                           // Number of steps for tracking lagged lam_M*S_M (needed for lag due to duration of sporogony)
    vector<vector<double>> lam_S_M_track;  // Lagged force of infection on moquitoes


    ////////////////////////////////
    // Seasonality paramters

    double dry_seas[N_spec];        // Proportion of dry season transmission compared to mean
    double kappa_seas[N_spec];      // Shape parameter for seasonality
    double t_peak_seas[N_spec];     // Offset for seasonal transmission
    double denom_seas[N_spec];      // Denominator for seasonality


    /////////////////////////////////
    // Larval paramters

    double d_E_larvae;              // Development time of early larval instars
    double d_L_larvae;              // Development time of late larval instars
    double d_pupae;                 // Development time of pupae
    double mu_E0;                   // Mortality rate of early larval instars (low density)
    double mu_L0;                   // Mortality rate of late larval instars (low density)
    double mu_P;                    // Mortality rate of pupae
    double beta_larvae;             // Number of eggs laid per day per mosquito
    double gamma_larvae;            // Effect of density dependence on late instars relative to early instars

    double omega_larvae[N_spec];    // Useful pre-calculated quantity

    double Karry[N_spec];           // Larval carry capacity

    double eps_max[N_spec];         // Number of eggs per day


    /////////////////////////////////
    // Treatment parameters (treatment as intervention)

    double MDA_BS_cover;            // Coverage of MDA with BS drugs
    double MDA_BS_BSeff;            // Efficacy of MDA with BS drugs
    double MDA_BS_BSproph;          // Duration of prophylaxis of BS drugs used for MDA

    double MDA_PQ_cover;            // Coverage of MDA with BS drugs and PQ
    double MDA_PQ_BSeff;            // Efficacy of MDA with BS drugs and PQ
    double MDA_PQ_PQeff;            // Efficacy of MDA with BS drugs and PQ
    double MDA_PQ_BSproph;          // Duration of prophylaxis of PQ used for MDA
    double MDA_PQ_PQproph;          // Duration of prophylaxis of PQ used for MDA
    int MDA_PQ_CYP2D6;              // Indicator for whether low CYP2D6 metabolization is a problem (yes=1 e.g. PQ; no=0, e.g. TQ)


    double BS_treat_cover;          // Coverage of first-line treatment with BS drugs
    double BS_treat_BSeff;          // Efficacy of first-line treatment with BS drugs
    double BS_treat_BSproph;        // Duration of prophylaxis with first-line BS drugs

    double PQ_treat_cover;          // Coverage of PQ treatment (as a proportion of those receiving BS treatment)
    double PQ_treat_PQcover;        // Coverage of PQ treatment (as a proportion of those receiving BS treatment)
    double PQ_treat_BSeff;          // Efficacy of PQ treatment
    double PQ_treat_PQeff;          // Efficacy of PQ treatment
    double PQ_treat_BSproph;        // Duration of PQ prophylaxis (i.e. for how long does PQ prevent new hypnozoites)
    double PQ_treat_PQproph;        // Duration of PQ prophylaxis (i.e. for how long does PQ prevent new hypnozoites)
    int PQ_treat_CYP2D6;            // Indicator for whether low CYP2D6 metabolization is a problem (yes=1 e.g. PQ; no=0, e.g. TQ)


    /////////////////////////////////////////////////////////
    // Baseline entomological parameters required for LLINs

    double LLIN_half_life;      // Half-life of loss of LLINs
    double P_LLIN_loss;         // Daily probability of losing LLIN - pre-calculated for efficiency
    double PYR_half_life;       // Half-life of pyrethroid decay on LLINs
    double PYR_decay;           // Daily pyrethroid decay - pre-calculated for efficiency

    double r_LLIN_0[N_spec];    // Probability mosquito repelled (with full insecticide activity)
    double r_LLIN_net[N_spec];  // Probability mosquito repelled due to barrier effect of net (no insecticide)
    double s_LLIN_0[N_spec];    // Probability mosquito feeds successfully
    double d_LLIN_0[N_spec];    // Probability mosquito dies during feeding attempt

    double IRS_half_life;       // Half-life of IRS insecticide decay
    double IRS_decay;           // Daily IRS insecticide decay - pre-calculated for efficiency

    double r_IRS_0[N_spec];     // Probability mosquito repelled by IRS (full insecticide activity)
    double d_IRS_0[N_spec];     // Probability mosquito killed by IRS (full insecticide activity)
    double s_IRS_0[N_spec];     // Probability mosquito survives feeding attempt with IRS (full insecticide activity)

    double Q_0[N_spec];         // Human Blood Index (proportion of blood meals taken on humans)
    double CHI_endo[N_spec];    // Endophily - proportion of mosquitoes resting indoors after feeding (no intervention)
    double PSI_indoors[N_spec]; // Proportion of bites taken on humans indoors
    double PSI_bed[N_spec];     // Proportion of bites taken on humans in bed

    double delta_1;             // Time spent foraging for a blood meal
    double delta_2;             // Time_spent_digesting_blood_meal
    double delta;               // Duration of gonotrophic cycle

    double p_1[N_spec];         // Daily? death probability when foraging for a blood meal
    double p_2[N_spec];         // Daily? death probability when digesting blood meal

    double rho_round_LLIN;      // Between round correlation of LLINS
    double rho_round_IRS;       // Between round correlation of IRS
    double rho_LLIN_IRS;        // Correlation between LLIN and IRS coverage

    double rho_round_MDA;       // Between round correlation of MDA
    double rho_MDA_LLIN;        // Correlation between MDA and LLINs
    double rho_MDA_IRS;         // Correlation between MDA and IRS

    double sig_round_LLIN;      // Derived parameter for correlation between rounds of LLINs
    double sig_round_IRS;       // Derived parameter for correlation between rounds of IRS
    double sig_round_MDA;       // Derived parameter for correlation between rounds of MDA

    double SSAT_sens;           // Sensitivity of serological screen and treat
    double SSAT_spec;           // Specificity of serological screen and treat

    float V_int[N_int][N_int];
    float V_int_dummy[N_int][N_int];


    //////////////////////////////////////////////////////
    // Pre-multiplication of quantities for efficiency

    double A_par_decay;         // Decay factor for BS immunity
    double A_clin_decay;        // Decay factor for clinical immunity
    double mat_decay;           // Decay factor for maternal immunity

    double age_0_inv;           // Inverse of age-dependent biting parameter

    double A_PCR_50pc_inv;      // Immune scalar for PCR-detectable infection
    double A_LM_50pc_inv;       // Immune scalar for LM-detectable infection
    double A_D_50pc_inv;        // Immune scalar for clinical disease

    double P_dead;              // Probability of dying in each time step
    double P_preg;              // Probability of woman 18-40 years becoming pregnant per day

    double P_PYR_decay;         // Proportional decay of pyrethroid per time step (pre-calculated for efficiency)
    double P_IRS_decay;         // Proportional decay of IRS insecticide per time step (pre-calculated for efficiency)


    ////////////////////////////////////
    // Probability of exiting states

    double S_out;
    double I_PCR_out;
    double I_LM_out;
    double I_D_out;
    double T_out;
    double P_out;


    ///////////////////////////////////////////////////
    // Vector of probabilities for competing hazards
    //
    // These are stored in the parameter structure for convenience.
    // It would perhaps be more natural to have them specific for each
    // individual, but that would require storing them N_pop times.

    double S_move[4];
    double I_PCR_move[5];
    double I_LM_move[4];
    double I_D_move[2];
    double T_move[2];
    double P_move[2];


    ////////////////////////////////////////
    // Matrices for hypnozoite transitions

    double D_MAT[K_max + 1][K_max + 1];
    double OD_MAT[K_max + 1][K_max + 1];
    double K_MAT[K_max + 1][K_max + 1];
    double L_MAT[K_max + 1][K_max + 1];
    double H_MAT[K_max + 1][K_max + 1];
};


/////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
// 0.2. Define a class for humans                                                      //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

class individual
{
public:

    ////////////////////////////////////////////////////
    // 0.2.1. Class constructor

    individual(double a, double zeta)
    {
        age = a;
        zeta_het = zeta;
    }


    ////////////////////////////////////////////////////
    // 0.2.2. Function declarations within the human class

    void state_mover(params theta, double lam_bite);
    void ager(params theta);
    void intervention_updater(params theta);


    ////////////////////////////////////////////////////
    // 0.2.3. Person-specific age and exposure classifiers

    double age;                      // Person's age
    double zeta_het;                 // Heterogeneity in exposure to mosquitoes

    bool gender;                     // 0 = male; 1 = female
    bool G6PD_def;                   // Is the person G6PD deficient? 0 = normal; 1 = deficient (homozygous); 2 = deficient (heterozygous - female only)
    bool CYP2D6;                     // Does the person have CYP2D6 phenotype? 0 = normal; 1 = low metabolizer


    ////////////////////////////////////////////////////
    //  0.2.4. Child-bearing age.
    //         Indicator for people of 18-22 years of age. New-born
    //         children acquire a proportion of the immunity of a
    //         women aged 20 years.

    bool preg_age;
    bool pregnant;
    double preg_timer;


    double lam_bite_lag;             // Lagged force of infection due to moquito bites
    vector<double> lam_bite_track;   // Tracking lagged force of infection due to moquito bites

    double lam_rel_lag;              // Lagged force of infection due to relapses
    vector<double> lam_rel_track;    // Tracking lagged force of infection due to relapses

    double lam_H_lag;                // Lagged total force of infection


    ////////////////////////////////////////////////////
    // 0.2.5. Indicators for compartments

    bool S;
    bool I_PCR;
    bool I_LM;
    bool I_D;
    bool T;
    bool P;

    /////////////////////////////////////////////////////////////////
    //  0.2.7. Number of batches of hypnozoites. Must be an integer.

    int Hyp;


    ////////////////////////////////////////////////////
    // Indicator for competing hazards move

    int CH_move;


    ////////////////////////////////////////////////////
    //  0.2.6. Indicators for new events

    bool I_PCR_new;    // New PCR-detectable infection (I_LM, I_D & T included here)
    bool I_LM_new;     // New LM-detectable infection (I_D & T included here)
    bool I_D_new;      // New clinical episode (treated or untreated)
    bool ACT_new;      //
    bool PQ_new;       //


    ////////////////////////////////////////////////////
    //  0.2.8. Person-specific levels of immunity and indicators
    //         for whether immunity is suppressed

    double A_par;
    double A_clin;

    double A_par_mat;
    double A_clin_mat;

    bool A_par_boost;
    bool A_clin_boost;

    double A_par_timer;
    double A_clin_timer;

    bool PQ_proph;
    double PQ_proph_timer;


    //////////////////////////////////////////////////////////
    // 0.2.9. Person-specific intervention access parameter

    double zz_int[N_int];


    ///////////////////////////////
    // LLINs

    bool LLIN;                 // Does the person have an LLIN?
    double LLIN_age;           // How old is the LLIN

    double r_LLIN[N_spec];     // repellency
    double d_LLIN[N_spec];     // killing effect
    double s_LLIN[N_spec];     // survival


    ///////////////////////////////
    // IRS

    bool IRS;                 // Is the person protected by IRS
    double IRS_age;           // How long ago was their house sprayed?

    double r_IRS[N_spec];     // repellency
    double d_IRS[N_spec];     // killing effect
    double s_IRS[N_spec];     // survival


    ///////////////////////////////////////////////////
    // Individual-level effect of vector control

    double z_VC[N_spec];      // probability of mosquito being repelled from this individual during a single feeding attempt
    double y_VC[N_spec];      // probability of mosquito feeding on this individual during a single attempt
    double w_VC[N_spec];      // probability of mosquito feeding and surviving on this individual during a single feeding attempt

private:
};


/////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
// 0.3. Define a structure for a population of individuals                             //
//                                                                                     //
//      Note that this object only stores information at a fixed point in time         //
//      and as such is memoryless                                                      //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

struct population
{
    ////////////////////////////////////////////////
    // 0.3.1. Human population

    int N_pop;                      // Population size - we have balanced demography at the moment so this will be effectively constant

    vector<individual> people;      // A vector of individuals

    vector<vector<double>> pi_n;    // Proportion of bites on humans that person n receives
    vector<vector<double>> lam_n;   // Biting rate of a single mosquito on person n


    ////////////////////////////////////////////////
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

    std::vector<std::pair<int, int>> prev_groups; // Contains age groups for prevalence summaries
    std::vector<std::pair<int, int>> incidence_groups; // Contains age groups for incidence summaries
    std::vector<std::vector<int>> prev_summaries;// Contains N_pop, PvPR_PCR, PvPR_LM, Pv_clin, PvHR, PvHR_batches,
    std::vector<std::vector<int>> incidence_summaries; // Contains new_PCR, new_LM, new_D, new_T


    double EIR_t;       // EIR
    int LLIN_cov_t;     // LLIN coverage
    int IRS_cov_t;      // IRS coverage
    int ACT_treat_t;    // Coverage with front-line treatment (ACT)
    int PQ_treat_t;     // Coverage with front-line treatment (primaquine or tafenoquine)
    int pregnant_t;     // Coverage with front-line treatment (primaquine or tafenoquine)

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


/////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
// 0.4. Define a structure for details of interventions                                //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

struct intervention
{
    ////////////////////////////////
    // LLINs

    vector<double> LLIN_year;
    vector<double> LLIN_cover;


    ////////////////////////////////
    // IRS

    vector<double> IRS_year;
    vector<double> IRS_cover;


    ////////////////////////////////
    // MDA - blood-stage drugs

    vector<double> MDA_BS_year;
    vector<double> MDA_BS_cover;
    vector<double> MDA_BS_BSeff;
    vector<double> MDA_BS_BSproph;


    ////////////////////////////////
    // MDA - blood-stage drugs + primaquine

    vector<double> MDA_PQ_year;
    vector<double> MDA_PQ_cover;
    vector<double> MDA_PQ_BSeff;
    vector<double> MDA_PQ_PQeff;
    vector<double> MDA_PQ_BSproph;
    vector<double> MDA_PQ_PQproph;
    vector<int>    MDA_PQ_CYP2D6;


    ////////////////////////////////
    // First-line treatment - blood-stage drugs

    vector<double> BS_treat_year_on;
    vector<double> BS_treat_year_off;
    vector<double> BS_treat_cover;
    vector<double> BS_treat_BSeff;
    vector<double> BS_treat_BSproph;


    ////////////////////////////////
    // First-line treatment - blood-stage drugs
    // plus primaquine

    vector<double> PQ_treat_year_on;
    vector<double> PQ_treat_year_off;
    vector<double> PQ_treat_cover;
    vector<double> PQ_treat_PQcover;
    vector<double> PQ_treat_BSeff;
    vector<double> PQ_treat_PQeff;
    vector<double> PQ_treat_BSproph;
    vector<double> PQ_treat_PQproph;
    vector<int>    PQ_treat_CYP2D6;
};


/////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
// 0.5. Define a structure for storing the output of a simulation                      //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

struct simulation
{
    //////////////////////////////////////////
    // 0.5.1. Vector of simulation times

    int N_time;

    vector<double> t_vec;


    //////////////////////////////////////////
    // 0.5.2. Tracking output

    vector<vector<int>> yH_t;
    vector<vector<vector<double>>> yM_t;


    vector<double> EIR_t;

    std::vector<std::pair<int, int>> prev_groups; // Contains age groups for prevalence summaries
    std::vector<std::pair<int, int>> incidence_groups; // Contains age groups for incidence summaries
    std::vector<std::vector<std::vector<int>>> prev_summaries;// Contains N_pop, PvPR_PCR, PvPR_LM, Pv_clin, PvHR, PvHR_batches,
    std::vector<std::vector<std::vector<int>>> incidence_summaries; // Contains new_PCR, new_LM, new_D, new_T

    /////////////////////////////////////////
    // 0.5.3. Tracking coverage over time

    vector<int> LLIN_cov_t;
    vector<int> IRS_cov_t;
    vector<int> ACT_treat_t;
    vector<int> PQ_treat_t;
    vector<int> pregnant_t;


    //////////////////////////////////////////
    // 0.5.4. Tracking immunity over time

    vector<double> A_par_mean_t;
    vector<double> A_clin_mean_t;
};


////////////////////////////////////////////////////////////
//                                                        //
// 0.6. Function declarations                             //
//                                                        //
////////////////////////////////////////////////////////////

void mosq_derivs(const double t, double(&yM)[N_spec][N_M_comp], double(&dyMdt)[N_spec][N_M_comp], params* theta, population* POP);
void mosq_rk4(const double t, const double t_step_mosq, double(&yM)[N_spec][N_M_comp], params* theta, population* POP);
void mosquito_step(double t, params* theta, population* POP);
void human_step(params* theta, population* POP);
void intervention_dist(double t, params* theta, population* POP, intervention* INTVEN);
void POP_summary(population* POP, simulation* SIM);
void model_simulator(params* theta, population* POP, intervention* INTVEN, simulation* SIM);
int CH_sample(double *xx, int nn);
double phi_inv(double pp, double mu, double sigma);
double gammln(const double xx);
Rcpp::DataFrame run_simulation(
    const char*,
    const char**,
    const char*,
    std::vector<int>,
    std::vector<int>,
    std::vector<int>,
    std::vector<int>
);
Rcpp::DataFrame run_simulation_from_path(
    string,
    string,
    string,
    string,
    string,
    std::vector<int>,
    std::vector<int>,
    std::vector<int>,
    std::vector<int>
);


///////////////////////////////////////////
// Functions only needed at start for setting up
// the population at equilibrium

void ludcmp(vector<vector<double>> &a, int n_dim, vector<int> &indx, double &d);
void lubksb(vector<vector<double>> &a, int n_dim, vector<int> &indx, vector<double> &b);
void matrix_inv(vector<vector<double>> &a, int n, vector<vector<double>> &a_inv);
void inv_MM_bb(vector<vector<double>> &MM, vector<double> &bb, vector<double> &xx, int n_dim);
void MM_ij(int i, int j, params* theta, population* POP, vector<vector<double>> &MM,
           vector<vector<double>> lam_eq, vector<vector<vector<double>>> phi_LM_eq,
           vector<vector<vector<double>>> phi_D_eq, vector<vector<vector<double>>> r_PCR_eq);
void gauher(population* POP, params* theta);
void equi_pop_setup(population* POP, params* theta);

#endif /* SRC_MODEL_H_ */
