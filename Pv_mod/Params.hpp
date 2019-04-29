/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///                                                                       ///
///  Individual-based Plasmodium vivax transmission model.                ///
///                                                                       ///
///  Dr Michael White                                                     ///
///  Institut Pasteur                                                     ///
///  michael.white@pasteur.fr                                             ///
///                                                                       ///
///  This contains the Params struct and parameter reading code.          ///
///                                                                       ///
///  This header is common to all model code.                             ///
///                                                                       ///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Include guard
#ifndef PVIVAX_MODEL_PARAMS
#define PVIVAX_MODEL_PARAMS

#include <vector>

using namespace std;


///////////////////////////////////////////////////////////////////
//                                                               //
// 0.0. Define model constants                                   //
//                                                               //
///////////////////////////////////////////////////////////////////


#define t_step 1            // Time step for update in humans
#define mosq_steps 20       // Number of mosquito steps per human step

#define N_H_comp 6          // Number of human compartments (indexed by p)
#define N_M_comp 6          // Number of mossquito compartments (indexed by p)

#define N_age 58            // Number of age categories for calculation of equilibrium set up (indexed by i)
#define N_het 21            // Number of heterogeneity categories for calculation of equilibrium set up (indexed by j)
#define K_max 30            // Maximum umber of hypnozoites (indexed by k)
#define N_int 6             // Number of interventions

#define N_spec_max 3        // Maximum number of mosquito species that can be modelled
#define N_spec 3            // Number of mosquito species to be modelled

const double CONST_LOG_2 = 0.6931471805599453094172321214581766L;


///////////////////////////////////////////////////////////////////
//                                                               //
// 0.1. Define structure for parameters                          //
//                                                               //
///////////////////////////////////////////////////////////////////
struct Params
{
    //////////////////////////////////////////////////////////////////////////
    //  Functions
    //////////////////////////////////////////////////////////////////////////
    
    
    /////////////////////////////////////
    // Read parameters from input files
    void read(const char *parameter_File, const char *mosquito_File[N_spec_max]);


    //////////////////////////////////////////////////////////////////////////
    //  Data
    //////////////////////////////////////////////////////////////////////////

    int N_pop;                     // human population size
    double time_start, time_end, burnin_time;
    
    
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

    double BS_treat_BScover_base;  // proportion of episodes of symptomatic disease treated (baseline)
    double BS_treat_BSeff_base;    // proportion of episodes of symptomatic disease treated (baseline)
    double BS_treat_BSproph_base;  // proportion of episodes of symptomatic disease treated (baseline)

    double treat_BScover;          // proportion of episodes of symptomatic disease treated (changing)
    double treat_BSeff;            // blood-stage treatment efficacy
    double treat_PQavail;          // availability of PQ


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


    ///////////////////////////////
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


    ////////////////////////////////
    // Treatment parameters (treatment as intervention)

    double BS_treat_BScover;        // Coverage of first-line treatment with BS drugs
    double BS_treat_BSeff;          // Efficacy of first-line treatment with BS drugs
    double BS_treat_BSproph;        // Duration of prophylaxis with first-line BS drugs

    double PQ_treat_BScover;        // Coverage of BS treatment in combined BS & PQ first-line regimen
    double PQ_treat_BSeff;          // Efficacy of BS treatment
    double PQ_treat_BSproph;        // Duration of BS prophylaxis 
    double PQ_treat_PQavail;        // Availability of PQ treatment (as a proportion of those receiving BS treatment)
    double PQ_treat_PQeff;          // Efficacy of PQ treatment
    double PQ_treat_PQproph;        // Duration of PQ prophylaxis (i.e. for how long does PQ prevent new hypnozoites)
    int PQ_treat_G6PD_risk;         // Risk in G6PD - deficient individuals
    int PQ_treat_CYP2D6_risk;       // Risk of not working in low CYP2D6 metabolizers
    int PQ_treat_preg_risk;         // Risk in pregnant women
    double PQ_treat_low_age;        // Lower age limit for treatment(in years)

    double MDA_BS_BScover;          // Coverage of blood - stage drugs
    double MDA_BS_BSeff;            // Efficacy of blood - stage drugs
    double MDA_BS_BSproph;          // Duration of blood - stage prophylaxis

    double MDA_PQ_BScover;          // Coverage of blood - stage drugs
    double MDA_PQ_BSeff;            // Efficacy of blood - stage drugs
    double MDA_PQ_BSproph;          // Duration of blood - stage prophylaxis
    double MDA_PQ_PQavail;          // Availability of primaquine
    double MDA_PQ_PQeff;            // Efficacy of primaquine
    double MDA_PQ_PQproph;          // Duration of primaquine prophylaxis(prevents new hypnozoites)
    int MDA_PQ_G6PD_risk;           // Risk in G6PD - deficient individuals
    int MDA_PQ_CYP2D6_risk;         // Risk of not working in low CYP2D6 metabolizers
    int MDA_PQ_preg_risk;           // Risk in pregnant women
    double MDA_PQ_low_age;          // Lower age limit for treatment(in years)

    double MSAT_PQ_RDT_PCR;         // What diagnostic tool : 1 for RDT(= LM); 2 for PCR
    double MSAT_PQ_sens;            // Sensitivity of diagnostic tool
    double MSAT_PQ_BScover;         // Coverage of blood - stage drugs
    double MSAT_PQ_BSeff;           // Efficacy of blood - stage drugs
    double MSAT_PQ_BSproph;         // Duration of blood - stage prophylaxis
    double MSAT_PQ_PQavail;         // Availability of primaquine
    double MSAT_PQ_PQeff;           // Efficacy of primaquine
    double MSAT_PQ_PQproph;         // Duration of primaquine prophylaxis(prevents new hypnozoites)
    int MSAT_PQ_G6PD_risk;          // Risk in G6PD - deficient individuals
    int MSAT_PQ_CYP2D6_risk;        // Risk of not working in low CYP2D6 metabolizers
    int MSAT_PQ_preg_risk;          // Risk in pregnant women
    double MSAT_PQ_low_age;         // Lower age limit for treatment(in years)

    double SSAT_PQ_sens;            // Sensitivity of diagnostic tool
    double SSAT_PQ_spec;            // Sensitivity of diagnostic tool
    double SSAT_PQ_BScover;         // Coverage of blood - stage drugs
    double SSAT_PQ_BSeff;           // Efficacy of blood - stage drugs
    double SSAT_PQ_BSproph;         // Duration of blood - stage prophylaxis
    double SSAT_PQ_PQavail;         // Availability of primaquine
    double SSAT_PQ_PQeff;           // Efficacy of primaquine
    double SSAT_PQ_PQproph;         // Duration of primaquine prophylaxis(prevents new hypnozoites)
    int SSAT_PQ_G6PD_risk;          // Risk in G6PD - deficient individuals
    int SSAT_PQ_CYP2D6_risk;        // Risk of not working in low CYP2D6 metabolizers
    int SSAT_PQ_preg_risk;          // Risk in pregnant women
    double SSAT_PQ_low_age;         // Lower age limit for treatment(in years)


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


    ////////////////////////////////////////
    // Matrices for hypnozoite transitions

    double D_MAT[K_max + 1][K_max + 1];
    double OD_MAT[K_max + 1][K_max + 1];
    double K_MAT[K_max + 1][K_max + 1];
    double L_MAT[K_max + 1][K_max + 1];
    double H_MAT[K_max + 1][K_max + 1];
};

#endif
