/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///                                                                       ///
///  Individual-based Plasmodium vivax transmission model.                ///
///                                                                       ///
///  Dr Michael White                                                     ///
///  Institut Pasteur                                                     ///
///  michael.white@pasteur.fr                                             ///
///                                                                       ///
///  With contributions from Dr Thomas Obadia                             /// 
///                                                                       ///
///  Please feel free to use and modify if you wish. However,             ///
///  please provide appropriate acknowledgement and get in touch          ///
///  if you have any questions. This is not necessarily the               ///
///  final, canonical version of the code - contact me to get that.       ///
///                                                                       ///
///  There is a prize of a pint for reporting any serious bugs or         ///
///  finding something new that results in >20% speed up.                 ///
///                                                                       ///
///                                                                       ///
///  The code below is structured as follows:                             ///
///                                                                       ///
///  0. SETTING UP STRUCTURES AND CLASSES                                 ///
///     All parameter values are stored in a structure called params.     ///
///     A class is created which stores all the information of a          ///
///     single individual.                                                ///
///	    A structure called population stores all individuals.             ///
///     The time-dependent output of the model is stored in a             ///
///     structure called simulation.                                      ///
///                                                                       ///
///  1. MAIN - SIMULATION                                                 ///
///     Here we read in parameters from files (model parameters,          ///
///     mosquito parameters, intervention parameters).                    ///
///     A population of individuals is created at equilibrium and         ///
///     then simulated.                                                   ///
///                                                                       ///
///  2. FUNCTIONS                                                         ///
///     Useful functions needed for model simulations.                    ///
///     Mosquitoes are simulated using a deterministic ODE solver.        ///
///     The functions for simulating a population of individuals are      ///
///     provided here - note that the model itself appears in Section 4.  ///
///                                                                       ///
///  3. EQUILIBRIUM SETUP                                                 ///
///     This set of functions calculates the equilibrium set up of the    ///
///     population. It is only called once while the population is        ///
///     being initialised.                                                ///
///                                                                       ///
///  4. INDIVIDUAL-BASED MODEL                                            ///
///     Details of the stochastic individual-based model for each         ///
///     person. Transitions occur with a fixed time step according to     ///
///     compting hazards                                                  ///
///                                                                       ///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <time.h>
#include "randlib.h"
#include <omp.h>
#include <vector>
#include <algorithm>
#include <regex>

using namespace std;

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

	double EIR_equil;			   // EIR at equilibrium


	/////////////////////////////////
	// Human demography

	double mu_H;                   // human death rate	 
	double age_mean;               // mean age of human population
	double age_max;                // maximum age in human population
	double het_max;                // maximum heterogeneity

	double G6PD_prev;              // Prevalence of G6PD deficiency
	double CYP2D6_prev;			   // Prevalence of CYP2D6 phenotype


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

	double Karry[N_spec];			// Larval carry capacity

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
	int	PQ_treat_G6PD_risk;         // Risk in G6PD - deficient individuals
	int	PQ_treat_CYP2D6_risk;       // Risk of not working in low CYP2D6 metabolizers
	int	PQ_treat_preg_risk;         // Risk in pregnant women
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
	double CHI_endo[N_spec];	// Endophily - proportion of mosquitoes resting indoors after feeding (no intervention)
	double PSI_indoors[N_spec]; // Proportion of bites taken on humans indoors
	double PSI_bed[N_spec];     // Proportion of bites taken on humans in bed

	double delta_1;	            // Time spent foraging for a blood meal
	double delta_2;	            // Time_spent_digesting_blood_meal
	double delta;	            // Duration of gonotrophic cycle

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


	///////////////////////////////////
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
	bool ACT_treat;    // 
	bool PQ_treat;     // 

	bool PQ_effective;
	bool PQ_overtreat;
	bool PQ_overtreat_9m;


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
	//
	// Note that while there are 9 interventions in total, 2 of these are
	// are related to first-line treatment. Therefore there are N_int = 6
	// 'pulsed' interventions, i.e. those distributed by campaigns where 
	// access will be an issue. 
	//
	// Of course there will be differential access to first-line treatment,
	// but that's a story for another day.

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


	///////////////////////////////
	// SSAT

	double T_last_BS;

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
	// First-line treatment - blood-stage drugs

	vector<double> BS_treat_year_on;
	vector<double> BS_treat_year_off;
	vector<double> BS_treat_BScover;
	vector<double> BS_treat_BSeff;
	vector<double> BS_treat_BSproph;


	////////////////////////////////
	// First-line treatment - blood-stage drugs
	// plus primaquine

	vector<double> PQ_treat_year_on;
	vector<double> PQ_treat_year_off;
	vector<double> PQ_treat_BScover;
	vector<double> PQ_treat_BSeff;
	vector<double> PQ_treat_BSproph;
	vector<double> PQ_treat_PQavail;
	vector<double> PQ_treat_PQeff;
	vector<double> PQ_treat_PQproph;
	vector<int>    PQ_treat_G6PD_risk;
	vector<int>    PQ_treat_CYP2D6_risk;
	vector<int>    PQ_treat_preg_risk;
	vector<double> PQ_treat_low_age;


	////////////////////////////////
	// MDA - blood-stage drugs

	vector<double> MDA_BS_year;
	vector<double> MDA_BS_BScover;
	vector<double> MDA_BS_BSeff;
	vector<double> MDA_BS_BSproph;


	////////////////////////////////
	// MDA - blood-stage drugs + primaquine

	vector<double> MDA_PQ_year;
	vector<double> MDA_PQ_BScover;
	vector<double> MDA_PQ_BSeff;
	vector<double> MDA_PQ_BSproph;
	vector<double> MDA_PQ_PQavail;
	vector<double> MDA_PQ_PQeff;
	vector<double> MDA_PQ_PQproph;
	vector<int>    MDA_PQ_G6PD_risk;
	vector<int>    MDA_PQ_CYP2D6_risk;
	vector<int>    MDA_PQ_preg_risk;
	vector<double> MDA_PQ_low_age;


	////////////////////////////////
	// MSAT - RDT or PCR with blood-stage drugs + primaquine

	vector<double> MSAT_PQ_year;
	vector<int>    MSAT_PQ_RDT_PCR;
	vector<double> MSAT_PQ_sens;
	vector<double> MSAT_PQ_spec;
	vector<double> MSAT_PQ_BScover;
	vector<double> MSAT_PQ_BSeff;
	vector<double> MSAT_PQ_BSproph;
	vector<double> MSAT_PQ_PQavail;
	vector<double> MSAT_PQ_PQeff;
	vector<double> MSAT_PQ_PQproph;
	vector<int>    MSAT_PQ_G6PD_risk;
	vector<int>    MSAT_PQ_CYP2D6_risk;
	vector<int>    MSAT_PQ_preg_risk;
	vector<double> MSAT_PQ_low_age;


	////////////////////////////////
	// SSAT - blood-stage drugs + primaquine

	vector<double> SSAT_PQ_year;
	vector<double> SSAT_PQ_sens;
	vector<double> SSAT_PQ_spec;
	vector<double> SSAT_PQ_BScover;
	vector<double> SSAT_PQ_BSeff;
	vector<double> SSAT_PQ_BSproph;
	vector<double> SSAT_PQ_PQavail;
	vector<double> SSAT_PQ_PQeff;
	vector<double> SSAT_PQ_PQproph;
	vector<int>    SSAT_PQ_G6PD_risk;
	vector<int>    SSAT_PQ_CYP2D6_risk;
	vector<int>    SSAT_PQ_preg_risk;
	vector<double> SSAT_PQ_low_age;
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

	vector<vector<int>> prev_all;   // Contains {N_pop, PvPR_PCR, PvPR_LM, Pv_clin, PvHR, PvHR_batches, new_PCR, new_LM, new_D, new_T} 
	vector<vector<int>> prev_U5;    // Contains {N_pop, PvPR_PCR, PvPR_LM, Pv_clin, PvHR, PvHR_batches, new_PCR, new_LM, new_D, new_T} 
	vector<vector<int>> prev_U10;   // Contains {N_pop, PvPR_PCR, PvPR_LM, Pv_clin, PvHR, PvHR_batches, new_PCR, new_LM, new_D, new_T} 


	////////////////////////////////////////
	// 0.5.3. Tracking coverage over time

	vector<int> LLIN_cov_t;
	vector<int> IRS_cov_t;
	vector<int> ACT_treat_t;
	vector<int> PQ_treat_t;
	vector<int> pregnant_t;

	vector<int> PQ_overtreat_t;
	vector<int> PQ_overtreat_9m_t;


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


////////////////////////////////////////////
////////////////////////////////////////////
//        //                              //
//   ##   //  #     #  ####  #### #   ##  //
//  ###   //  ##   ## ##  ##  ##  ##  ##  //
//   ##   //  ####### ######  ##  ### ##  //
//   ##   //  ## # ## ##  ##  ##  ## ###  //
//  ####  //  ##   ## ##  ## #### ##  ##  //
//        //                              //
////////////////////////////////////////////
////////////////////////////////////////////

int main(int argc, char** argv)
{
	setall(time(NULL), 7);

	clock_t clock_time;
	clock_time = clock();


	////////////////////////////////////////////
	//                                        //  
	//  1.1 Read in file names                //
	//                                        //
	////////////////////////////////////////////

	// do we have the correct command line?
	if (argc != 7)
	{
		std::cout << "Incorrect command line.\n";
		return 0;
	}

	char* parameter_File = argv[1];

	char* mosquito_File[N_spec_max];
	for (int g = 0; g < N_spec_max; g++)
	{
		mosquito_File[g] = argv[2 + g];
	}

	char* coverage_File = argv[5];
	char* output_File = argv[6];


	////////////////////////////////////////////
	//                                        //
	// 1.2. Initialise objects                //
	//                                        //
	////////////////////////////////////////////

	population PNG_pop;
	params Pv_mod_par;

	int time_start, time_end, burnin_time, N_time;


	////////////////////////////////////////////
	//                                        //
	// 1.3. Read in model parameters          //
	//                                        //
	////////////////////////////////////////////

	cout << "Reading in parameter file............." << endl;
	cout << endl;

	string discard;

	std::ifstream parameter_Stream(parameter_File);

	if (parameter_Stream.fail())
	{
		std::cout << "Failure reading in data." << endl;
	}


	//////////////////////////////////////////////////
	// Population size and simulation time

	parameter_Stream >> discard >> PNG_pop.N_pop >> discard;               // Number of participants

	parameter_Stream >> discard >> Pv_mod_par.EIR_equil >> discard;        // EIR at equilibrium

	for (int g = 0; g < N_spec_max; g++)
	{
		parameter_Stream >> discard >> Pv_mod_par.Prop_mosq[g] >> discard; // Proportions of each mosquito species
	}

	parameter_Stream >> discard >> time_start >> discard;                  // Start time for simulation
	parameter_Stream >> discard >> time_end >> discard;                    // End time for simulation
	parameter_Stream >> discard >> burnin_time >> discard;                 // End time for simulation

	N_time = (1 / t_step)*(burnin_time + time_end - time_start) * 365;     // Number of time steps for simulation


	/////////////////////////////////
	// Human demography
	// Note that the mean age isn't technically the mean 
	// due to the effects of truncation at maximum age

	parameter_Stream >> discard >> Pv_mod_par.age_mean >> discard;    // mean age of human population
	parameter_Stream >> discard >> Pv_mod_par.age_max >> discard;     // maximum age in human population

	Pv_mod_par.mu_H = 1.0 / Pv_mod_par.age_mean;                      // human death rate

	Pv_mod_par.het_max = 100.0;                                       // Maximum relative heterogeneity value


	////////////////////////////////////////////////////////
	// Heterogeneity in exposure and age-dependent biting

	parameter_Stream >> discard >> Pv_mod_par.age_0 >> discard;       // age-dependent biting parameter
	parameter_Stream >> discard >> Pv_mod_par.rho_age >> discard;     // age-dependent biting parameter
	parameter_Stream >> discard >> Pv_mod_par.sig_het >> discard;     // heterogeneity in exposure - standard deviation on a log scale


	/////////////////////////////////
	// Transmission probabilities

	parameter_Stream >> discard >> Pv_mod_par.bb >> discard;          // mosquito -> human transmission probability
	parameter_Stream >> discard >> Pv_mod_par.c_PCR >> discard;       // human -> mosquito transmission probability (PCR detectable)
	parameter_Stream >> discard >> Pv_mod_par.c_LM >> discard;        // human -> mosquito transmission probability (LM detectable)
	parameter_Stream >> discard >> Pv_mod_par.c_D >> discard;         // human -> mosquito transmission probability (clinical disease)
	parameter_Stream >> discard >> Pv_mod_par.c_T >> discard;         // human -> mosquito transmission probability (during treatment)


	////////////////////////////////
	// Human recovery paramters

	parameter_Stream >> discard >> Pv_mod_par.d_latent >> discard;             // latent period in liver	
	parameter_Stream >> discard >> Pv_mod_par.r_LM >> discard;                 // rate of recovery from LM detectable infection
	parameter_Stream >> discard >> Pv_mod_par.r_D >> discard;                  // rate of recovery from symptomatic disease
	parameter_Stream >> discard >> Pv_mod_par.r_T >> discard;                  // rate of progression through treatment 
	parameter_Stream >> discard >> Pv_mod_par.d_PCR_min >> discard;            // minimum duration of PCR-detectable infection - full immunity
	parameter_Stream >> discard >> Pv_mod_par.d_PCR_max >> discard;            // maximum duration of PCR-detectable infection - no immunity

	parameter_Stream >> discard >> Pv_mod_par.BS_treat_BScover_base >> discard;      // proportion of episodes of symptomatic disease treated (baseline)
	parameter_Stream >> discard >> Pv_mod_par.BS_treat_BSeff_base >> discard;      // efficacy of front-line treatment (baseline)
	parameter_Stream >> discard >> Pv_mod_par.BS_treat_BSproph_base >> discard;  // duration of prophylaxis of front-line treatment (baseline)

	parameter_Stream >> discard >> Pv_mod_par.A_PCR_50pc >> discard;           // PCR_detectable infection scale parameter
	parameter_Stream >> discard >> Pv_mod_par.K_PCR >> discard;                // PCR_detectable infection shape parameter


	Pv_mod_par.H_track = int(Pv_mod_par.d_latent / t_step);                       // Number of time steps for duration of latency

	Pv_mod_par.treat_BScover = Pv_mod_par.BS_treat_BScover_base;                  // Treatment coverage
	Pv_mod_par.treat_BSeff   = Pv_mod_par.BS_treat_BSeff_base;                    // Efficacy of treatment
	Pv_mod_par.r_P           = 1.0 / Pv_mod_par.BS_treat_BSproph_base;            // rate of recovery from prophylaxis
	Pv_mod_par.treat_PQavail = Pv_mod_par.MDA_PQ_PQavail;                         // PQ availability

	Pv_mod_par.PQ_treat_BScover     = 0.0;
	Pv_mod_par.PQ_treat_BSeff       = 0.0;
	Pv_mod_par.PQ_treat_BSproph     = 10.0;
	Pv_mod_par.PQ_treat_PQavail     = 0.0;
	Pv_mod_par.PQ_treat_PQeff       = 0.0;
	Pv_mod_par.PQ_treat_PQproph     = 0.0;
	Pv_mod_par.PQ_treat_G6PD_risk   = 1;
	Pv_mod_par.PQ_treat_CYP2D6_risk = 1;
	Pv_mod_par.PQ_treat_preg_risk   = 1;
	Pv_mod_par.PQ_treat_low_age     = 180.0;


	//////////////////
	// temporary setting treatment cov to zero
	//
	//Pv_mod_par.treat_cov_eff = 0.0;


	/////////////////////////////////
	// Blood-stage immunity paramters

	parameter_Stream >> discard >> Pv_mod_par.u_par >> discard;         // scale paramter for acquisition of blood-stage immunity
	parameter_Stream >> discard >> Pv_mod_par.r_par >> discard;         // rate of decay of blood-stage immunity

	parameter_Stream >> discard >> Pv_mod_par.phi_LM_max >> discard;    // probability of blood-stage infection with no immunity 
	parameter_Stream >> discard >> Pv_mod_par.phi_LM_min >> discard;    // probability of blood-stage infection with maximum immunity
	parameter_Stream >> discard >> Pv_mod_par.A_LM_50pc >> discard;     // blood-stage immunity scale parameter
	parameter_Stream >> discard >> Pv_mod_par.K_LM >> discard;          // blood-stage immunity shape parameter


	/////////////////////////////////
	// Clinical immunity paramters

	parameter_Stream >> discard >> Pv_mod_par.u_clin >> discard;       // scale paramter for acquisition of blood-stage immunity
	parameter_Stream >> discard >> Pv_mod_par.r_clin >> discard;       // rate of decay of clinical immunity

	parameter_Stream >> discard >> Pv_mod_par.phi_D_max >> discard;    // probability of clinical episode with no immunity
	parameter_Stream >> discard >> Pv_mod_par.phi_D_min >> discard;    // probability of clinical episode with maximum immunity
	parameter_Stream >> discard >> Pv_mod_par.A_D_50pc >> discard;     // clinical immunity scale parameter
	parameter_Stream >> discard >> Pv_mod_par.K_D >> discard;          // clinical immunity shape parameter


	/////////////////////////////////////////////
	// maternal immunity

	parameter_Stream >> discard >> Pv_mod_par.P_mat >> discard;        // New-born immunity relative to mother's
	parameter_Stream >> discard >> Pv_mod_par.d_mat >> discard;        // Inverse of decay rate of maternal immunity


	/////////////////////////////////
	// Relapse paramters

	parameter_Stream >> discard >> Pv_mod_par.ff >> discard;           // relapse rate
	parameter_Stream >> discard >> Pv_mod_par.gamma_L >> discard;      // liver clearance rate


	////////////////////////////////
	// Human genotype prevalences

	parameter_Stream >> discard >> Pv_mod_par.G6PD_prev >> discard;    // prevalence of G6PD deficiency
	parameter_Stream >> discard >> Pv_mod_par.CYP2D6_prev >> discard;  // prevalence of CYP2D6 phenotype 


	////////////////////////////////////////////////////////
	// Intervention distribution parameters

	parameter_Stream >> discard >> Pv_mod_par.rho_round_LLIN >> discard;
	parameter_Stream >> discard >> Pv_mod_par.rho_round_IRS >> discard;
	parameter_Stream >> discard >> Pv_mod_par.rho_round_MDA >> discard;

	parameter_Stream >> discard >> Pv_mod_par.rho_LLIN_IRS >> discard;
	parameter_Stream >> discard >> Pv_mod_par.rho_MDA_IRS >> discard;
	parameter_Stream >> discard >> Pv_mod_par.rho_MDA_LLIN >> discard;


	Pv_mod_par.sig_round_LLIN = sqrt((1.0 - Pv_mod_par.rho_round_LLIN) / Pv_mod_par.rho_round_LLIN);
	Pv_mod_par.sig_round_IRS = sqrt((1.0 - Pv_mod_par.rho_round_IRS) / Pv_mod_par.rho_round_IRS);
	Pv_mod_par.sig_round_MDA = sqrt((1.0 - Pv_mod_par.rho_round_MDA) / Pv_mod_par.rho_round_MDA);


	//////////////////////////////////////////////
	// TO DO - no need to have different possible correlations
	//         for the various MDA interventions - can assume access
	//         to individuals is all the same

	Pv_mod_par.V_int[0][0] = 1.0;
	Pv_mod_par.V_int[0][1] = Pv_mod_par.rho_LLIN_IRS;
	Pv_mod_par.V_int[0][2] = Pv_mod_par.rho_MDA_LLIN;
	Pv_mod_par.V_int[0][3] = Pv_mod_par.rho_MDA_LLIN;
	Pv_mod_par.V_int[0][4] = Pv_mod_par.rho_MDA_LLIN;
	Pv_mod_par.V_int[0][5] = Pv_mod_par.rho_MDA_LLIN;

	Pv_mod_par.V_int[1][0] = Pv_mod_par.rho_LLIN_IRS;
	Pv_mod_par.V_int[1][1] = 1.0;
	Pv_mod_par.V_int[1][2] = Pv_mod_par.rho_MDA_IRS;
	Pv_mod_par.V_int[1][3] = Pv_mod_par.rho_MDA_IRS;
	Pv_mod_par.V_int[1][4] = Pv_mod_par.rho_MDA_IRS;
	Pv_mod_par.V_int[1][5] = Pv_mod_par.rho_MDA_IRS;

	Pv_mod_par.V_int[2][0] = Pv_mod_par.rho_MDA_LLIN;
	Pv_mod_par.V_int[2][1] = Pv_mod_par.rho_MDA_IRS;
	Pv_mod_par.V_int[2][2] = 1.0;
	Pv_mod_par.V_int[2][3] = Pv_mod_par.rho_round_MDA;
	Pv_mod_par.V_int[2][4] = Pv_mod_par.rho_round_MDA;
	Pv_mod_par.V_int[2][5] = Pv_mod_par.rho_round_MDA;

	Pv_mod_par.V_int[3][0] = Pv_mod_par.rho_MDA_LLIN;
	Pv_mod_par.V_int[3][1] = Pv_mod_par.rho_MDA_IRS;
	Pv_mod_par.V_int[3][2] = Pv_mod_par.rho_round_MDA;
	Pv_mod_par.V_int[3][3] = 1.0;
	Pv_mod_par.V_int[3][4] = Pv_mod_par.rho_round_MDA;
	Pv_mod_par.V_int[3][5] = Pv_mod_par.rho_round_MDA;

	Pv_mod_par.V_int[4][0] = Pv_mod_par.rho_MDA_LLIN;
	Pv_mod_par.V_int[4][1] = Pv_mod_par.rho_MDA_IRS;
	Pv_mod_par.V_int[4][2] = Pv_mod_par.rho_round_MDA;
	Pv_mod_par.V_int[4][3] = Pv_mod_par.rho_round_MDA;
	Pv_mod_par.V_int[4][4] = 1.0;
	Pv_mod_par.V_int[4][5] = Pv_mod_par.rho_round_MDA;

	Pv_mod_par.V_int[5][0] = Pv_mod_par.rho_MDA_LLIN;
	Pv_mod_par.V_int[5][1] = Pv_mod_par.rho_MDA_IRS;
	Pv_mod_par.V_int[5][2] = Pv_mod_par.rho_round_MDA;
	Pv_mod_par.V_int[5][3] = Pv_mod_par.rho_round_MDA;
	Pv_mod_par.V_int[5][4] = Pv_mod_par.rho_round_MDA;
	Pv_mod_par.V_int[5][5] = 1.0;


	/////////////////////////////////////////////////////////
	// We need to make a dummy covariance matrix for intervention
	// distribution because of the way genmn works

	for (int p = 0; p<N_int; p++)
	{
		for (int q = 0; q<N_int; q++)
		{
			Pv_mod_par.V_int_dummy[p][q] = Pv_mod_par.V_int[p][q];
		}
	}

	parameter_Stream.close();

	cout << "Parameter values read in from file!" << endl;
	cout << endl;


	////////////////////////////////////////////
	//                                        //
	// 1.4. Read in mosquito parameters       //
	//                                        //
	////////////////////////////////////////////

	cout << "Reading in mosquito files............." << endl;
	cout << endl;

	for (int g = 0; g < N_spec; g++)
	{
		std::ifstream mosquito_Stream(mosquito_File[g]);

		if (mosquito_Stream.fail())
		{
			std::cout << "Failure reading in mosquito parameters." << endl;
		}


		/////////////////////////////////
		// Death rate and duration of sporogony

		mosquito_Stream >> discard >> Pv_mod_par.mu_M[g] >> discard;      // mosquito death rate
		mosquito_Stream >> discard >> Pv_mod_par.tau_M[g] >> discard;     // duration of sporogony

		Pv_mod_par.M_track = int(Pv_mod_par.tau_M[g] * mosq_steps / t_step);


		/////////////////////////////////
		// Larval paramters	
		// From White et al (2011) P&V

		mosquito_Stream >> discard >> Pv_mod_par.d_E_larvae >> discard;      // Development time of early larval instars
		mosquito_Stream >> discard >> Pv_mod_par.d_L_larvae >> discard;      // Development time of late larval instars
		mosquito_Stream >> discard >> Pv_mod_par.d_pupae >> discard;         // Development time of pupae
		mosquito_Stream >> discard >> Pv_mod_par.mu_E0 >> discard;           // Mortality rate of early larval instars (low density)
		mosquito_Stream >> discard >> Pv_mod_par.mu_L0 >> discard;           // Mortality rate of late larval instars (low density)
		mosquito_Stream >> discard >> Pv_mod_par.mu_P >> discard;            // Mortality rate of pupae
		mosquito_Stream >> discard >> Pv_mod_par.beta_larvae >> discard;     // Number of eggs laid per day per mosquito
		mosquito_Stream >> discard >> Pv_mod_par.gamma_larvae >> discard;    // Effect of density dependence on late instars relative to early instars

		Pv_mod_par.omega_larvae[g] = Pv_mod_par.gamma_larvae*Pv_mod_par.mu_L0 / Pv_mod_par.mu_E0 - Pv_mod_par.d_E_larvae / Pv_mod_par.d_L_larvae + (Pv_mod_par.gamma_larvae - 1.0)*Pv_mod_par.mu_L0*Pv_mod_par.d_E_larvae;
		Pv_mod_par.omega_larvae[g] = - 0.5*Pv_mod_par.omega_larvae[g] + sqrt(0.25*Pv_mod_par.omega_larvae[g] * Pv_mod_par.omega_larvae[g] + 0.5*Pv_mod_par.gamma_larvae*Pv_mod_par.beta_larvae*Pv_mod_par.mu_L0*Pv_mod_par.d_E_larvae /
			                           (Pv_mod_par.mu_E0*Pv_mod_par.mu_M[g] * Pv_mod_par.d_L_larvae*(1.0 + Pv_mod_par.d_pupae*Pv_mod_par.mu_P)));


		/////////////////////////////////
		// Seasonality paramters
		// Denominator for seasonality - see Griffin (2015) PLoS Comp Biol

		mosquito_Stream >> discard >> Pv_mod_par.dry_seas[g] >> discard;      // Proportion of dry season transmission compared to mean 
		mosquito_Stream >> discard >> Pv_mod_par.kappa_seas[g] >> discard;    // Shape parameter for seasonality
		mosquito_Stream >> discard >> Pv_mod_par.t_peak_seas[g] >> discard;   // Timing of peak for seasonal transmission

		Pv_mod_par.denom_seas[g] = exp(gammln(0.5) + gammln(Pv_mod_par.kappa_seas[g] + 0.5) - gammln(Pv_mod_par.kappa_seas[g] + 1.0)) / 3.14159265359;


		//////////////////////////////////////////////
		// Entomology paramters

		mosquito_Stream >> discard >> Pv_mod_par.Q_0[g] >> discard;           // Human Blood Index (proportion of blood meals taken on humans)
		mosquito_Stream >> discard >> Pv_mod_par.CHI_endo[g] >> discard;	  // Endophily - proportion of mosquitoes resting indoors after feeding (no intervention)
		mosquito_Stream >> discard >> Pv_mod_par.PSI_indoors[g] >> discard;   // Proportion of bites taken on humans indoors
		mosquito_Stream >> discard >> Pv_mod_par.PSI_bed[g] >> discard;       // Proportion of bites taken on humans in bed

		mosquito_Stream >> discard >> Pv_mod_par.delta_1 >> discard;	      // Time spent foraging for a blood meal
		mosquito_Stream >> discard >> Pv_mod_par.delta >> discard;            // Duration of gonotrophic cycle

		Pv_mod_par.delta_2 = Pv_mod_par.delta - Pv_mod_par.delta_1;

		Pv_mod_par.p_1[g] = exp(-Pv_mod_par.mu_M[g] * Pv_mod_par.delta_1);
		Pv_mod_par.p_2[g] = exp(-Pv_mod_par.mu_M[g] * Pv_mod_par.delta_2);

		Pv_mod_par.aa[g] = Pv_mod_par.Q_0[g] / (Pv_mod_par.delta_1 + Pv_mod_par.delta_2);

		Pv_mod_par.eps_max[g] = Pv_mod_par.beta_larvae*(exp(Pv_mod_par.delta*Pv_mod_par.mu_M[g]) - 1.0) / Pv_mod_par.mu_M[g];


		//////////////////////////////////////////////
		// LLIN paramters

		mosquito_Stream >> discard >> Pv_mod_par.LLIN_half_life >> discard;

		mosquito_Stream >> discard >> Pv_mod_par.PYR_half_life >> discard;

		mosquito_Stream >> discard >> Pv_mod_par.r_LLIN_0[g] >> discard;       // Probability mosquito repelled (with full insecticide activity)
		mosquito_Stream >> discard >> Pv_mod_par.r_LLIN_net[g] >> discard;     // Probability mosquito repelled due to barrier effect of net (no insecticide)
		mosquito_Stream >> discard >> Pv_mod_par.d_LLIN_0[g] >> discard;       // Probability mosquito dies during feeding attempt


		Pv_mod_par.P_LLIN_loss = 1.0 - exp(-t_step*log(2.0) / Pv_mod_par.LLIN_half_life);   // Probability of losing LLIN in a time step
		Pv_mod_par.PYR_decay = log(2.0) / Pv_mod_par.PYR_half_life;                            // Rate of pyrethroid decay

		Pv_mod_par.s_LLIN_0[g] = 1.0 - Pv_mod_par.r_LLIN_0[g] - Pv_mod_par.d_LLIN_0[g];     // Probability mosquito feeds successfully


		////////////////////////////////////////////////////////
		// IRS parameters

		mosquito_Stream >> discard >> Pv_mod_par.IRS_half_life >> discard;    // IRS insecticide half-life

		mosquito_Stream >> discard >> Pv_mod_par.r_IRS_0[g] >> discard;          // IRS repellency
		mosquito_Stream >> discard >> Pv_mod_par.d_IRS_0[g] >> discard;          // IRS death

		Pv_mod_par.IRS_decay = log(2.0) / Pv_mod_par.IRS_half_life;         // IRS decay rate

		Pv_mod_par.s_IRS_0[g] = 1.0 - Pv_mod_par.d_IRS_0[g] - Pv_mod_par.s_IRS_0[g];   // Feeding success of mosquito on IRS protected person
	}

	cout << "Mosquito parameter values read in from file!" << endl;
	cout << endl;


	//////////////////////////////////////////////
	// Normalise relative proprotions of 
	// different mosquito species

	double Prop_mosq_denom = 0.0;

	for (int g = 0; g < N_spec; g++)
	{
		Prop_mosq_denom = Prop_mosq_denom + Pv_mod_par.Prop_mosq[g];
	}

	for (int g = 0; g < N_spec; g++)
	{
		Pv_mod_par.Prop_mosq[g] = Pv_mod_par.Prop_mosq[g] / Prop_mosq_denom;
	}


	///////////////////////////////////////////////////////////
	//                                                       //
	// 1.5. Pre-multiplication of quantities for efficiency  //
	//                                                       //
	///////////////////////////////////////////////////////////

	Pv_mod_par.A_par_decay  = exp(-Pv_mod_par.r_par*t_step);
	Pv_mod_par.A_clin_decay = exp(-Pv_mod_par.r_clin*t_step);
	Pv_mod_par.mat_decay    = exp(-Pv_mod_par.d_mat*t_step);

	Pv_mod_par.age_0_inv = 1.0 / Pv_mod_par.age_0;                 // Inverse of age-dependent biting parameter

	Pv_mod_par.A_PCR_50pc_inv = log2 / Pv_mod_par.A_PCR_50pc;      // Immune scalar for clearance of infection
	Pv_mod_par.A_LM_50pc_inv  = 1.0 / Pv_mod_par.A_LM_50pc;        // Immune scalar for BS infection
	Pv_mod_par.A_D_50pc_inv   = 1.0 / Pv_mod_par.A_D_50pc;         // Immune scalar for clinical disease

	Pv_mod_par.P_dead = 1.0 - exp(-t_step*Pv_mod_par.mu_H);
	Pv_mod_par.P_preg = 0.0014189;


	Pv_mod_par.P_PYR_decay = exp(-Pv_mod_par.PYR_decay*t_step);
	Pv_mod_par.P_IRS_decay = exp(-Pv_mod_par.IRS_decay*t_step);


	///////////////////////////////////////////////////////////
	//                                                       //
	// 1.6. Fill out hypnozoite transition matrices          //
	//                                                       //
	///////////////////////////////////////////////////////////

	for (int k1 = 0; k1<(K_max + 1); k1++)
	{
		for (int k2 = 0; k2 < (K_max + 1); k2++)
		{
			Pv_mod_par.D_MAT[k1][k2] = 0.0;
			Pv_mod_par.OD_MAT[k1][k2] = 0.0;
			Pv_mod_par.K_MAT[k1][k2] = 0.0;
			Pv_mod_par.L_MAT[k1][k2] = 0.0;
			Pv_mod_par.H_MAT[k1][k2] = 0.0;
		}
	}

	for (int k = 0; k < (K_max + 1); k++)
	{
		Pv_mod_par.D_MAT[k][k] = 1.0;
	}

	for (int k = 0; k < K_max; k++)
	{
		Pv_mod_par.OD_MAT[k + 1][k] = 1.0;
	}
	Pv_mod_par.OD_MAT[K_max][K_max] = 1.0;

	for (int k = 0; k < (K_max + 1); k++)
	{
		Pv_mod_par.K_MAT[k][k] = (double)(k);
	}

	for (int k = 0; k < K_max; k++)
	{
		Pv_mod_par.L_MAT[k][k + 1] = +(double)(k + 1);
		Pv_mod_par.L_MAT[k + 1][k + 1] = -(double)(k + 1);
	}

	for (int k = 0; k < K_max; k++)
	{
		Pv_mod_par.H_MAT[k][k] = -1.0;
		Pv_mod_par.H_MAT[k + 1][k] = +1.0;
	}


	/////////////////////////////////////////////////////////////////////////
	//                                                                     //
	// 1.7. Read in intervention coverage                                  //
	//                                                                     //
	/////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////
	// 1.7.1. Read in matrix of years and coverages

	cout << "Read in intervention coverage file............" << endl;

	std::ifstream coverage_Stream(coverage_File);

	if (coverage_Stream.fail())
	{
		std::cout << "Failure reading in data." << endl;
	}


	/////////////////////////////////////////////////////////////////////////
	// 1.7.2. Note that the matrix we read in may have variable size.
	//        We first read in the first line of the file, and then matching
	//        the spaces between columns to get the number of interventions 
	//        rounds.
	//
	//		  There's very likely a much more effective way to do this.

	string line_read;
	int N_cov_rounds = 0;

	getline(coverage_Stream, line_read);


	std::string str(line_read);
	std::string str2(" ");

	std::size_t space_find = str.find(str2);


	// Go through the string finding spaces (" ")

	while (space_find < 10000)
	{
		space_find = line_read.find(str2);

		line_read.erase(0, space_find + 1);

		N_cov_rounds = N_cov_rounds + 1;
	}

	N_cov_rounds = N_cov_rounds - 2;

	cout << "Number of intervention rounds: " << N_cov_rounds << endl;

	coverage_Stream.close();


	/////////////////////////////////////////////////////////////////////////
	// 1.7.3. Note that the matrix we read in may have variable size.

	std::ifstream coverage_Stream2(coverage_File);

	vector<vector<double>> coverage;
	coverage.resize(0);

	vector<double> cov_read(N_cov_rounds + 1);

	for (int i = 0; i < 55; i++)
	{
		coverage_Stream2 >> discard;

		for (int j = 1; j<(N_cov_rounds + 2); j++)
		{
			coverage_Stream2 >> cov_read[j-1];
		}

		coverage.push_back(cov_read);
	}


	/////////////////////////////////////////////////////////////////////////
	// 1.7.2. Fill out intervention structure

	intervention PNG_intven;

	for (int j = 0; j<N_cov_rounds; j++)
	{
		//////////////////////////////////////////////////////////////
		// LLINs

		if ((coverage[0][j] > -0.5) && (coverage[1][j] > -0.5))
		{
			PNG_intven.LLIN_year.push_back(  coverage[0][j]*365.0 );
			PNG_intven.LLIN_cover.push_back( coverage[1][j] );
		}


		//////////////////////////////////////////////////////////////
		// IRS

		if ((coverage[0][j] > -0.5) && (coverage[2][j] > -0.5))
		{
			PNG_intven.IRS_year.push_back(  coverage[0][j]*365.0 );
			PNG_intven.IRS_cover.push_back( coverage[2][j] );
		}


		//////////////////////////////////////////////////////////////
		// Front-line treatment - blood-stage drugs

		if ((coverage[0][j] > -0.5) && (coverage[4][j] > -0.5))
		{
			PNG_intven.BS_treat_year_on.push_back(  coverage[0][j]*365.0 );
			PNG_intven.BS_treat_year_off.push_back( coverage[3][j]*365.0 );
			PNG_intven.BS_treat_BScover.push_back(  coverage[4][j] );
			PNG_intven.BS_treat_BSeff.push_back(    coverage[5][j] );
			PNG_intven.BS_treat_BSproph.push_back(  coverage[6][j] );
		}


		//////////////////////////////////////////////////////////////
		// Front-line treatment - primaquine

		if ((coverage[0][j] > -0.5) && (coverage[8][j] > -0.5))
		{
			PNG_intven.PQ_treat_year_on.push_back(     coverage[0][j]*365.0 );
			PNG_intven.PQ_treat_year_off.push_back(    coverage[7][j]*365.0 );
			PNG_intven.PQ_treat_BScover.push_back(     coverage[8][j] );
			PNG_intven.PQ_treat_BSeff.push_back(       coverage[9][j] );
			PNG_intven.PQ_treat_BSproph.push_back(     coverage[10][j] );
			PNG_intven.PQ_treat_PQavail.push_back(     coverage[11][j] );
			PNG_intven.PQ_treat_PQeff.push_back(       coverage[12][j] );
			PNG_intven.PQ_treat_PQproph.push_back(     coverage[13][j] );
			PNG_intven.PQ_treat_G6PD_risk.push_back(   (int)(coverage[14][j]) );
			PNG_intven.PQ_treat_CYP2D6_risk.push_back( (int)(coverage[15][j]) );
			PNG_intven.PQ_treat_preg_risk.push_back(   (int)(coverage[16][j]) );
			PNG_intven.PQ_treat_low_age.push_back(     coverage[17][j] );
		}


		//////////////////////////////////////////////////////////////
		// MDA - blood-stage drugs

		if ((coverage[0][j] > -0.5) && (coverage[18][j] > -0.5))
		{
			PNG_intven.MDA_BS_year.push_back(    coverage[0][j]*365.0 );
			PNG_intven.MDA_BS_BScover.push_back( coverage[18][j] );
			PNG_intven.MDA_BS_BSeff.push_back(   coverage[19][j] );
			PNG_intven.MDA_BS_BSproph.push_back( coverage[20][j] );
		}


		//////////////////////////////////////////////////////////////
		// MDA - blood-stage drugs plus primaquine

		if ((coverage[0][j] > -0.5) && (coverage[21][j] > -0.5))
		{
			PNG_intven.MDA_PQ_year.push_back(        coverage[0][j]*365.0 );
			PNG_intven.MDA_PQ_BScover.push_back(     coverage[21][j] );
			PNG_intven.MDA_PQ_BSeff.push_back(       coverage[22][j] );
			PNG_intven.MDA_PQ_BSproph.push_back(     coverage[23][j] );
			PNG_intven.MDA_PQ_PQavail.push_back(     coverage[24][j] );
			PNG_intven.MDA_PQ_PQeff.push_back(       coverage[25][j] );
			PNG_intven.MDA_PQ_PQproph.push_back(     coverage[26][j] );
			PNG_intven.MDA_PQ_G6PD_risk.push_back(   (int)(coverage[27][j]) );
			PNG_intven.MDA_PQ_CYP2D6_risk.push_back( (int)(coverage[28][j]) );
			PNG_intven.MDA_PQ_preg_risk.push_back(   (int)(coverage[29][j]) );
			PNG_intven.MDA_PQ_low_age.push_back(     coverage[30][j] );
		}


		//////////////////////////////////////////////////////////////
		// MSAT - blood-stage drugs plus primaquine

		if ((coverage[0][j] > -0.5) && (coverage[31][j] > -0.5))
		{
			PNG_intven.MSAT_PQ_year.push_back(        coverage[0][j]*365.0 );
			PNG_intven.MSAT_PQ_BScover.push_back(     coverage[31][j]);
			PNG_intven.MSAT_PQ_RDT_PCR.push_back(     coverage[32][j] );
			PNG_intven.MSAT_PQ_sens.push_back(        coverage[33][j] );
			PNG_intven.MSAT_PQ_BSeff.push_back(       coverage[34][j] );
			PNG_intven.MSAT_PQ_BSproph.push_back(     coverage[35][j] );
			PNG_intven.MSAT_PQ_PQavail.push_back(     coverage[36][j] );
			PNG_intven.MSAT_PQ_PQeff.push_back(       coverage[37][j] );
			PNG_intven.MSAT_PQ_PQproph.push_back(     coverage[38][j] );
			PNG_intven.MSAT_PQ_G6PD_risk.push_back(   (int)(coverage[39][j]) );
			PNG_intven.MSAT_PQ_CYP2D6_risk.push_back( (int)(coverage[40][j]) );
			PNG_intven.MSAT_PQ_preg_risk.push_back(   (int)(coverage[41][j]) );
			PNG_intven.MSAT_PQ_low_age.push_back(     coverage[42][j] );
		}

		//////////////////////////////////////////////////////////////
		// SSAT - blood-stage drugs plus primaquine

		if ((coverage[0][j] > -0.5) && (coverage[43][j] > -0.5))
		{
			PNG_intven.SSAT_PQ_year.push_back(        coverage[0][j]*365.0 );
			PNG_intven.SSAT_PQ_BScover.push_back(     coverage[43][j] );
			PNG_intven.SSAT_PQ_sens.push_back(        coverage[44][j] );
			PNG_intven.SSAT_PQ_spec.push_back(        coverage[45][j] );
			PNG_intven.SSAT_PQ_BSeff.push_back(       coverage[46][j] );
			PNG_intven.SSAT_PQ_BSproph.push_back(     coverage[47][j] );
			PNG_intven.SSAT_PQ_PQavail.push_back(     coverage[48][j] );
			PNG_intven.SSAT_PQ_PQeff.push_back(       coverage[49][j] );
			PNG_intven.SSAT_PQ_PQproph.push_back(     coverage[50][j] );
			PNG_intven.SSAT_PQ_G6PD_risk.push_back(   (int)(coverage[51][j]) );
			PNG_intven.SSAT_PQ_CYP2D6_risk.push_back( (int)(coverage[52][j]) );
			PNG_intven.SSAT_PQ_preg_risk.push_back(   (int)(coverage[53][j]) );
			PNG_intven.SSAT_PQ_low_age.push_back(     coverage[54][j] );
		}

	}


	///////////////////////////////////////////////////////////////////////////
	//                                                                       //
	// 1.8. Initialise Population of individuals                             //
	//      Note that they begin with exponential age distribution           //
	//      and susceptible without immunity                                 //
	//                                                                       // 
	///////////////////////////////////////////////////////////////////////////

	cout << "Initialise population of individuals for simulation at equilbirium EIR of " << 365.0*Pv_mod_par.EIR_equil << endl;
	cout << endl;

	equi_pop_setup(&PNG_pop, &Pv_mod_par);

	cout << "Population of size " << PNG_pop.N_pop << " initialised!" << endl;
	cout << endl;


	/////////////////////////////////////////////////////////////////////////
	//                                                                     //
	// 1.9. Create simulation object                                       //
	//                                                                     //
	/////////////////////////////////////////////////////////////////////////

	simulation PNG_sim;


	/////////////////////////////////////////////////////////////////////////
	// 1.9.1. Vector of simulation times

	PNG_sim.N_time = N_time;

	for (int i = 0; i<N_time; i++)
	{
		PNG_sim.t_vec.push_back((double)(time_start * 365 - burnin_time * 365 + i*t_step));
	}


	/////////////////////////////////////////////////////////////////////////
	// 1.9.2. Create storage for output

	PNG_sim.yH_t.resize(N_time);
	for (int i = 0; i<N_time; i++)
	{
		PNG_sim.yH_t[i].resize(N_H_comp);
	}


	PNG_sim.yM_t.resize(N_time);
	for (int i = 0; i<N_time; i++)
	{
		PNG_sim.yM_t[i].resize(N_spec);
		for (int g = 0; g < N_spec; g++)
		{
			PNG_sim.yM_t[i][g].resize(N_M_comp);
		}
	}


	PNG_sim.prev_all.resize(N_time);
	for (int i = 0; i<N_time; i++)
	{
		PNG_sim.prev_all[i].resize(11);
	}

	PNG_sim.prev_U5.resize(N_time);
	for (int i = 0; i<N_time; i++)
	{
		PNG_sim.prev_U5[i].resize(11);
	}

	PNG_sim.prev_U10.resize(N_time);
	for (int i = 0; i<N_time; i++)
	{
		PNG_sim.prev_U10[i].resize(11);
	}

	PNG_sim.EIR_t.resize(N_time);

	PNG_sim.LLIN_cov_t.resize(N_time);
	PNG_sim.IRS_cov_t.resize(N_time);
	PNG_sim.ACT_treat_t.resize(N_time);
	PNG_sim.PQ_treat_t.resize(N_time);
	PNG_sim.pregnant_t.resize(N_time);

	PNG_sim.PQ_overtreat_t.resize(N_time);
	PNG_sim.PQ_overtreat_9m_t.resize(N_time);

	PNG_sim.A_par_mean_t.resize(N_time);
	PNG_sim.A_clin_mean_t.resize(N_time);


	//////////////////////////////////////////////////////
	//                                                  //
	// 1.10. Begin stochastic simulations               //
	//                                                  //
	////////////////////////////////////////////////////// 

	cout << "Starting model simulations......." << endl;

	model_simulator(&Pv_mod_par, &PNG_pop, &PNG_intven, &PNG_sim);

	cout << "Model simulations completed....." << endl;
	cout << endl;


	//////////////////////////////////////////////////////
	//                                                  //
	// 1.11. Output to file                             //
	//                                                  //
	////////////////////////////////////////////////////// 

	cout << "Start writing output to file......" << endl;
	cout << endl;

	ofstream output_Stream(output_File);

	for (int i = (int) (1/t_step)*(burnin_time)*365; i<N_time; i++)
	{
		output_Stream << PNG_sim.t_vec[i] << "\t";

		for (int k = 0; k<N_H_comp; k++)
		{
			output_Stream << PNG_sim.yH_t[i][k] << "\t";
		}

		for (int g = 0; g < N_spec; g++)
		{
			for (int k = 3; k < N_M_comp; k++)
			{
				output_Stream << PNG_sim.yM_t[i][g][k] << "\t";
			}
		}

		for (int k = 0; k<10; k++)
		{
			output_Stream << PNG_sim.prev_all[i][k] << "\t";
		}

		output_Stream << PNG_sim.EIR_t[i] << "\t";
		output_Stream << PNG_sim.LLIN_cov_t[i] << "\t";
		output_Stream << PNG_sim.IRS_cov_t[i] << "\t";
		output_Stream << PNG_sim.ACT_treat_t[i] << "\t";
		output_Stream << PNG_sim.PQ_treat_t[i] << "\t";

		output_Stream << PNG_sim.PQ_overtreat_t[i] << "\t";
		output_Stream << PNG_sim.PQ_overtreat_9m_t[i] << "\t";

		output_Stream << endl;
	}

	output_Stream.close();


	cout << "Output successfully written to file......" << endl;
	cout << endl;


	cout << "Time taken: " << ((double)clock() - clock_time) / ((double)CLOCKS_PER_SEC) << " seconds" << endl;


	return 0;
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//          //                                                              //
//   ####   //  ##### ##  ## #   ##  ####  ###### ####  ####  #   ##  ###   //
//  ##  ##  //  ##    ##  ## ##  ## ##  ##   ##    ##  ##  ## ##  ## ##     // 
//     ##   //  ####  ##  ## ### ## ##       ##    ##  ##  ## ### ##  ###   //
//    ##    //  ##    ##  ## ## ### ##  ##   ##    ##  ##  ## ## ###    ##  //
//   #####  //  ##     ####  ##  ##  ####    ##   ####  ####  ##  ##  ###   //
//          //                                                              //
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.1. Derivatives of mosquito ODE model                                  //
//                                                                          //
//  0.01721421 = 2*pi/365                                                   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void mosq_derivs(const double t, double(&yM)[N_spec][N_M_comp], double(&dyMdt)[N_spec][N_M_comp], params* theta, population* POP)
{
	double Karry_seas_inv[N_spec];

	for (int g = 0; g < N_spec; g++)
	{
		Karry_seas_inv[g] = 1.0 / (theta->Karry[g] * (theta->dry_seas[g] + (1 - theta->dry_seas[g])*pow(0.5 + 0.5*cos(0.01721421*(t - theta->t_peak_seas[g])), theta->kappa_seas[g]) / theta->denom_seas[g]));

		//Karry_seas_inv[g] = 1.0/theta->Karry[g];

		dyMdt[g][0] = POP->beta_VC[g] * (yM[g][3] + yM[g][4] + yM[g][5]) - yM[g][0] / theta->d_E_larvae - yM[g][0] * theta->mu_E0*(1.0 + (yM[g][0] + yM[g][1])*Karry_seas_inv[g]);
		dyMdt[g][1] = yM[g][0] / theta->d_E_larvae - yM[g][1] / theta->d_L_larvae - yM[g][1] * theta->mu_L0*(1.0 + theta->gamma_larvae*(yM[g][0] + yM[g][1])*Karry_seas_inv[g]);
		dyMdt[g][2] = yM[g][1] / theta->d_L_larvae - yM[g][2] / theta->d_pupae - yM[g][2] * theta->mu_P;
		dyMdt[g][3] = 0.5*yM[g][2] / theta->d_pupae - theta->lam_M[g] * yM[g][3] - POP->mu_M_VC[g] * yM[g][3];
		dyMdt[g][4] = +theta->lam_M[g] * yM[g][3] - theta->lam_S_M_track[g][0] * POP->exp_muM_tauM_VC[g] - POP->mu_M_VC[g] * yM[g][4];
		dyMdt[g][5] = +theta->lam_S_M_track[g][0] * POP->exp_muM_tauM_VC[g] - POP->mu_M_VC[g] * yM[g][5];
	}
}

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.2. Runge-Kutta 4 step updater for mosquito model                      //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void mosq_rk4(const double t, const double t_step_mosq, double(&yM)[N_spec][N_M_comp], params* theta, population* POP)
{
	double k1_yM[N_spec][N_M_comp], k2_yM[N_spec][N_M_comp], k3_yM[N_spec][N_M_comp], k4_yM[N_spec][N_M_comp], yM_temp[N_spec][N_M_comp];


	//////////////////////////
	// step 1

	mosq_derivs(t, yM, k1_yM, theta, POP);


	//////////////////////////
	// step 2

	for (int g = 0; g < N_spec; g++)
	{
		for (int k = 0; k < N_M_comp; k++)
		{
			yM_temp[g][k] = yM[g][k] + 0.5*t_step_mosq*k1_yM[g][k];
		}
	}

	mosq_derivs(t + 0.5*t_step_mosq, yM_temp, k2_yM, theta, POP);


	//////////////////////////
	// step 3

	for (int g = 0; g < N_spec; g++)
	{
		for (int k = 0; k < N_M_comp; k++)
		{
			yM_temp[g][k] = yM[g][k] + 0.5*t_step_mosq*k2_yM[g][k];
		}
	}

	mosq_derivs(t + 0.5*t_step_mosq, yM_temp, k3_yM, theta, POP);


	//////////////////////////
	// step 4

	for (int g = 0; g < N_spec; g++)
	{
		for (int k = 0; k < N_M_comp; k++)
		{
			yM_temp[g][k] = yM[g][k] + t_step_mosq*k3_yM[g][k];
		}
	}

	mosq_derivs(t + t_step_mosq, yM_temp, k4_yM, theta, POP);


	//////////////////////////
	// output

	for (int g = 0; g < N_spec; g++)
	{
		for (int k = 0; k < N_M_comp; k++)
		{
			yM[g][k] = yM[g][k] + t_step_mosq*k1_yM[g][k] / 6.0 + t_step_mosq*k2_yM[g][k] / 3.0 + t_step_mosq*k3_yM[g][k] / 3.0 + t_step_mosq*k4_yM[g][k] / 6.0;
		}
	}

}


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  2.3. Update step for mosquitoes                                          //
//                                                                           // 
//       For every human step we take mosq_steps (=10) steps for mosquitoes  //
//       The smaller step size ensures that the ODE solver works smoothly.   //
//       Especially an issue for the larval stages                           //   
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void mosquito_step(double t, params* theta, population* POP)
{
	//////////////////////////////////
	// Set up mosquito state vector

	double yM[N_spec][N_M_comp];

	for (int g = 0; g < N_spec; g++)
	{
		for (int k = 0; k<N_M_comp; k++)
		{
			yM[g][k] = POP->yM[g][k];
		}
	}


	double t_step_mosq = (double(t_step)) / (double(mosq_steps));


	//////////////////////////////////
	// Force of infection on mosquitoes

	for (int g = 0; g < N_spec; g++)
	{
		theta->lam_M[g] = 0.0;
	}

	for (int n = 0; n < POP->N_pop; n++)
	{
		for (int g = 0; g < N_spec; g++)
		{
			theta->lam_M[g] = theta->lam_M[g] + POP->lam_n[n][g] * (theta->c_PCR*POP->people[n].I_PCR + theta->c_LM*POP->people[n].I_LM +
				              theta->c_D*POP->people[n].I_D + theta->c_T*POP->people[n].T);
		}
	}


	//////////////////////////////////////
	// Carry out the mosq_steps

	for (int j = 0; j<mosq_steps; j++)
	{
		mosq_rk4(t, t_step_mosq, yM, theta, POP);

		for (int g = 0; g < N_spec; g++)
		{
			theta->lam_S_M_track[g].push_back(theta->lam_M[g] * yM[g][3]);
			theta->lam_S_M_track[g].erase(theta->lam_S_M_track[g].begin());
		}
	}

	for (int g = 0; g < N_spec; g++)
	{
		for (int k = 0; k < N_M_comp; k++)
		{
			POP->yM[g][k] = yM[g][k];
		}
	}

}


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.4. Update the vector of human classes                                 //
//                                                                          // 
//       THINK CAREFULLY ABOUT THE ORDERING OF EVENTS                       //
//////////////////////////////////////////////////////////////////////////////

void human_step(params* theta, population* POP)
{

	//////////////////////////////////////////////////////////////////////////
	// 2.4.1 Temporary objects for setting up individuals' intervention   
	//       access characteristics 

	float GMN_parm[(N_int)*(N_int + 3) / 2 + 1];
	float GMN_work[N_int];
	float GMN_zero[N_int];
	float zz_GMN[N_int];

	for (int k = 0; k<N_int; k++)
	{
		GMN_zero[k] = 0.0;
	}


	///////////////////////////////////////////////
	// 2.4.2. Apply ageing

	for (int n = 0; n<POP->N_pop; n++)
	{
		POP->people[n].ager(*theta);
	}


	///////////////////////////////////////////////
	// 2.4.3. Deaths
	//
	// Look again at how things are erased from vectors.

	int N_dead = 0;

	for (int n = 0; n<POP->people.size(); n++)
	{
		/////////////////////////////////////////////
		// Everyone has an equal probability of dying

		if (theta->P_dead > genunf(0, 1))
		{
			POP->people.erase(POP->people.begin() + n);

			POP->pi_n.erase(POP->pi_n.begin() + n);
			POP->lam_n.erase(POP->lam_n.begin() + n);

			N_dead = N_dead + 1;
			n = n - 1;      // If we erase something, the next one moves into it's place so we don't want to step forward.
		}
		else {

			///////////////////////////////////////////
			// People die once they reach the maximum age

			if (POP->people[n].age > theta->age_max)
			{
				POP->people.erase(POP->people.begin() + n);

				POP->pi_n.erase(POP->pi_n.begin() + n);
				POP->lam_n.erase(POP->lam_n.begin() + n);

				N_dead = N_dead + 1;
				n = n - 1;       // If we erase something, the next one moves into it's place so we don't want to step forward.
			}
		}
	}

	/////////////////////////////////////////////////////////////
	// 2.4.4. Births - set up to ensure balanced population.
	//        Can be adjusted to account for changing demography.

	double zeta_start, het_dif_track, q_rand;

	vector<double> zero_push(N_spec);
	for (int g = 0; g < N_spec; g++)
	{
		zero_push[g] = 0.0;
	}


	for (int n = 0; n<N_dead; n++)
	{
		zeta_start = exp(gennor(-0.5*theta->sig_het*theta->sig_het, theta->sig_het));

		while (zeta_start > theta->het_max)
		{
			zeta_start = exp(gennor(-0.5*theta->sig_het*theta->sig_het, theta->sig_het));
		}

		individual HH(0.0, zeta_start);

		HH.S = 1;
		HH.I_PCR = 0;
		HH.I_LM = 0;
		HH.I_D = 0;
		HH.T = 0;
		HH.P = 0;


		HH.A_par = 0.0;
		HH.A_clin = 0.0;

		HH.A_par_boost = 0;
		HH.A_clin_boost = 0;

		HH.A_par_timer = -1.0;
		HH.A_clin_timer = -1.0;

		HH.PQ_proph = 0;
		HH.PQ_proph_timer = -1.0;

		HH.Hyp = 0;

		if (genunf(0.0, 1.0) < 0.5)
		{
			HH.gender = 0;
		}
		else {
			HH.gender = 1;
		}

		if (HH.gender == 0)
		{
			if (genunf(0.0, 1.0) < theta->G6PD_prev)
			{
				HH.G6PD_def = 1;
			}
			else {
				HH.G6PD_def = 0;
			}
		}
		else {

			q_rand = genunf(0.0, 1.0);

			if (q_rand <= theta->G6PD_prev*theta->G6PD_prev)
			{
				HH.G6PD_def = 2;
			}

			if ((q_rand > theta->G6PD_prev*theta->G6PD_prev) && (q_rand <= theta->G6PD_prev*theta->G6PD_prev + 2 * theta->G6PD_prev*(1.0 - theta->G6PD_prev)))
			{
				HH.G6PD_def = 1;
			}

			if (q_rand > (theta->G6PD_prev*theta->G6PD_prev + 2 * theta->G6PD_prev*(1.0 - theta->G6PD_prev)))
			{
				HH.G6PD_def = 0;
			}
		}

		if (genunf(0.0, 1.0) < theta->CYP2D6_prev)
		{
			HH.CYP2D6 = 1;
		}
		else {
			HH.CYP2D6 = 0;
		}

		HH.preg_age = 0;
		HH.pregnant = 0;
		HH.preg_timer = 0.0;


		/////////////////////////////////
		// Assign levels of maternally-acquired immunity
		// by finding women of child-bearing age with the
		// closest level of heterogeneity

		HH.A_par_mat = 0.0;
		HH.A_clin_mat = 0.0;

		het_dif_track = 1e10;

		for (int j = 0; j<POP->people.size(); j++)
		{
			if (POP->people[j].preg_age == 1)
			{
				if (abs(HH.zeta_het - POP->people[j].zeta_het) < het_dif_track)
				{
					HH.A_par_mat = theta->P_mat*POP->people[j].A_par_mat;
					HH.A_clin_mat = theta->P_mat*POP->people[j].A_clin_mat;

					het_dif_track = (HH.zeta_het - POP->people[j].zeta_het)*(HH.zeta_het - POP->people[j].zeta_het);
				}
			}
		}


		///////////////////////////////////////////////////
		// Lagged exposure equals zero - they're not born yet!

		for (int k = 0; k<theta->H_track; k++)
		{
			HH.lam_bite_track.push_back(0.0);
			HH.lam_rel_track.push_back(0.0);
		}


		///////////////////////////////////////////////////
		// Assign intervention access scores

		for (int p = 0; p<N_int; p++)
		{
			for (int q = 0; q<N_int; q++)
			{
				theta->V_int_dummy[p][q] = theta->V_int[p][q];
			}
		}

		setgmn(GMN_zero, *theta->V_int_dummy, N_int, GMN_parm);

		genmn(GMN_parm, zz_GMN, GMN_work);

		for (int k = 0; k<N_int; k++)
		{
			HH.zz_int[k] = zz_GMN[k];
		}


		///////////////////////////////////////////////////
		// Born with no interventions

		HH.LLIN = 0;
		HH.IRS = 0;

		for (int g = 0; g < N_spec; g++)
		{
			HH.w_VC[g] = 1.0;
			HH.y_VC[g] = 1.0;
			HH.z_VC[g] = 0.0;
		}


		/////////////////////////////////////////////////////////////////
		// 2.4.5. Push the created individual onto the vector of people

		POP->people.push_back(HH);

		POP->pi_n.push_back(zero_push);
		POP->lam_n.push_back(zero_push);
	}



	///////////////////////////////////////////////////
	// 2.4.6. Update individual-level vector control

	for (int n = 0; n<POP->N_pop; n++)
	{
		POP->people[n].intervention_updater(*theta);
	}


	///////////////////////////////////////////////////
	// 2.4.7. Update proportion of bites
	//
	//        Note the ordering of n and g loops. Need to 
	//        check if this makes a difference for speed.
	//
	//        Should be able to make this quicker


	for (int n = 0; n<POP->N_pop; n++)
	{
		for (int g = 0; g < N_spec; g++)
		{
			POP->pi_n[n][g] = POP->people[n].zeta_het*(1.0 - theta->rho_age*exp(-POP->people[n].age*theta->age_0_inv));

			//POP->pi_n[n][g] = POP->people[n].zeta_het - (POP->people[n].zeta_het - POP->people[n].zeta_het)*POP->P_age_bite;   // Slightly quicker - no calling of exponentials
		}
	}

	double SIGMA_PI[N_spec];
	for (int g = 0; g < N_spec; g++)
	{
		SIGMA_PI[g] = 0.0;
	}

	for (int n = 0; n < POP->N_pop; n++)
	{
		for (int g = 0; g < N_spec; g++)
		{
			SIGMA_PI[g] = SIGMA_PI[g] + POP->pi_n[n][g];
		}
	}

	for (int g = 0; g < N_spec; g++)
	{
		SIGMA_PI[g] = 1.0 / SIGMA_PI[g];
	}

	for (int n = 0; n < POP->N_pop; n++)
	{
		for (int g = 0; g < N_spec; g++)
		{
			POP->pi_n[n][g] = POP->pi_n[n][g] * SIGMA_PI[g];
		}
	}


	///////////////////////////////////////////////////
	// 2.4.8 Update population-level vector control quantities

	for (int g = 0; g < N_spec; g++)
	{
		POP->SUM_pi_w[g] = 0;
	}

	for (int n = 0; n < POP->N_pop; n++)
	{
		for (int g = 0; g < N_spec; g++)
		{
			POP->SUM_pi_w[g] = POP->SUM_pi_w[g] + POP->pi_n[n][g] * POP->people[n].w_VC[g];
		}
	}


	for (int g = 0; g < N_spec; g++)
	{
		POP->W_VC[g] = 1.0 - theta->Q_0[g] + theta->Q_0[g] * POP->SUM_pi_w[g];
		POP->Z_VC[g] = theta->Q_0[g] * POP->SUM_pi_z[g];

		POP->delta_1_VC[g] = theta->delta_1 / (1.0 - POP->Z_VC[g]);
		POP->delta_VC[g] = POP->delta_1_VC[g] + theta->delta_2;

		POP->p_1_VC[g] = theta->p_1[g] * POP->W_VC[g] / (1.0 - POP->Z_VC[g] * theta->p_1[g]);

		POP->mu_M_VC[g] = -log(POP->p_1_VC[g] * theta->p_2[g]) / POP->delta_VC[g];

		POP->Q_VC[g] = 1.0 - (1.0 - theta->Q_0[g]) / POP->W_VC[g];

		POP->aa_VC[g] = POP->Q_VC[g] / POP->delta_VC[g];

		POP->exp_muM_tauM_VC[g] = exp(-POP->mu_M_VC[g] * theta->tau_M[g]);
		POP->beta_VC[g] = theta->eps_max[g] * POP->mu_M_VC[g] / (exp(POP->delta_VC[g] * POP->mu_M_VC[g]) - 1.0);
	}


	///////////////////////////////////////////////////
	// 2.4.9. Update individual-level force of infection on humans

	for (int n = 0; n < POP->N_pop; n++)
	{
		for (int g = 0; g < N_spec; g++)
		{
			POP->lam_n[n][g] = POP->aa_VC[g] * POP->pi_n[n][g] * POP->people[n].w_VC[g] / POP->SUM_pi_w[g];
		}
	}


	///////////////////////////////////////////////////
	// 2.4.10. Implement moves between compartments
	//
	// TO DO: Can take some multiplications out of the loop.

	double lam_bite_base[N_spec];
	double lam_bite_n;     // better notation (this is lam_bite)

	for (int g = 0; g < N_spec; g++)
	{
		lam_bite_base[g] = (double(POP->N_pop))*theta->bb*POP->yM[g][5];
	}

	for (int n = 0; n<POP->N_pop; n++)
	{
		lam_bite_n = 0.0;

		for (int g = 0; g < N_spec; g++)
		{
			lam_bite_n = lam_bite_n + POP->lam_n[n][g] * lam_bite_base[g];
		}

		POP->people[n].state_mover(*theta, lam_bite_n);
	}

}


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.5. Summarise the output from the population                           //
//                                                                          // 
//////////////////////////////////////////////////////////////////////////////

void POP_summary(population* POP, simulation* SIM)
{
	for (int k = 0; k<N_H_comp; k++)
	{
		POP->yH[k] = 0.0;
	}

	for (int k = 0; k<10; k++)
	{
		POP->prev_all[k] = 0.0;
		POP->prev_U5[k] = 0.0;
		POP->prev_U10[k] = 0.0;
	}


	for (int n = 0; n<POP->N_pop; n++)
	{
		////////////////////////////////////////
		// Numbers in each compartment

		POP->yH[0] = POP->yH[0] + POP->people[n].S;
		POP->yH[1] = POP->yH[1] + POP->people[n].I_PCR;
		POP->yH[2] = POP->yH[2] + POP->people[n].I_LM;
		POP->yH[3] = POP->yH[3] + POP->people[n].I_D;
		POP->yH[4] = POP->yH[4] + POP->people[n].T;
		POP->yH[5] = POP->yH[5] + POP->people[n].P;


		//////////////////////////////////////////////
		//////////////////////////////////////////////
		// Summary - full population

		////////////////////////////////////////
		// Prevalence

		POP->prev_all[0] = POP->prev_all[0] + 1;                                                                        // Numbers - denominator
		POP->prev_all[1] = POP->prev_all[1] + POP->people[n].I_PCR + POP->people[n].I_LM +
			                                + POP->people[n].I_D + POP->people[n].T;                                      // PCR detectable infections
		POP->prev_all[2] = POP->prev_all[2] + POP->people[n].I_LM + POP->people[n].I_D + POP->people[n].T;                // LM detectable infections
		POP->prev_all[3] = POP->prev_all[3] + POP->people[n].I_D + POP->people[n].T;                                      // Clinical episodes

		if (POP->people[n].Hyp > 0)
		{
			POP->prev_all[4] = POP->prev_all[4] + 1;                     // Hypnozoite positive

			POP->prev_all[5] = POP->prev_all[5] + POP->people[n].Hyp;    // Number of batches of hypnozoites
		}


		////////////////////////////////////////
		// Incidence

		POP->prev_all[6]  = POP->prev_all[6]  + POP->people[n].I_PCR_new;
		POP->prev_all[7]  = POP->prev_all[7]  + POP->people[n].I_LM_new;
		POP->prev_all[8]  = POP->prev_all[8]  + POP->people[n].I_D_new;
		POP->prev_all[9]  = POP->prev_all[9]  + POP->people[n].ACT_treat;
		POP->prev_all[10] = POP->prev_all[10] + POP->people[n].PQ_treat;


		//////////////////////////////////////////////
		//////////////////////////////////////////////
		// Summary - under 5's

		if (POP->people[n].age < 1825.0)
		{
			////////////////////////////////////////
			// Prevalence

			POP->prev_U5[0] = POP->prev_U5[0] + 1;                                                                // Numbers - denominator
			POP->prev_U5[1] = POP->prev_U5[1] + POP->people[n].I_PCR + POP->people[n].I_LM
				                              + POP->people[n].I_D + POP->people[n].T;                              // PCR detectable infections
			POP->prev_U5[2] = POP->prev_U5[2] + POP->people[n].I_LM + POP->people[n].I_D + POP->people[n].T;        // LM detectable infections
			POP->prev_U5[3] = POP->prev_U5[3] + POP->people[n].I_D + POP->people[n].T;                              // Clinical episodes

			if (POP->people[n].Hyp > 0)
			{
				POP->prev_U5[4] = POP->prev_U5[4] + 1;                     // Hypnozoite positive

				POP->prev_U5[5] = POP->prev_U5[5] + POP->people[n].Hyp;    // Number of batches of hypnozoites
			}


			////////////////////////////////////////
			// Incidence

			POP->prev_U5[6] = POP->prev_U5[6] + POP->people[n].I_PCR_new;
			POP->prev_U5[7] = POP->prev_U5[7] + POP->people[n].I_LM_new;
			POP->prev_U5[8] = POP->prev_U5[8] + POP->people[n].I_D_new;
			POP->prev_U5[9] = POP->prev_U5[9] + POP->people[n].ACT_treat;
			POP->prev_U5[10] = POP->prev_U5[10] + POP->people[n].PQ_treat;
		}

		//////////////////////////////////////////////
		//////////////////////////////////////////////
		// Summary - under 10's

		if (POP->people[n].age < 3650.0)
		{
			////////////////////////////////////////
			// Prevalence

			POP->prev_U10[0] = POP->prev_U10[0] + 1;                                                            // Numbers - denominator
			POP->prev_U10[1] = POP->prev_U10[1] + POP->people[n].I_PCR + POP->people[n].I_LM
				                                + POP->people[n].I_D + POP->people[n].T;                          // PCR detectable infections
			POP->prev_U10[2] = POP->prev_U10[2] + POP->people[n].I_LM + POP->people[n].I_D + POP->people[n].T;    // LM detectable infections
			POP->prev_U10[3] = POP->prev_U10[3] + POP->people[n].I_D + POP->people[n].T;                          // Clinical episodes

			if (POP->people[n].Hyp > 0)
			{
				POP->prev_U10[4] = POP->prev_U10[4] + 1;                     // Hypnozoite positive

				POP->prev_U10[5] = POP->prev_U10[5] + POP->people[n].Hyp;    // Number of batches of hypnozoites
			}


			////////////////////////////////////////
			// Incidence

			POP->prev_U10[6]  = POP->prev_U10[6]  + POP->people[n].I_PCR_new;
			POP->prev_U10[7]  = POP->prev_U10[7]  + POP->people[n].I_LM_new;
			POP->prev_U10[8]  = POP->prev_U10[8]  + POP->people[n].I_D_new;
			POP->prev_U10[9]  = POP->prev_U10[9]  + POP->people[n].ACT_treat;
			POP->prev_U10[10] = POP->prev_U10[10] + POP->people[n].PQ_treat;
		}
	}


	//////////////////////////////
	// Intervention coverage

	POP->LLIN_cov_t = 0;
	POP->IRS_cov_t = 0;
	POP->ACT_treat_t = 0;
	POP->PQ_treat_t = 0;
	POP->pregnant_t = 0;

	POP->PQ_overtreat_t = 0;
	POP->PQ_overtreat_9m_t = 0;


	for (int n = 0; n<POP->N_pop; n++)
	{
		POP->LLIN_cov_t  = POP->LLIN_cov_t  + POP->people[n].LLIN;
		POP->IRS_cov_t   = POP->IRS_cov_t   + POP->people[n].IRS;
		POP->ACT_treat_t = POP->ACT_treat_t + POP->people[n].ACT_treat;
		POP->PQ_treat_t  = POP->PQ_treat_t  + POP->people[n].PQ_treat;
		POP->pregnant_t  = POP->pregnant_t  + POP->people[n].pregnant;

		POP->PQ_overtreat_t    = POP->PQ_overtreat_t    + POP->people[n].PQ_overtreat;
		POP->PQ_overtreat_9m_t = POP->PQ_overtreat_9m_t + POP->people[n].PQ_overtreat_9m;
	}


	//////////////////////////////
	// Immunity

	double A_par_mean = 0.0, A_clin_mean = 0.0;

	for (int n = 0; n<POP->N_pop; n++)
	{
		A_par_mean = A_par_mean + POP->people[n].A_par;
		A_clin_mean = A_clin_mean + POP->people[n].A_clin;
	}

	POP->A_par_mean_t = A_par_mean / ((double)POP->N_pop);
	POP->A_clin_mean_t = A_clin_mean / ((double)POP->N_pop);
}



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.6. Simulate the model and store the output in SIM                     //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void model_simulator(params* theta, population* POP, intervention* INTVEN, simulation* SIM)
{

	for (int i = 0; i<SIM->N_time; i++)
	{
		if (SIM->t_vec[i] / 365.0 - floor(SIM->t_vec[i] / 365.0) < 0.5*t_step / 365.0)
		{
			cout << "time = " << SIM->t_vec[i] / 365.0 << "\t" << 100.0*(SIM->t_vec[i] - SIM->t_vec[0]) / (double(t_step*SIM->N_time)) << "% complete" << endl;
		}

		human_step(theta, POP);

		mosquito_step(SIM->t_vec[i], theta, POP);

		intervention_dist(SIM->t_vec[i], theta, POP, INTVEN);

		POP_summary(POP, SIM);

		//////////////////////////////////////
		// Fill out simulation object

		for (int k = 0; k<N_H_comp; k++)
		{
			SIM->yH_t[i][k] = POP->yH[k];
		}

		for (int k = 0; k<N_M_comp; k++)
		{
			for (int g = 0; g < N_spec; g++)
			{
				SIM->yM_t[i][g][k] = POP->yM[g][k];
			}
		}

		for (int k = 0; k<11; k++)
		{
			SIM->prev_all[i][k] = POP->prev_all[k];
			SIM->prev_U5[i][k] = POP->prev_U5[k];
			SIM->prev_U10[i][k] = POP->prev_U10[k];
		}


		SIM->LLIN_cov_t[i] = POP->LLIN_cov_t;
		SIM->IRS_cov_t[i] = POP->IRS_cov_t;
		SIM->ACT_treat_t[i] = POP->ACT_treat_t;
		SIM->PQ_treat_t[i] = POP->PQ_treat_t;
		SIM->pregnant_t[i] = POP->pregnant_t;

		SIM->PQ_overtreat_t[i] = POP->PQ_overtreat_t;
		SIM->PQ_overtreat_9m_t[i] = POP->PQ_overtreat_9m_t;


		SIM->EIR_t[i] = 0.0;
		for (int g = 0; g < N_spec; g++)
		{
			SIM->EIR_t[i] = SIM->EIR_t[i] + POP->aa_VC[g] * POP->yM[g][5];
		}

		SIM->A_par_mean_t[i] = POP->A_par_mean_t;
		SIM->A_clin_mean_t[i] = POP->A_clin_mean_t;

	}

}



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.7. Vector control distribution                                        //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////


void intervention_dist(double t, params* theta, population* POP, intervention* INTVEN)
{
	double QQ;

	bool BS_effective;
	bool PQ_treat;
	bool PQ_effective;

	bool MSAT_pos;
	bool SSAT_pos;

	//////////////////////////////////////////////////////////
	// Intervention 1: LLINS

	for (int m = 0; m<INTVEN->LLIN_year.size(); m++)
	{
		if ((t > INTVEN->LLIN_year[m] - 0.5*t_step) &&
			(t < INTVEN->LLIN_year[m] + 0.51*t_step))
		{
			cout << "LLIN distribution" << endl;

			QQ = phi_inv(INTVEN->LLIN_cover[m], 0.0, sqrt(1.0 + theta->sig_round_LLIN*theta->sig_round_LLIN));

			for (int n = 0; n<POP->N_pop; n++)
			{
				if (gennor(POP->people[n].zz_int[0], theta->sig_round_LLIN) < QQ)
				{
					POP->people[n].LLIN = 1;
					POP->people[n].LLIN_age = 0.0;

					for (int g = 0; g < N_spec; g++)
					{
						POP->people[n].d_LLIN[g] = theta->d_LLIN_0[g];
						POP->people[n].r_LLIN[g] = theta->r_LLIN_0[g];
						POP->people[n].s_LLIN[g] = 1.0 - POP->people[n].d_LLIN[g] - POP->people[n].r_LLIN[g];
					}
				}
			}
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 2: IRS

	for (int m = 0; m<INTVEN->IRS_year.size(); m++)
	{
		if ((t > INTVEN->IRS_year[m] - 0.5*t_step) &&
			(t < INTVEN->IRS_year[m] + 0.51*t_step))
		{
			cout << "IRS distribution" << endl;

			QQ = phi_inv(INTVEN->IRS_cover[m], 0.0, sqrt(1.0 + theta->sig_round_IRS*theta->sig_round_IRS));

			for (int n = 0; n<POP->N_pop; n++)
			{
				if (gennor(POP->people[n].zz_int[1], theta->sig_round_IRS) < QQ)
				{
					POP->people[n].IRS = 1;
					POP->people[n].IRS_age = 0.0;

					for (int g = 0; g < N_spec; g++)
					{
						POP->people[n].d_IRS[g] = theta->d_IRS_0[g];
						POP->people[n].r_IRS[g] = theta->r_IRS_0[g];
						POP->people[n].s_IRS[g] = 1.0 - POP->people[n].d_IRS[g] - POP->people[n].r_IRS[g];
					}
				}
			}
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 3: first-line treatment (blood-stage)

	for (int m = 0; m<INTVEN->BS_treat_year_on.size(); m++)
	{
		if ((t > INTVEN->BS_treat_year_on[m] - 0.5*t_step) &&
			(t < INTVEN->BS_treat_year_on[m] + 0.51*t_step))
		{
			cout << "New front-line BS treatment" << endl;

			theta->BS_treat_BScover = INTVEN->BS_treat_BScover[m];
			theta->BS_treat_BSeff   = INTVEN->BS_treat_BSeff[m];
			theta->BS_treat_BSproph = INTVEN->BS_treat_BSproph[m];

			theta->treat_BScover = theta->BS_treat_BScover;
			theta->treat_BSeff   = theta->BS_treat_BSeff;
			theta->treat_PQavail = 0.0;
			theta->r_P           = 1.0 / theta->BS_treat_BSproph;
		}
	}


	//////////////////////////////////////////////////////////
	// Switching back to baseline.

	for (int m = 0; m<INTVEN->BS_treat_year_on.size(); m++)
	{
		if ((t > INTVEN->BS_treat_year_off[m] - 0.5*t_step) &&
			(t < INTVEN->BS_treat_year_off[m] + 0.51*t_step))
		{
			cout << "End of changing front-line BS treatment" << endl;

			theta->treat_BScover = theta->BS_treat_BScover_base;
			theta->treat_BSeff   = theta->BS_treat_BSeff_base;
			theta->treat_PQavail = 0.0;
			theta->r_P           = 1.0 / theta->BS_treat_BSproph_base;
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 4: first-line treatment (primaquine)

	for (int m = 0; m<INTVEN->PQ_treat_year_on.size(); m++)
	{
		if ((t > INTVEN->PQ_treat_year_on[m] - 0.5*t_step) &&
			(t < INTVEN->PQ_treat_year_on[m] + 0.51*t_step))
		{
			cout << "New front-line PQ treatment" << endl;

			theta->PQ_treat_BScover     = INTVEN->PQ_treat_BScover[m];
			theta->PQ_treat_BSeff       = INTVEN->PQ_treat_BSeff[m];
			theta->PQ_treat_BSproph     = INTVEN->PQ_treat_BSproph[m];
			theta->PQ_treat_PQavail     = INTVEN->PQ_treat_PQavail[m];
			theta->PQ_treat_PQeff       = INTVEN->PQ_treat_PQeff[m];
			theta->PQ_treat_PQproph     = INTVEN->PQ_treat_PQproph[m];
			theta->PQ_treat_G6PD_risk   = INTVEN->PQ_treat_G6PD_risk[m];
			theta->PQ_treat_CYP2D6_risk = INTVEN->PQ_treat_CYP2D6_risk[m];
			theta->PQ_treat_preg_risk   = INTVEN->PQ_treat_preg_risk[m];
			theta->PQ_treat_low_age     = INTVEN->PQ_treat_low_age[m];

			theta->treat_BScover = theta->PQ_treat_BScover;
			theta->treat_BSeff   = theta->PQ_treat_BSeff;
			theta->treat_PQavail = theta->PQ_treat_PQavail;
			theta->r_P           = 1.0 / theta->PQ_treat_BSproph;
		}
	}


	//////////////////////////////////////////////////////////
	// Switching back to baseline.

	for (int m = 0; m<INTVEN->PQ_treat_year_on.size(); m++)
	{
		if ((t > INTVEN->PQ_treat_year_off[m] - 0.5*t_step) &&
			(t < INTVEN->PQ_treat_year_off[m] + 0.51*t_step))
		{
			cout << "End of changing front-line PQ treatment" << endl;

			theta->PQ_treat_BScover     = 0.0;
			theta->PQ_treat_BSeff       = 0.0;
			theta->PQ_treat_BSproph     = 10.0;
			theta->PQ_treat_PQavail     = 0.0;
			theta->PQ_treat_PQeff       = 0.0;
			theta->PQ_treat_PQproph     = 10.0;
			theta->PQ_treat_G6PD_risk   = 1;
			theta->PQ_treat_CYP2D6_risk = 1;
			theta->PQ_treat_preg_risk   = 1;
			theta->PQ_treat_low_age     = 180.0;

			theta->treat_BScover = theta->BS_treat_BScover_base;
			theta->treat_BSeff   = theta->BS_treat_BSeff_base;
			theta->treat_PQavail = 0.0;
			theta->r_P           = 1.0 / theta->BS_treat_BSproph_base;
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 5: MDA (blood-stage)

	for (int m = 0; m<INTVEN->MDA_BS_year.size(); m++)
	{
		if ((t > INTVEN->MDA_BS_year[m] - 0.5*t_step) &&
			(t < INTVEN->MDA_BS_year[m] + 0.51*t_step))
		{
			cout << "MDA (BS) distribution" << endl;

			theta->MDA_BS_BScover = INTVEN->MDA_BS_BScover[m];
			theta->MDA_BS_BSeff   = INTVEN->MDA_BS_BSeff[m];
			theta->MDA_BS_BSproph = INTVEN->MDA_BS_BSproph[m];

			QQ = phi_inv(theta->MDA_BS_BScover, 0.0, sqrt(1.0 + theta->sig_round_MDA*theta->sig_round_MDA));

			for (int n = 0; n<POP->N_pop; n++)
			{
				if (gennor(POP->people[n].zz_int[2], theta->sig_round_MDA) < QQ)
				{
					POP->people[n].ACT_treat = 1;

					if (genunf(0.0, 1.0) < theta->MDA_BS_BSeff)
					{
						if (POP->people[n].S == 1    ) { POP->people[n].S = 0;     POP->people[n].P = 1; }
						if (POP->people[n].I_PCR == 1) { POP->people[n].I_PCR = 0; POP->people[n].P = 1; }
						if (POP->people[n].I_LM == 1 ) { POP->people[n].I_LM = 0;  POP->people[n].P = 1; }
						if (POP->people[n].I_D == 1  ) { POP->people[n].I_D = 0;   POP->people[n].T = 1; }
					}
				}
			}
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 6: MDA (blood-stage and liver-stage)

	for (int m = 0; m<INTVEN->MDA_PQ_year.size(); m++)
	{
		if ((t > INTVEN->MDA_PQ_year[m] - 0.5*t_step) &&
			(t < INTVEN->MDA_PQ_year[m] + 0.51*t_step))
		{
			cout << "MDA (BS+PQ) distribution" << endl;

			theta->MDA_PQ_BScover     = INTVEN->MDA_PQ_BScover[m];
			theta->MDA_PQ_BSeff       = INTVEN->MDA_PQ_BSeff[m];
			theta->MDA_PQ_BSproph     = INTVEN->MDA_PQ_BSproph[m];
			theta->MDA_PQ_PQavail     = INTVEN->MDA_PQ_PQavail[m];
			theta->MDA_PQ_PQeff       = INTVEN->MDA_PQ_PQeff[m];
			theta->MDA_PQ_PQproph     = INTVEN->MDA_PQ_PQproph[m];
			theta->MDA_PQ_G6PD_risk   = INTVEN->MDA_PQ_G6PD_risk[m];
			theta->MDA_PQ_CYP2D6_risk = INTVEN->MDA_PQ_CYP2D6_risk[m];
			theta->MDA_PQ_preg_risk   = INTVEN->MDA_PQ_preg_risk[m];
			theta->MDA_PQ_low_age     = INTVEN->MDA_PQ_low_age[m];

			QQ = phi_inv(theta->MDA_PQ_BScover, 0.0, sqrt(1.0 + theta->sig_round_MDA*theta->sig_round_MDA));

			for (int n = 0; n<POP->N_pop; n++)
			{
				if (gennor(POP->people[n].zz_int[4], theta->sig_round_MDA) < QQ)
				{
					/////////////////////////////////////////////////////
					// Blood-stage treatment is always administered
					// Is blood-stage treatment effective

					BS_effective = 0;

					if (genunf(0.0, 1.0) < theta->MDA_PQ_BSeff)
					{
						BS_effective = 1;
					}


					/////////////////////////////////////////////////////////////////////
					// Is PQ administered?

					PQ_treat = 0;

					if( genunf(0.0, 1.0) < theta->MDA_PQ_PQavail )
					{
						PQ_treat = 1;
					}
					

					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of G6PD deficiency

					if( (theta->MDA_PQ_G6PD_risk == 1) && (POP->people[n].G6PD_def == 1) )
					{
						PQ_treat = 0;
					}


					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of pregancy

					if( (theta->MDA_PQ_preg_risk == 1) && (POP->people[n].pregnant == 1) )
					{
						PQ_treat = 0;
					}


					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of young age

					if (POP->people[n].age < theta->MDA_PQ_low_age)
					{
						PQ_treat = 0;
					}

					if (PQ_treat == 1)
					{
						POP->people[n].PQ_treat = 1;
					}

					/////////////////////////////////////////////////////////////////////
					// Is PQ effective?

					PQ_effective = 1;

					if( genunf(0.0, 1.0) > theta->MDA_PQ_PQeff )
					{
						PQ_effective = 0;
					}

					if( (theta->MDA_PQ_CYP2D6_risk == 1) && (POP->people[n].CYP2D6 == 1) )
					{
						PQ_effective = 0;
					}


					/////////////////////////////////////////////////////////////////////
					// Was there PQ overtreatment?

					if( (PQ_treat == 1) && (POP->people[n].Hyp == 0) )
					{
						POP->people[n].PQ_overtreat = 1;
					}

					if( (PQ_treat == 1) && (POP->people[n].T_last_BS > 270.0) )
					{
						POP->people[n].PQ_overtreat_9m = 1;
					}


					/////////////////////////////////////////////////////////////////////
					// ACTION: administer blood-stage drug

					POP->people[n].ACT_treat = 1;

					if (BS_effective == 1)
					{
						if (POP->people[n].S == 1) {     POP->people[n].S = 0;     POP->people[n].P = 1; }
						if (POP->people[n].I_PCR == 1) { POP->people[n].I_PCR = 0; POP->people[n].P = 1; }
						if (POP->people[n].I_LM == 1) {  POP->people[n].I_LM = 0;  POP->people[n].P = 1; }
						if (POP->people[n].I_D == 1) {   POP->people[n].I_D = 0;   POP->people[n].T = 1; }
					}


					/////////////////////////////////////////////////////////////////////
					// ACTION: administer primaquine

					if ((PQ_treat == 1) && (PQ_effective == 1))
					{
						POP->people[n].Hyp = 0;

						POP->people[n].PQ_proph = 1;
						POP->people[n].PQ_proph_timer = theta->SSAT_PQ_PQproph;
					}
				}
			}
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 7: MSAT (blood-stage and liver-stage)

	for (int m = 0; m < INTVEN->MSAT_PQ_year.size(); m++)
	{
		if( (t > INTVEN->MSAT_PQ_year[m] - 0.5*t_step) &&
			(t < INTVEN->MSAT_PQ_year[m] + 0.51*t_step) )
		{
			cout << "MSAT (BS+PQ) distribution" << endl;

			theta->MSAT_PQ_BScover     = INTVEN->MSAT_PQ_BScover[m];
			theta->MSAT_PQ_RDT_PCR     = INTVEN->MSAT_PQ_RDT_PCR[m];
			theta->MSAT_PQ_sens        = INTVEN->MSAT_PQ_sens[m];
			theta->MSAT_PQ_BSeff       = INTVEN->MSAT_PQ_BSeff[m];
			theta->MSAT_PQ_BSproph     = INTVEN->MSAT_PQ_BSproph[m];
			theta->MSAT_PQ_PQavail     = INTVEN->MSAT_PQ_PQavail[m];
			theta->MSAT_PQ_PQeff       = INTVEN->MSAT_PQ_PQeff[m];
			theta->MSAT_PQ_PQproph     = INTVEN->MSAT_PQ_PQproph[m];
			theta->MSAT_PQ_G6PD_risk   = INTVEN->MSAT_PQ_G6PD_risk[m];
			theta->MSAT_PQ_CYP2D6_risk = INTVEN->MSAT_PQ_CYP2D6_risk[m];
			theta->MSAT_PQ_preg_risk   = INTVEN->MSAT_PQ_preg_risk[m];
			theta->MSAT_PQ_low_age     = INTVEN->MSAT_PQ_low_age[m];

			QQ = phi_inv(theta->MSAT_PQ_BScover, 0.0, sqrt(1.0 + theta->sig_round_MDA*theta->sig_round_MDA));

			for (int n = 0; n < POP->N_pop; n++)
			{
				if (gennor(POP->people[n].zz_int[5], theta->sig_round_MDA) < QQ)
				{
					/////////////////////////////////////////////////////
					// Blood-stage treatment is always administered

					MSAT_pos = 0;

					////////////////////////////////////////////////
					// Diagnosis by RDT, assumed the same as LM 	

					if (theta->MSAT_PQ_RDT_PCR == 1)
					{
						if ((POP->people[n].I_LM == 1) || (POP->people[n].I_D == 1) || (POP->people[n].T == 1))
						{
							if (genunf(0.0, 1.0) < theta->MSAT_PQ_sens)
							{
								MSAT_pos = 1;
							}
						}
					}


					////////////////////////////////////////////////
					// Diagnosis by PCR 	

					if (theta->MSAT_PQ_RDT_PCR == 2)
					{
						if ((POP->people[n].I_PCR == 1) || (POP->people[n].I_LM == 1) || (POP->people[n].I_D == 1) || (POP->people[n].T == 1))
						{
							if (genunf(0.0, 1.0) < theta->MSAT_PQ_sens)
							{
								MSAT_pos = 1;
							}
						}
					}


					/////////////////////////////////////////////////////
					// Is blood-stage treatment effective

					BS_effective = 0;

					if (genunf(0.0, 1.0) < theta->MSAT_PQ_BSeff)
					{
						BS_effective = 1;
					}


					/////////////////////////////////////////////////////////////////////
					// Is PQ administered?

					PQ_treat = 0;

					if (MSAT_pos == 1)
					{
						if (genunf(0.0, 1.0) < theta->MSAT_PQ_PQavail)
						{
							PQ_treat = 1;
						}
					}


					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of G6PD deficiency

					if ( (theta->MSAT_PQ_G6PD_risk == 1) && (POP->people[n].G6PD_def == 1) )
					{
						PQ_treat = 0;
					}

					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of pregancy

					if ( (theta->MSAT_PQ_preg_risk == 1) && (POP->people[n].pregnant == 1) )
					{
						PQ_treat = 0;
					}

					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of young age

					if (POP->people[n].age < theta->MSAT_PQ_low_age)
					{
						PQ_treat = 0;
					}

					if (PQ_treat == 1)
					{
						POP->people[n].PQ_treat = 1;
					}

					/////////////////////////////////////////////////////////////////////
					// Is PQ effective?

					PQ_effective = 1;

					if (genunf(0.0, 1.0) > theta->MSAT_PQ_PQeff)
					{
						PQ_effective = 0;
					}

					if ((theta->MSAT_PQ_CYP2D6_risk == 1) && (POP->people[n].CYP2D6 == 1))
					{
						PQ_effective = 0;
					}


					/////////////////////////////////////////////////////////////////////
					// Was there PQ overtreatment?

					if ((PQ_treat == 1) && (POP->people[n].Hyp == 0))
					{
						POP->people[n].PQ_overtreat = 1;
					}

					if ((PQ_treat == 1) && (POP->people[n].T_last_BS > 270.0))
					{
						POP->people[n].PQ_overtreat_9m = 1;
					}


					/////////////////////////////////////////////////////////////////////
					// ACTION: administer blood-stage drug

					if (MSAT_pos == 1)
					{
						POP->people[n].ACT_treat = 1;

						if (BS_effective == 1)
						{
							if (POP->people[n].S == 1) {     POP->people[n].S = 0;     POP->people[n].P = 1; }
							if (POP->people[n].I_PCR == 1) { POP->people[n].I_PCR = 0; POP->people[n].P = 1; }
							if (POP->people[n].I_LM == 1) {  POP->people[n].I_LM = 0;  POP->people[n].P = 1; }
							if (POP->people[n].I_D == 1) {   POP->people[n].I_D = 0;   POP->people[n].T = 1; }
						}
					}

					
					/////////////////////////////////////////////////////////////////////
					// ACTION: administer primaquine

					if ((PQ_treat == 1) && (PQ_effective == 1))
					{
						POP->people[n].Hyp = 0;

						POP->people[n].PQ_proph = 1;
						POP->people[n].PQ_proph_timer = theta->SSAT_PQ_PQproph;
					}
				}
			}
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 8: SSAT (blood-stage and liver-stage)

	for (int m = 0; m<INTVEN->SSAT_PQ_year.size(); m++)
	{
		if ((t > INTVEN->SSAT_PQ_year[m] - 0.5*t_step) &&
			(t < INTVEN->SSAT_PQ_year[m] + 0.51*t_step))
		{
			cout << "SSAT (BS+PQ) distribution" << endl;

			theta->SSAT_PQ_BScover     = INTVEN->SSAT_PQ_BScover[m];
			theta->SSAT_PQ_sens        = INTVEN->SSAT_PQ_sens[m];
			theta->SSAT_PQ_spec        = INTVEN->SSAT_PQ_spec[m];
			theta->SSAT_PQ_BSeff       = INTVEN->SSAT_PQ_BSeff[m];
			theta->SSAT_PQ_BSproph     = INTVEN->SSAT_PQ_BSproph[m];
			theta->SSAT_PQ_PQavail     = INTVEN->SSAT_PQ_PQavail[m];
			theta->SSAT_PQ_PQeff       = INTVEN->SSAT_PQ_PQeff[m];
			theta->SSAT_PQ_PQproph     = INTVEN->SSAT_PQ_PQproph[m];
			theta->SSAT_PQ_G6PD_risk   = INTVEN->SSAT_PQ_G6PD_risk[m];
			theta->SSAT_PQ_CYP2D6_risk = INTVEN->SSAT_PQ_CYP2D6_risk[m];
			theta->SSAT_PQ_preg_risk   = INTVEN->SSAT_PQ_preg_risk[m];
			theta->SSAT_PQ_low_age     = INTVEN->SSAT_PQ_low_age[m];

			QQ = phi_inv(theta->SSAT_PQ_BScover, 0.0, sqrt(1.0 + theta->sig_round_MDA*theta->sig_round_MDA));

			for (int n = 0; n < POP->N_pop; n++)
			{
				if (gennor(POP->people[n].zz_int[6], theta->sig_round_MDA) < QQ)
				{
					/////////////////////////////////////////////////////
					// Blood-stage treatment is always administered

					POP->people[n].ACT_treat = 1;


					/////////////////////////////////////////////////////
					// Is blood-stage treatment effective

					BS_effective = 0;

					if( genunf(0.0, 1.0) < theta->SSAT_PQ_BSeff )
					{
						BS_effective = 1;
					}

					/////////////////////////////////////////////////////////////////////
					// SSAT screening for blood-stage infection in the last 9 months

					SSAT_pos = 0;

					if( (POP->people[n].T_last_BS <= 270.0) && (genunf(0.0, 1.0) < theta->SSAT_PQ_sens) )
					{
						SSAT_pos = 1;
					}

					if( (POP->people[n].T_last_BS > 270.0) && (genunf(0.0, 1.0) > theta->SSAT_PQ_spec) )
					{
						SSAT_pos = 1;
					}

					/////////////////////////////////////////////////////////////////////
					// Is PQ administered?

					PQ_treat = 0;

					if( SSAT_pos == 1 )
					{
						if( genunf(0.0, 1.0) < theta->SSAT_PQ_PQavail )
						{
							PQ_treat = 1;
						}
					}


					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of G6PD deficiency

					if( (theta->SSAT_PQ_G6PD_risk == 1) && (POP->people[n].G6PD_def == 1) )
					{
						PQ_treat = 0;
					}

					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of pregancy

					if( (theta->SSAT_PQ_preg_risk == 1) && (POP->people[n].pregnant == 1) )
					{
						PQ_treat = 0;
					}

					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of young age

					if( POP->people[n].age < theta->SSAT_PQ_low_age )
					{
						PQ_treat = 0;
					}

					if( PQ_treat == 1 )
					{
						POP->people[n].PQ_treat = 1;
					}

					/////////////////////////////////////////////////////////////////////
					// Is PQ effective?

					PQ_effective = 0;

					if( genunf(0.0, 1.0) < theta->SSAT_PQ_PQeff )
					{
						PQ_effective = 1;
					}

					if( (theta->SSAT_PQ_CYP2D6_risk == 1) && (POP->people[n].CYP2D6 == 1) )
					{
						PQ_effective = 0;
					}


					/////////////////////////////////////////////////////////////////////
					// Was there PQ overtreatment?

					if ((PQ_treat == 1) && (POP->people[n].Hyp == 0))
					{
						POP->people[n].PQ_overtreat = 1;
					}

					if ((PQ_treat == 1) && (POP->people[n].T_last_BS > 270.0))
					{
						POP->people[n].PQ_overtreat_9m = 1;
					}


					/////////////////////////////////////////////////////////////////////
					// ACTION: administer blood-stage drug

					if( BS_effective == 1 )
					{
						if (POP->people[n].S     == 1) { POP->people[n].S = 0;     POP->people[n].P = 1; }
						if (POP->people[n].I_PCR == 1) { POP->people[n].I_PCR = 0; POP->people[n].P = 1; }
						if (POP->people[n].I_LM  == 1) { POP->people[n].I_LM = 0;  POP->people[n].P = 1; }
						if (POP->people[n].I_D   == 1) { POP->people[n].I_D = 0;   POP->people[n].T = 1; }
					}

					/////////////////////////////////////////////////////////////////////
					// ACTION: administer primaquine

					if( (PQ_treat == 1) && (PQ_effective == 1) )
					{
						POP->people[n].Hyp = 0;

						POP->people[n].PQ_proph = 1;
						POP->people[n].PQ_proph_timer = theta->SSAT_PQ_PQproph;
					}
				}
			}


		}
	}

}


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.8. Competing hazards sampler                                          //
//                                                                          //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

int CH_sample(double *xx, int nn)
{
	vector<double> xx_cum(nn);

	xx_cum[0] = xx[0];

	for (int k = 1; k<nn; k++)
	{
		xx_cum[k] = xx_cum[k - 1] + xx[k];
	}

	int index = 0;
	double unif = genunf(0, 1);

	if (unif < xx_cum[0])
	{
		return index;
	}

	for (int k = 1; k<nn; k++)
	{
		if ((unif > xx_cum[k - 1]) & (unif < xx_cum[k]))
		{
			index = k;

			return index;
		}
	}

	return index;
}


//////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////// 
//        //                                                                            // 
//  2.9.  //  Inverse of the cumulative normal distribution function required           //
//        //  for implementing correlated intervention coverage.                        //
////////////                                                                            //
////////////  The code is based on the following website                                //
////////////  http://www.johndcook.com/blog/cpp_phi_inverse/                            //
////////////  which is based the algorithm in Abramowitz and Stegun formula 26.2.23.    // 
////////////                                                                            //
////////////  The absolute value of the error should be less than 4.5 e-4 and it tests  //
////////////  out nicely in R.                                                          //
////////////                                                                            //
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

double phi_inv(double pp, double mu, double sigma)
{
	if (pp <= 0.0 || pp >= 1.0)
	{
		throw("bad vlaue of pp (coverage) in erfinv");
	}

	double cc[] = { 2.515517, 0.802853, 0.010328 };
	double dd[] = { 1.432788, 0.189269, 0.001308 };

	double tt, temp;

	if (pp < 0.5)
	{
		tt = sqrt(-2.0*log(pp));

		temp = -(tt - ((cc[2] * tt + cc[1])*tt + cc[0]) / (((dd[2] * tt + dd[1])*tt + dd[0])*tt + 1.0));

	}
	else {

		tt = sqrt(-2.0*log(1.0 - pp));

		temp = (tt - ((cc[2] * tt + cc[1])*tt + cc[0]) / (((dd[2] * tt + dd[1])*tt + dd[0])*tt + 1.0));
	}

	return mu + sigma*temp;
}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//         //                                                        //
//  2.10.  //  Log gamma function, based on gamma.h from NRC3        //
//         //                                                        //
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

double gammln(const double xx)
{
	int j;
	double x, tmp, y, ser;
	static const double cof[14] = { 57.1562356658629235, -59.5979603554754912,
		14.1360979747417471, -0.491913816097620199, 0.339946499848118887e-4,
		0.465236289270485756e-4, -0.983744753048795646e-4, 0.158088703224912494e-3,
		-0.210264441724104883e-3, 0.217439618115212643e-3, -0.164318106536763890e-3,
		0.844182239838527433e-4, -0.261908384015814087e-4, 0.368991826595316234e-5 };
	if (xx <= 0) throw("bad arg in gammln");
	y = x = xx;
	tmp = x + 5.242187500000000;
	tmp = (x + 0.5)*log(tmp) - tmp;
	ser = 0.999999999999997092;
	for (j = 0; j<14; j++) ser += cof[j] / ++y;
	return tmp + log(2.5066282746310005*ser / x);
}


/////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////// 
//          //                                                         // 
//   ####   //  #####  ####  ##  ## #### ##      #####   ####  #####   //
//  ##  ##  //  ##    ##  ## ##  ##  ##  ##      ##  ## ##  ## ##  ##  //
//     ##   //  ####  ##  ## ##  ##  ##  ##      #####  ##  ## #####   //
//  ##  ##  //  ##     ####  ##  ##  ##  ##      ##     ##  ## ##      //
//   ####   //  #####     ##  ####  #### #####   ##      ####  ##      //
//          //                                                         //
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
//       //                                                                                           //
//  3.1  //  LU decomposition of a matrix                                                             // 
//       //  Based on ludcmp.cpp from Numerical Recipes in C++                                        //
//       //                                                                                           //
////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                    //
//  Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise       //
//  permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;    //
//  indx[1..n] is an output vector that records the row permutation effected by the partial           //
//  pivoting; d is output as ±1 depending on whether the number of row interchanges was even          //
//  or odd, respectively. This routine is used in combination with lubksb to solve linear equations   //
//  or invert a matrix                                                                                //
//                                                                                                    // 
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

void ludcmp(vector<vector<double>> &a, int n_dim, vector<int> &indx, double &d)
{
	const double TINY = 1.0e-20;
	int i, imax, j, k;
	double big, dum, sum, temp;


	vector<double> vv(n_dim);
	d = 1.0;
	for (i = 0; i<n_dim; i++) {
		big = 0.0;
		for (j = 0; j<n_dim; j++)
			if ((temp = fabs(a[i][j])) > big) big = temp;
		if (big == 0.0) throw("Singular matrix in routine ludcmp");
		vv[i] = 1.0 / big;
	}
	for (j = 0; j<n_dim; j++) {
		for (i = 0; i<j; i++) {
			sum = a[i][j];
			for (k = 0; k<i; k++) sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for (i = j; i<n_dim; i++) {
			sum = a[i][j];
			for (k = 0; k<j; k++) sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			if ((dum = vv[i] * fabs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 0; k<n_dim; k++) {
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			d = -d;
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a[j][j] == 0.0) a[j][j] = TINY;
		if (j != n_dim - 1) {
			dum = 1.0 / (a[j][j]);
			for (i = j + 1; i<n_dim; i++) a[i][j] *= dum;
		}
	}
}


///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
//       //                                                      // 
//  3.2  //  Matrix back substitution                            //
//       //  Based on lubksb.cpp from Numerical Recipes in C++   //
//       //                                                      //
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

void lubksb(vector<vector<double>> &a, int n_dim, vector<int> &indx, vector<double> &b)
{
	int i, ii = 0, ip, j;
	double sum;


	for (i = 0; i<n_dim; i++) {
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii != 0)
			for (j = ii - 1; j<i; j++) sum -= a[i][j] * b[j];
		else if (sum != 0.0)
			ii = i + 1;
		b[i] = sum;
	}
	for (i = n_dim - 1; i >= 0; i--) {
		sum = b[i];
		for (j = i + 1; j<n_dim; j++) sum -= a[i][j] * b[j];
		b[i] = sum / a[i][i];
	}
}


///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
//       //                                                      //
//  3.3  //  Matrix inversion and calculation of determinant.    //
//       //  Based on ludcmp.cpp from Numerical Recipes in C++   //
//       //                                                      //
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

void matrix_inv(vector<vector<double>> &a, int n, vector<vector<double>> &a_inv)
{
	vector<int> a_index(n);
	vector<double> col(n);
	double d;

	ludcmp(a, n, a_index, d);

	for (int j = 0; j<n; j++)
	{
		for (int i = 0; i<n; i++)
		{
			col[i] = 0.0;
		}
		col[j] = 1.0;

		lubksb(a, n, a_index, col);

		for (int i = 0; i<n; i++)
		{
			a_inv[i][j] = col[i];
		}
	}

}


///////////////////////////////////////////////////
///////////////////////////////////////////////////
//       //                                      // 
//  3.4  //  Matrix multiplication               //
//       //  Calculates inv(MM*)xx               //
//       //                                      //
///////////////////////////////////////////////////
///////////////////////////////////////////////////

void inv_MM_bb(vector<vector<double>> &MM, vector<double> &bb, vector<double> &xx, int n_dim)
{
	///////////////////////////////////////////////
	// 3.4.1. calculate inv(MM)

	vector<vector<double>> MM_inv;
	MM_inv.resize(n_dim);
	for (int k = 0; k < n_dim; k++)
	{
		MM_inv[k].resize(n_dim);
	}

	matrix_inv(MM, n_dim, MM_inv);


	///////////////////////////////////////////////
	// 3.4.2. calculate xx = MM_inv*bb

	for (int i = 0; i<n_dim; i++)
	{
		xx[i] = 0.0;

		for (int j = 0; j<n_dim; j++)
		{
			xx[i] = xx[i] + MM_inv[i][j] * bb[j];
		}
	}

}


///////////////////////////////////////////////////
///////////////////////////////////////////////////
//       //                                      // 
//  3.5  //  Equilibrium matrix                  //
//       //                                      //
///////////////////////////////////////////////////
///////////////////////////////////////////////////

void MM_ij(int i, int j, params* theta, population* POP, vector<vector<double>> &MM,
	vector<vector<double>> lam_eq, vector<vector<vector<double>>> phi_LM_eq,
	vector<vector<vector<double>>> phi_D_eq, vector<vector<vector<double>>> r_PCR_eq)
{
	//////////////////////////////////////////////
	// 4.5.1. Initialise matrix with zeros

	for (int i1 = 0; i1<(N_H_comp * (K_max + 1)); i1++)
	{
		for (int j1 = 0; j1<(N_H_comp * (K_max + 1)); j1++)
		{
			MM[i1][j1] = 0.0;
		}
	}


	//////////////////////////////////////////////
	// 4.5.2. Fill out non-zero elements

	for (int k1 = 0; k1<(K_max + 1); k1++)
	{
		for (int k2 = 0; k2<(K_max + 1); k2++)
		{
			MM[0 * (K_max + 1) + k1][0 * (K_max + 1) + k2] = - lam_eq[i][j]*theta->D_MAT[k1][k2] - theta->ff*theta->K_MAT[k1][k2]
				                                             + theta->gamma_L*theta->L_MAT[k1][k2] - theta->mu_H*theta->D_MAT[k1][k2] - POP->r_age[i] * theta->D_MAT[k1][k2];
			MM[0 * (K_max + 1) + k1][1 * (K_max + 1) + k2] = + r_PCR_eq[i][j][k2]*theta->D_MAT[k1][k2];
			MM[0 * (K_max + 1) + k1][5 * (K_max + 1) + k2] = +theta->r_P*theta->D_MAT[k1][k2];

			MM[1 * (K_max + 1) + k1][0 * (K_max + 1) + k2] = + lam_eq[i][j]*(1.0 - phi_LM_eq[i][j][k2])*theta->OD_MAT[k1][k2] + theta->ff*(1.0 - phi_LM_eq[i][j][k2])*theta->K_MAT[k1][k2];
			MM[1 * (K_max + 1) + k1][1 * (K_max + 1) + k2] = - lam_eq[i][j]*theta->D_MAT[k1][k2] - theta->ff*theta->K_MAT[k1][k2] - r_PCR_eq[i][j][k2]*theta->D_MAT[k1][k2]
				                                             + lam_eq[i][j]*(1.0 - phi_LM_eq[i][j][k2])*theta->OD_MAT[k1][k2] + theta->ff*(1.0 - phi_LM_eq[i][j][k2])*theta->K_MAT[k1][k2]
				                                             + theta->gamma_L*theta->L_MAT[k1][k2] - theta->mu_H*theta->D_MAT[k1][k2] - POP->r_age[i] * theta->D_MAT[k1][k2];
			MM[1 * (K_max + 1) + k1][2 * (K_max + 1) + k2] = + theta->r_LM*theta->D_MAT[k1][k2];

			MM[2 * (K_max + 1) + k1][0 * (K_max + 1) + k2] = + lam_eq[i][j]*phi_LM_eq[i][j][k2]*(1.0 - phi_D_eq[i][j][k2])*theta->OD_MAT[k1][k2] + theta->ff*phi_LM_eq[i][j][k2]*(1.0 - phi_D_eq[i][j][k2])*theta->K_MAT[k1][k2];
			MM[2 * (K_max + 1) + k1][1 * (K_max + 1) + k2] = + lam_eq[i][j]*phi_LM_eq[i][j][k2]*(1.0 - phi_D_eq[i][j][k2])*theta->OD_MAT[k1][k2] + theta->ff*phi_LM_eq[i][j][k2] * (1.0 - phi_D_eq[i][j][k2])*theta->K_MAT[k1][k2];
			MM[2 * (K_max + 1) + k1][2 * (K_max + 1) + k2] = - lam_eq[i][j]*theta->D_MAT[k1][k2] - theta->ff*theta->K_MAT[k1][k2] - theta->r_LM*theta->D_MAT[k1][k2]
				                                             + lam_eq[i][j]*(1.0 - phi_D_eq[i][j][k2])*theta->OD_MAT[k1][k2] + theta->ff*(1.0 - phi_D_eq[i][j][k2])*theta->K_MAT[k1][k2]
				                                             + theta->gamma_L*theta->L_MAT[k1][k2] - theta->mu_H*theta->D_MAT[k1][k2] - POP->r_age[i] * theta->D_MAT[k1][k2];
			MM[2 * (K_max + 1) + k1][3 * (K_max + 1) + k2] = + theta->r_D*theta->D_MAT[k1][k2];

			MM[3 * (K_max + 1) + k1][0 * (K_max + 1) + k2] = + lam_eq[i][j]*phi_LM_eq[i][j][k2]*phi_D_eq[i][j][k2]*(1.0 - theta->treat_BScover*theta->treat_BSeff)*theta->OD_MAT[k1][k2] 
				                                             + theta->ff*phi_LM_eq[i][j][k2]*phi_D_eq[i][j][k2]*(1.0 - theta->treat_BScover*theta->treat_BSeff)*theta->K_MAT[k1][k2];
			MM[3 * (K_max + 1) + k1][1 * (K_max + 1) + k2] = + lam_eq[i][j]*phi_LM_eq[i][j][k2]*phi_D_eq[i][j][k2] * (1.0 - theta->treat_BScover*theta->treat_BSeff)*theta->OD_MAT[k1][k2] 
				                                             + theta->ff*phi_LM_eq[i][j][k2]*phi_D_eq[i][j][k2]*(1.0 - theta->treat_BScover*theta->treat_BSeff)*theta->K_MAT[k1][k2];
			MM[3 * (K_max + 1) + k1][2 * (K_max + 1) + k2] = + lam_eq[i][j]*phi_D_eq[i][j][k2]*(1.0 - theta->treat_BScover*theta->treat_BSeff)*theta->OD_MAT[k1][k2] + theta->ff*phi_D_eq[i][j][k2]*(1.0 - theta->treat_BScover*theta->treat_BSeff)*theta->K_MAT[k1][k2];
			MM[3 * (K_max + 1) + k1][3 * (K_max + 1) + k2] = - lam_eq[i][j]*theta->D_MAT[k1][k2] - theta->r_D*theta->D_MAT[k1][k2] + lam_eq[i][j]*theta->OD_MAT[k1][k2]
				                                             + theta->gamma_L*theta->L_MAT[k1][k2] - theta->mu_H*theta->D_MAT[k1][k2] - POP->r_age[i]*theta->D_MAT[k1][k2];

			MM[4 * (K_max + 1) + k1][0 * (K_max + 1) + k2] = + lam_eq[i][j]*phi_LM_eq[i][j][k2]*phi_D_eq[i][j][k2]*theta->treat_BScover*theta->treat_BSeff*theta->OD_MAT[k1][k2] 
				                                             + theta->ff*phi_LM_eq[i][j][k2]*phi_D_eq[i][j][k2]*theta->treat_BScover*theta->treat_BSeff*theta->K_MAT[k1][k2];
			MM[4 * (K_max + 1) + k1][1 * (K_max + 1) + k2] = + lam_eq[i][j]*phi_LM_eq[i][j][k2]*phi_D_eq[i][j][k2]*theta->treat_BScover*theta->treat_BSeff*theta->OD_MAT[k1][k2] 
				                                             + theta->ff*phi_LM_eq[i][j][k2]*phi_D_eq[i][j][k2]*theta->treat_BScover*theta->treat_BSeff*theta->K_MAT[k1][k2];
			MM[4 * (K_max + 1) + k1][2 * (K_max + 1) + k2] = + lam_eq[i][j]*phi_D_eq[i][j][k2]*theta->treat_BScover*theta->treat_BSeff*theta->OD_MAT[k1][k2] + theta->ff*phi_D_eq[i][j][k2]*theta->treat_BScover*theta->treat_BSeff*theta->K_MAT[k1][k2];
			MM[4 * (K_max + 1) + k1][4 * (K_max + 1) + k2] = - lam_eq[i][j]*theta->D_MAT[k1][k2] - theta->r_T*theta->D_MAT[k1][k2] + lam_eq[i][j] * theta->OD_MAT[k1][k2]
				                                             + theta->gamma_L*theta->L_MAT[k1][k2] - theta->mu_H*theta->D_MAT[k1][k2] - POP->r_age[i] * theta->D_MAT[k1][k2];

			MM[5 * (K_max + 1) + k1][4 * (K_max + 1) + k2] = + theta->r_T*theta->D_MAT[k1][k2];
			MM[5 * (K_max + 1) + k1][5 * (K_max + 1) + k2] = - lam_eq[i][j]*theta->D_MAT[k1][k2] - theta->r_P*theta->D_MAT[k1][k2] + lam_eq[i][j] * theta->OD_MAT[k1][k2]
				                                             + theta->gamma_L*theta->L_MAT[k1][k2] - theta->mu_H*theta->D_MAT[k1][k2] - POP->r_age[i] * theta->D_MAT[k1][k2];
		}
	}

}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//       //                                                 //
//  3.6  //  Guass-Hermite weights for Gaussian quadrature  //
//       //  integration with Normal distribution.          //
//       //  Code is adapted from gauher.cpp from NR3       // 
//       //                                                 //
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void gauher(population* POP, params* theta)
{
	double x[N_het];
	double w[N_het];

	/////////////////////////////
	// PART 1

	const double EPS = 1.0e-14, PIM4 = 0.7511255444649425;
	const int MAXIT = 10;
	int i, its, j, m;
	double p1, p2, p3, pp, z, z1;

	m = (N_het + 1) / 2;
	for (i = 0; i<m; i++) {
		if (i == 0) {
			z = sqrt(double(2 * N_het + 1)) - 1.85575*pow(double(2 * N_het + 1), -0.16667);
		}
		else if (i == 1) {
			z -= 1.14*pow(double(N_het), 0.426) / z;
		}
		else if (i == 2) {
			z = 1.86*z - 0.86*x[0];
		}
		else if (i == 3) {
			z = 1.91*z - 0.91*x[1];
		}
		else {
			z = 2.0*z - x[i - 2];
		}
		for (its = 0; its<MAXIT; its++) {
			p1 = PIM4;
			p2 = 0.0;
			for (j = 0; j<N_het; j++) {
				p3 = p2;
				p2 = p1;
				p1 = z*sqrt(2.0 / (j + 1))*p2 - sqrt(double(j) / (j + 1))*p3;
			}
			pp = sqrt(double(2 * N_het))*p2;
			z1 = z;
			z = z1 - p1 / pp;
			if (fabs(z - z1) <= EPS) break;
		}
		if (its >= MAXIT) throw("too many iterations in gauher");
		x[i] = z;
		x[N_het - 1 - i] = -z;
		w[i] = 2.0 / (pp*pp);
		w[N_het - 1 - i] = w[i];
	}


	/////////////////////////////
	// PART 2

	double w_sum = 0.0;

	for (int j = 0; j<N_het; j++)
	{
		w_sum = w_sum + w[j];
	}

	////////////////////////////////
	// Note that N_het-1-j here instead of j to ensure x_het increasing

	for (int j = 0; j<N_het; j++)
	{
		POP->x_het[j] = exp(theta->sig_het*x[N_het - 1 - j] * sqrt(2.0) - 0.5*theta->sig_het*theta->sig_het);
		POP->w_het[j] = w[j] / w_sum;
	}


	////////////////////////////
	// temporary for N_het = 1

	if (N_het == 1)
	{
		POP->x_het[0] = 1.0;
		POP->w_het[0] = 1.0;
	}

	for (int i = 0; i<N_age; i++)
	{
		for (int j = 0; j<N_het; j++)
		{
			POP->x_age_het[i][j] = POP->age_bite[i] * POP->x_het[j];
			POP->w_age_het[i][j] = POP->age_demog[i] * POP->w_het[j];
		}
	}

	////////////////////////////////
	// Boundaries of heterogeneity compartments

	POP->x_het_bounds[0] = 0.0;

	for (int j = 1; j<N_het; j++)
	{
		POP->x_het_bounds[j] = exp(0.5*(log(POP->x_het[j - 1]) + log(POP->x_het[j])));
	}

	POP->x_het_bounds[N_het] = theta->het_max;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//       //                                                                   // 
//  3.7  //  Initialise a population of individuals at equilibrium            //
//       //  This calculates the non-seasonal equilibrium solution from a     //
//       //  comparable deterministic model. This will give an approximately  // 
//       //  correct solution for a seasonal setting.                         //           
//       //                                                                   //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void equi_pop_setup(population* POP, params* theta)
{
	//////////////////////////////////////////////////////
	// 3.7.1. Set up age and heterogeneity compartments

	////////////////////////////////////////
	// 3.7.1.1. Bounds of age bins

	//	double age_bounds[N_age+1] = {0.0*365.0, 20.0*365.0, 40.0*365.0, 60.0*365.0, 80.0*365.0};

	double age_bounds[N_age + 1] = { 0.0*365.0, 0.2*365.0, 0.4*365.0, 0.6*365.0, 0.8*365.0, 1.0*365.0,
		                             1.2*365.0, 1.4*365.0, 1.6*365.0, 1.8*365.0, 2.0*365.0,
		                             2.2*365.0, 2.4*365.0, 2.6*365.0, 2.8*365.0, 3.0*365.0,
		                             3.4*365.0, 3.8*365.0, 4.2*365.0, 4.6*365.0, 5.0*365.0,
		                             5.5*365.0, 6.0*365.0, 6.5*365.0, 7.0*365.0, 7.5*365.0, 8.0*365.0, 8.5*365.0, 9.0*365.0, 9.5*365.0, 10.0*365.0,
		                             11.0*365.0, 12.0*365.0, 13.0*365.0, 14.0*365.0, 15.0*365.0, 16.0*365.0, 17.0*365.0, 18.0*365.0, 19.0*365.0, 20.0*365.0,
		                             22.0*365.0, 24.0*365.0, 26.0*365.0, 28.0*365.0, 30.0*365.0, 32.0*365.0, 34.0*365.0, 36.0*365.0, 38.0*365.0, 40.0*365.0,
		                             45.0*365.0, 50.0*365.0, 55.0*365.0, 60.0*365.0, 65.0*365.0, 70.0*365.0, 75.0*365.0, 80.0*365.0 };


	////////////////////////////////////////////
	// 3.7.1.2 Proportion in each age bin

	for (int i = 0; i<(N_age - 1); i++)
	{
		POP->age_demog[i] = exp(-theta->mu_H*age_bounds[i]) - exp(-theta->mu_H*age_bounds[i + 1]);
	}

	POP->age_demog[N_age - 1] = 1.0;

	for (int i = 0; i<(N_age - 1); i++)
	{
		POP->age_demog[N_age - 1] = POP->age_demog[N_age - 1] - POP->age_demog[i];
	}


	////////////////////////////////////////
	// 3.7.1.3. Ageing rates - formula below ensures
	//          balanced demography

	POP->r_age[0] = theta->mu_H*(1.0 - POP->age_demog[0]) / POP->age_demog[0];

	for (int i = 1; i<(N_age - 1); i++)
	{
		POP->r_age[i] = (POP->r_age[i - 1] * POP->age_demog[i - 1] - theta->mu_H*POP->age_demog[i]) / POP->age_demog[i];
	}

	POP->r_age[N_age - 1] = 0.0;


	////////////////////////////////////////
	// 3.7.1.4. Age-dependent mosquito biting rates

	for (int i = 0; i<N_age; i++)
	{
		POP->age_mids[i] = 0.5*(age_bounds[i] + age_bounds[i + 1]);
	}

	for (int i = 0; i<N_age; i++)
	{
		POP->age_bite[i] = 1.0 - theta->rho_age*exp(-POP->age_mids[i] / theta->age_0);
	}

	POP->P_age_bite = exp(-t_step / theta->age_0);


	///////////////////////////////////////
	// Ensure total bites are normalised

	double omega_age = 0.0;

	for (int i = 0; i<N_age; i++)
	{
		omega_age = omega_age + POP->age_demog[i] * POP->age_bite[i];
	}

	omega_age = 1 / omega_age;


	for (int i = 0; i<N_age; i++)
	{
		POP->age_bite[i] = omega_age*POP->age_bite[i];
	}


	//////////////////////////////////////////////////////
	// 3.7.1.5. Find age category closest to 20 yr old woman

	POP->index_age_20 = 0;

	double age_diff = (POP->age_mids[0] - 20.0*365.0)*(POP->age_mids[0] - 20.0*365.0);

	for (int i = 1; i<N_age; i++)
	{
		if ((POP->age_mids[i] - 20.0*365.0)*(POP->age_mids[i] - 20.0*365.0) < age_diff)
		{
			age_diff = (POP->age_mids[i] - 20.0*365.0)*(POP->age_mids[i] - 20.0*365.0);
			POP->index_age_20 = i;
		}
	}


	////////////////////////////////////////
	// 3.7.1.6. Heterogeneity in expsoure and age demographics.
	//          Weights for heterogeneity calculate using Gaussian-Hemite quadrature.
	//          The function also fills out the N_age*N_het matrix

	gauher(POP, theta);


	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	// 3.7.2. Calculate equilibrium of model in humans

	//////////////////////////////////////////////////////////////
	// 3.7.2.1. Object for storing equilibrium solution

	vector<vector<vector<vector<double>>>> yH_eq;
	yH_eq.resize(N_age);
	for (int i = 0; i<N_age; i++)
	{
		yH_eq[i].resize(N_het);
		for (int j = 0; j<N_het; j++)
		{
			yH_eq[i][j].resize(K_max + 1);
			for (int k = 0; k < (K_max + 1); k++)
			{
				yH_eq[i][j][k].resize(N_H_comp);
			}
		}
	}


	vector<vector<double>> lam_eq;
	lam_eq.resize(N_age);
	for (int i = 0; i < N_age; i++)
	{
		lam_eq[i].resize(N_het);
	}

	vector<vector<vector<double>>> A_par_eq;
	A_par_eq.resize(N_age);
	for (int i = 0; i < N_age; i++)
	{
		A_par_eq[i].resize(N_het);
		for (int j = 0; j < N_het; j++)
		{
			A_par_eq[i][j].resize(K_max + 1);
		}
	}

	vector<vector<double>> A_par_eq_mean;
	A_par_eq_mean.resize(N_age);
	for (int i = 0; i < N_age; i++)
	{
		A_par_eq_mean[i].resize(N_het);
	}

	vector<vector<vector<double>>> A_clin_eq;
	A_clin_eq.resize(N_age);
	for (int i = 0; i < N_age; i++)
	{
		A_clin_eq[i].resize(N_het);
		for (int j = 0; j < N_het; j++)
		{
			A_clin_eq[i][j].resize(K_max + 1);
		}
	}

	vector<vector<double>> A_clin_eq_mean;
	A_clin_eq_mean.resize(N_age);
	for (int i = 0; i < N_age; i++)
	{
		A_clin_eq_mean[i].resize(N_het);
	}

	vector<vector<vector<double>>> phi_LM_eq;
	phi_LM_eq.resize(N_age);
	for (int i = 0; i < N_age; i++)
	{
		phi_LM_eq[i].resize(N_het);
		for (int j = 0; j < N_het; j++)
		{
			phi_LM_eq[i][j].resize(K_max + 1);
		}
	}

	vector<vector<vector<double>>> phi_D_eq;
	phi_D_eq.resize(N_age);
	for (int i = 0; i < N_age; i++)
	{
		phi_D_eq[i].resize(N_het);
		for (int j = 0; j < N_het; j++)
		{
			phi_D_eq[i][j].resize(K_max + 1);
		}
	}

	vector<vector<vector<double>>> r_PCR_eq;
	r_PCR_eq.resize(N_age);
	for (int i = 0; i < N_age; i++)
	{
		r_PCR_eq[i].resize(N_het);
		for (int j = 0; j < N_het; j++)
		{
			r_PCR_eq[i][j].resize(K_max + 1);
		}
	}


	vector<vector<double>> MM;
	MM.resize(N_H_comp * (K_max + 1));
	for (int l = 0; l < N_H_comp * (K_max + 1); l++)
	{
		MM[l].resize(N_H_comp * (K_max + 1));
	}

	vector<double> bb(N_H_comp * (K_max + 1));
	vector<double> xx(N_H_comp * (K_max + 1));


	//////////////////////////////////////////////////////////////
	// 3.7.2.2. Equilibrium force of infection
	//
	//   Only the contribution from mosquito bites

	for (int i = 0; i<N_age; i++)
	{
		for (int j = 0; j<N_het; j++)
		{
			lam_eq[i][j] = theta->EIR_equil*theta->bb*POP->x_age_het[i][j];
		}
	}


	//////////////////////////////////////////////////////////////
	// 3.7.2.5. Equilibrium number of batches of relapses

	vector<vector<vector<double>>> HH_eq;
	HH_eq.resize(N_age);
	for (int i = 0; i < N_age; i++)
	{
		HH_eq[i].resize(N_het);
		for (int j = 0; j < N_het; j++)
		{
			HH_eq[i][j].resize(K_max + 1);
		}
	}

	vector<double> HH_bb(K_max + 1);
	vector<double> HH_xx(K_max + 1);

	vector<vector<double>> HH_mat;
	HH_mat.resize(K_max + 1);
	for (int k = 0; k < (K_max + 1); k++)
	{
		HH_mat[k].resize(K_max + 1);
	}

	double HH_denom;


	for (int j = 0; j < N_het; j++)
	{
		/////////////////////////////
		// Youngest age category

		for (int k = 0; k < (K_max + 1); k++)
		{
			HH_bb[k] = 0.0;
		}
		HH_bb[0] = POP->w_het[j] * theta->mu_H;


		for (int k1 = 0; k1 < (K_max + 1); k1++)
		{
			for (int k2 = 0; k2 < (K_max + 1); k2++)
			{
				HH_mat[k1][k2] = lam_eq[0][j] * theta->H_MAT[k1][k2] + theta->gamma_L*theta->L_MAT[k1][k2] - theta->mu_H*theta->D_MAT[k1][k2] - POP->r_age[0] * theta->D_MAT[k1][k2];
			}
		}


		inv_MM_bb(HH_mat, HH_bb, HH_xx, K_max + 1);


		for (int k = 0; k < (K_max + 1); k++)
		{
			HH_eq[0][j][k] = -HH_xx[k];
		}


		///////////////////////////////////
		// Older age categories

		for (int i = 1; i < N_age; i++)
		{
			for (int k = 0; k < (K_max + 1); k++)
			{
				HH_bb[k] = -POP->r_age[i - 1] * HH_xx[k];
			}

			for (int k1 = 0; k1 < (K_max + 1); k1++)
			{
				for (int k2 = 0; k2 < (K_max + 1); k2++)
				{
					HH_mat[k1][k2] = lam_eq[i][j] * theta->H_MAT[k1][k2] + theta->gamma_L*theta->L_MAT[k1][k2] - theta->mu_H*theta->D_MAT[k1][k2] - POP->r_age[i] * theta->D_MAT[k1][k2];
				}
			}

			inv_MM_bb(HH_mat, HH_bb, HH_xx, K_max + 1);

			for (int k = 0; k < (K_max + 1); k++)
			{
				HH_eq[i][j][k] = -HH_xx[k];
			}
		}
	}


	//////////////////////////////////////////////////////////////
	// 3.7.2.7. Equilibrium levels of immunity

	double w_HH[K_max + 1];
	vector<double> A_par_vec(K_max + 1);
	vector<double> A_clin_vec(K_max + 1);

	vector<double> ODE_eq_vec(K_max + 1);

	vector<vector<double>> ODE_eq_MAT;
	ODE_eq_MAT.resize(K_max + 1);
	for (int k = 0; k < (K_max + 1); k++)
	{
		ODE_eq_MAT[k].resize(K_max + 1);
	}


	double G_VEC[K_max + 1];
	double LAM_MAT[K_max + 1][K_max + 1];
	double GAM_MAT[K_max + 1][K_max + 1];


	//////////////////////////////////////////////////////////////
	// 3.7.2.7.1. ANTI-PARASITE IMMUNITY

	for (int j = 0; j<N_het; j++)
	{
		/////////////////////////////
		//                         //
		//  Youngest age category  // 
		//                         // 
		/////////////////////////////

		/////////////////////////////
		// Vector part of ODE

		G_VEC[0] = 0.0;
		for (int k = 1; k < (K_max + 1); k++)
		{
			G_VEC[k] = lam_eq[0][j] * HH_eq[0][j][k - 1] / HH_eq[0][j][k] + theta->ff*((double)k);
		}
		G_VEC[K_max] = G_VEC[K_max] + lam_eq[0][j];

		for (int k = 0; k < (K_max + 1); k++)
		{
			G_VEC[k] = G_VEC[k] / (G_VEC[k] * theta->u_par + 1.0);
		}


		for (int k = 0; k < (K_max + 1); k++)
		{
			ODE_eq_vec[k] = G_VEC[k];
		}


		/////////////////////////////
		// Matrix part of ODE

		for (int k1 = 0; k1 < (K_max + 1); k1++)
		{
			for (int k2 = 0; k2 < (K_max + 1); k2++)
			{
				LAM_MAT[k1][k2] = 0.0;
				GAM_MAT[k1][k2] = 0.0;
			}
		}

		for (int k = 0; k<K_max; k++)
		{
			LAM_MAT[k][k] = -1.0;
			LAM_MAT[k + 1][k] = HH_eq[0][j][k] / HH_eq[0][j][k + 1];
		}

		for (int k = 1; k <(K_max + 1); k++)
		{
			GAM_MAT[k][k] = -((double)k);
			GAM_MAT[k - 1][k] = ((double)k)*HH_eq[0][j][k] / HH_eq[0][j][k - 1];
		}

		for (int k1 = 0; k1 < (K_max + 1); k1++)
		{
			for (int k2 = 0; k2 < (K_max + 1); k2++)
			{
				ODE_eq_MAT[k1][k2] = -(lam_eq[0][j] * LAM_MAT[k1][k2] + theta->gamma_L*GAM_MAT[k1][k2] -
					theta->r_par*theta->D_MAT[k1][k2] - theta->mu_H*theta->D_MAT[k1][k2] - POP->r_age[0] * theta->D_MAT[k1][k2]);
			}
		}


		/////////////////////////////
		// Solve matrix equation

		inv_MM_bb(ODE_eq_MAT, ODE_eq_vec, A_par_vec, K_max + 1);

		for (int k = 0; k < (K_max + 1); k++)
		{
			A_par_eq[0][j][k] = A_par_vec[k];
		}


		/////////////////////////////
		//                         //
		//  Older age categories   // 
		//                         // 
		/////////////////////////////

		for (int i = 1; i < N_age; i++)
		{
			/////////////////////////////
			// Vector part of ODE

			G_VEC[0] = 0.0;
			for (int k = 1; k < (K_max + 1); k++)
			{
				G_VEC[k] = lam_eq[i][j] * HH_eq[i][j][k - 1] / HH_eq[i][j][k] + theta->ff*((double)k);
			}
			G_VEC[K_max] = G_VEC[K_max] + lam_eq[i][j];

			for (int k = 0; k < (K_max + 1); k++)
			{
				G_VEC[k] = G_VEC[k] / (G_VEC[k] * theta->u_par + 1.0);
			}


			for (int k = 0; k < (K_max + 1); k++)
			{
				ODE_eq_vec[k] = G_VEC[k] + POP->r_age[i - 1] * A_par_eq[i - 1][j][k] * HH_eq[i - 1][j][k] / HH_eq[i][j][k];
			}


			/////////////////////////////
			// Matrix part of ODE

			for (int k1 = 0; k1 < (K_max + 1); k1++)
			{
				for (int k2 = 0; k2 < (K_max + 1); k2++)
				{
					LAM_MAT[k1][k2] = 0.0;
					GAM_MAT[k1][k2] = 0.0;
				}
			}

			for (int k = 0; k<K_max; k++)
			{
				LAM_MAT[k][k] = -1.0;
				LAM_MAT[k + 1][k] = HH_eq[i][j][k] / HH_eq[i][j][k + 1];
			}

			for (int k = 1; k <(K_max + 1); k++)
			{
				GAM_MAT[k][k] = -((double)k);
				GAM_MAT[k - 1][k] = ((double)k)*HH_eq[i][j][k] / HH_eq[i][j][k - 1];
			}

			for (int k1 = 0; k1 < (K_max + 1); k1++)
			{
				for (int k2 = 0; k2 < (K_max + 1); k2++)
				{
					ODE_eq_MAT[k1][k2] = -(lam_eq[i][j] * LAM_MAT[k1][k2] + theta->gamma_L*GAM_MAT[k1][k2] -
						theta->r_par*theta->D_MAT[k1][k2] - theta->mu_H*theta->D_MAT[k1][k2] - POP->r_age[i] * theta->D_MAT[k1][k2]);
				}
			}


			/////////////////////////////
			// Solve matrix equation

			inv_MM_bb(ODE_eq_MAT, ODE_eq_vec, A_par_vec, K_max + 1);

			for (int k = 0; k < (K_max + 1); k++)
			{
				A_par_eq[i][j][k] = A_par_vec[k];
			}

		}
	}


	//////////////////////////////////////////////////////////////
	// 3.7.2.7.2. MEAN ANTI-PARASITE IMMUNITY

	for (int i = 0; i < N_age; i++)
	{
		for (int j = 0; j < N_het; j++)
		{
			/////////////////////////////
			// Set up weights and rates

			HH_denom = 0.0;

			for (int k = 0; k < (K_max + 1); k++)
			{
				w_HH[k] = HH_eq[i][j][k];
				HH_denom = HH_denom + HH_eq[i][j][k];
			}

			for (int k = 0; k < (K_max + 1); k++)
			{
				w_HH[k] = w_HH[k] / HH_denom;
			}

			/////////////////////////////
			// Average level of immunity

			A_par_eq_mean[i][j] = 0.0;

			for (int k = 0; k < (K_max + 1); k++)
			{
				A_par_eq_mean[i][j] = A_par_eq_mean[i][j] + A_par_eq[i][j][k] * w_HH[k];
			}
		}
	}


	//////////////////////////////////////////////////////////////
	// 3.7.2.7.3. ANTI-CLINICAL IMMUNITY

	for (int j = 0; j<N_het; j++)
	{
		/////////////////////////////
		//                         //
		//  Youngest age category  // 
		//                         // 
		/////////////////////////////

		/////////////////////////////
		// Vector part of ODE

		G_VEC[0] = 0.0;
		for (int k = 1; k < (K_max + 1); k++)
		{
			G_VEC[k] = lam_eq[0][j] * (HH_eq[0][j][k - 1] / HH_eq[0][j][k]) + theta->ff*((double)k);
		}
		G_VEC[K_max] = G_VEC[K_max] + lam_eq[0][j];

		for (int k = 0; k < (K_max + 1); k++)
		{
			G_VEC[k] = G_VEC[k] / (G_VEC[k] * theta->u_clin + 1.0);
		}


		for (int k = 0; k < (K_max + 1); k++)
		{
			ODE_eq_vec[k] = G_VEC[k];
		}


		/////////////////////////////
		// Matrix part of ODE

		for (int k1 = 0; k1 < (K_max + 1); k1++)
		{
			for (int k2 = 0; k2 < (K_max + 1); k2++)
			{
				LAM_MAT[k1][k2] = 0.0;
				GAM_MAT[k1][k2] = 0.0;
			}
		}

		for (int k = 0; k<K_max; k++)
		{
			LAM_MAT[k][k] = -1.0;
			LAM_MAT[k + 1][k] = HH_eq[0][j][k] / HH_eq[0][j][k + 1];
		}

		for (int k = 1; k <(K_max + 1); k++)
		{
			GAM_MAT[k][k] = -((double)k);
			GAM_MAT[k - 1][k] = ((double)k)*HH_eq[0][j][k] / HH_eq[0][j][k - 1];
		}

		for (int k1 = 0; k1 < (K_max + 1); k1++)
		{
			for (int k2 = 0; k2 < (K_max + 1); k2++)
			{
				ODE_eq_MAT[k1][k2] = -(lam_eq[0][j] * LAM_MAT[k1][k2] + theta->gamma_L*GAM_MAT[k1][k2] -
					theta->r_clin*theta->D_MAT[k1][k2] - theta->mu_H*theta->D_MAT[k1][k2] - POP->r_age[0] * theta->D_MAT[k1][k2]);
			}
		}


		/////////////////////////////
		// Solve matrix equation

		inv_MM_bb(ODE_eq_MAT, ODE_eq_vec, A_clin_vec, K_max + 1);

		for (int k = 0; k < (K_max + 1); k++)
		{
			A_clin_eq[0][j][k] = A_clin_vec[k];
		}


		/////////////////////////////
		//                         //
		//  Older age categories   // 
		//                         // 
		/////////////////////////////

		for (int i = 1; i < N_age; i++)
		{
			/////////////////////////////
			// Vector part of ODE

			G_VEC[0] = 0.0;
			for (int k = 1; k < (K_max + 1); k++)
			{
				G_VEC[k] = lam_eq[i][j] * HH_eq[i][j][k - 1] / HH_eq[i][j][k] + theta->ff*((double)k);
			}
			G_VEC[K_max] = G_VEC[K_max] + lam_eq[i][j];

			for (int k = 0; k < (K_max + 1); k++)
			{
				G_VEC[k] = G_VEC[k] / (G_VEC[k] * theta->u_clin + 1.0);
			}


			for (int k = 0; k < (K_max + 1); k++)
			{
				ODE_eq_vec[k] = G_VEC[k] + POP->r_age[i - 1] * A_clin_eq[i - 1][j][k] * HH_eq[i - 1][j][k] / HH_eq[i][j][k];
			}


			/////////////////////////////
			// Matrix part of ODE

			for (int k1 = 0; k1 < (K_max + 1); k1++)
			{
				for (int k2 = 0; k2 < (K_max + 1); k2++)
				{
					LAM_MAT[k1][k2] = 0.0;
					GAM_MAT[k1][k2] = 0.0;
				}
			}

			for (int k = 0; k<K_max; k++)
			{
				LAM_MAT[k][k] = -1.0;
				LAM_MAT[k + 1][k] = HH_eq[i][j][k] / HH_eq[i][j][k + 1];
			}

			for (int k = 1; k <(K_max + 1); k++)
			{
				GAM_MAT[k][k] = -((double)k);
				GAM_MAT[k - 1][k] = ((double)k)*HH_eq[i][j][k] / HH_eq[i][j][k - 1];
			}

			for (int k1 = 0; k1 < (K_max + 1); k1++)
			{
				for (int k2 = 0; k2 < (K_max + 1); k2++)
				{
					ODE_eq_MAT[k1][k2] = -(lam_eq[i][j] * LAM_MAT[k1][k2] + theta->gamma_L*GAM_MAT[k1][k2] -
						theta->r_clin*theta->D_MAT[k1][k2] - theta->mu_H*theta->D_MAT[k1][k2] - POP->r_age[i] * theta->D_MAT[k1][k2]);
				}
			}


			/////////////////////////////
			// Solve matrix equation

			inv_MM_bb(ODE_eq_MAT, ODE_eq_vec, A_clin_vec, K_max + 1);

			for (int k = 0; k < (K_max + 1); k++)
			{
				A_clin_eq[i][j][k] = A_clin_vec[k];
			}

		}
	}


	//////////////////////////////////////////////////////////////
	// 3.7.2.7.4. MEAN ANTI-CLINICAL IMMUNITY

	for (int i = 0; i < N_age; i++)
	{
		for (int j = 0; j < N_het; j++)
		{
			/////////////////////////////
			// Set up weights and rates

			HH_denom = 0.0;

			for (int k = 0; k < (K_max + 1); k++)
			{
				w_HH[k] = HH_eq[i][j][k];
				HH_denom = HH_denom + HH_eq[i][j][k];
			}

			for (int k = 0; k < (K_max + 1); k++)
			{
				w_HH[k] = w_HH[k] / HH_denom;
			}

			/////////////////////////////
			// Average level of immunity

			A_clin_eq_mean[i][j] = 0.0;

			for (int k = 0; k < (K_max + 1); k++)
			{
				A_clin_eq_mean[i][j] = A_clin_eq_mean[i][j] + A_clin_eq[i][j][k] * w_HH[k];
			}
		}
	}


	//////////////////////////////////////////////////////////////
	// 3.7.2.7.5. ADD IN MATERNAL IMMUNITY

	for (int j = 0; j < N_het; j++)
	{
		for (int i = 0; i < N_age; i++)
		{
			for (int k = 0; k < (K_max + 1); k++)
			{
				A_par_eq[i][j][k] = A_par_eq[i][j][k] + A_par_eq_mean[POP->index_age_20][j] * theta->P_mat*exp(-POP->age_mids[i] / theta->d_mat);
				A_clin_eq[i][j][k] = A_clin_eq[i][j][k] + A_clin_eq_mean[POP->index_age_20][j] * theta->P_mat*exp(-POP->age_mids[i] / theta->d_mat);
			}
		}
	}


	//////////////////////////////////////////////////////////////
	// 3.7.2.7.6. EFFECTS OF IMMUNITY

	for (int j = 0; j < N_het; j++)
	{
		for (int i = 0; i < N_age; i++)
		{
			for (int k = 0; k < (K_max + 1); k++)
			{
				phi_LM_eq[i][j][k] = theta->phi_LM_min + (theta->phi_LM_max - theta->phi_LM_min) / (1.0 + pow(A_par_eq[i][j][k] / theta->A_LM_50pc, theta->K_LM));
				phi_D_eq[i][j][k] = theta->phi_D_min + (theta->phi_D_max - theta->phi_D_min) / (1.0 + pow(A_clin_eq[i][j][k] / theta->A_D_50pc, theta->K_D));
				r_PCR_eq[i][j][k] = 1.0 / (theta->d_PCR_min + (theta->d_PCR_max - theta->d_PCR_min) / (1.0 + pow(A_par_eq[i][j][k] / theta->A_PCR_50pc, theta->K_PCR)));
			}
		}
	}


	///////////////////////////////////////////////
	// 3.7.2.8.. Equilibrium states

	for (int j = 0; j<N_het; j++)
	{
		///////////////////////////////////////////////
		// 3.7.2.8.1. Youngest age category

		MM_ij(0, j, theta, POP, MM,
			lam_eq, phi_LM_eq, phi_D_eq, r_PCR_eq);

		for (int k = 0; k<N_H_comp*(K_max + 1); k++)
		{
			bb[k] = 0.0;
		}

		bb[0] = -POP->w_het[j] * theta->mu_H;

		inv_MM_bb(MM, bb, xx, N_H_comp*(K_max + 1));


		/////////////////////////////////////
		// Fill out equilibrium levels

		for (int c = 0; c < N_H_comp; c++)
		{
			for (int k = 0; k < (K_max + 1); k++)
			{
				yH_eq[0][j][k][c] = xx[c*(K_max + 1) + k];
			}
		}


		///////////////////////////////////////////////
		// 3.7.2.8.2. Older age categories

		for (int i = 1; i<N_age; i++)
		{
			MM_ij(i, j, theta, POP, MM,
				lam_eq, phi_LM_eq, phi_D_eq, r_PCR_eq);

			for (int c = 0; c<N_H_comp*(K_max + 1); c++)
			{
				bb[c] = -POP->r_age[i - 1] * xx[c];
			}

			inv_MM_bb(MM, bb, xx, N_H_comp * (K_max + 1));


			/////////////////////////////////////
			// Fill out equilibrium levels

			for (int c = 0; c < N_H_comp; c++)
			{
				for (int k = 0; k < (K_max + 1); k++)
				{
					yH_eq[i][j][k][c] = xx[c*(K_max + 1) + k];
				}
			}

		}
	}

	/*

	cout << "Testing yH sum......" << endl;

	double yH_sum = 0.0;

	for (int i = 0; i < N_age; i++)
	{
	for (int j = 0; j < N_het; j++)
	{
	for (int k = 0; k < (K_max + 1); k++)
	{
	for (int c = 0; c < N_H_comp; c++)
	{
	yH_sum = yH_sum + yH_eq[i][j][k][c];
	}
	}
	}
	}


	cout << "yH_sum = " << yH_sum << endl;
	system("PAUSE");

	*/

	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	// 3.7.3. Calculate equilibrium of model in mosquitoes

	for (int g = 0; g < N_spec; g++)
	{
		theta->lam_M[g] = 0.0;

		for (int i = 0; i<N_age; i++)
		{
			for (int j = 0; j<N_het; j++)
			{
				for (int k = 0; k < (K_max + 1); k++)
				{
					theta->lam_M[g] = theta->lam_M[g] + POP->x_age_het[i][j] * theta->aa[g] * (theta->c_PCR*yH_eq[i][j][k][1] + theta->c_LM*yH_eq[i][j][k][2] +
						theta->c_D*yH_eq[i][j][k][3] + theta->c_T*yH_eq[i][j][k][4]);
				}
			}
		}
	}


	double I_M_star[N_spec];
	for (int g = 0; g < N_spec; g++)
	{
		I_M_star[g] = theta->lam_M[g] * exp(-theta->mu_M[g] * theta->tau_M[g]) / (theta->lam_M[g] + theta->mu_M[g]);
	}

	double a_I_M_sum;

	if (theta->Prop_mosq[0] > 0.0)
	{
		a_I_M_sum = theta->aa[0] * I_M_star[0];

		for (int g = 1; g < N_spec; g++)
		{
			a_I_M_sum = a_I_M_sum + (theta->Prop_mosq[g] / theta->Prop_mosq[0])*theta->aa[g] * I_M_star[g];
		}

		theta->mm_0[0] = theta->EIR_equil / a_I_M_sum;

		for (int g = 1; g < N_spec; g++)
		{
			theta->mm_0[g] = (theta->Prop_mosq[g] / theta->Prop_mosq[0])*theta->mm_0[0];
		}
	}

	if ((theta->Prop_mosq[0] < 1.0e-10) && (theta->Prop_mosq[1] > 0.0))
	{
		theta->mm_0[0] = 0.0;

		a_I_M_sum = theta->aa[1] * I_M_star[1];

		for (int g = 2; g < N_spec; g++)
		{
			a_I_M_sum = a_I_M_sum + (theta->Prop_mosq[g] / theta->Prop_mosq[1])*theta->aa[g] * I_M_star[g];
		}

		theta->mm_0[1] = theta->EIR_equil / a_I_M_sum;

		for (int g = 2; g < N_spec; g++)
		{
			theta->mm_0[g] = (theta->Prop_mosq[g] / theta->Prop_mosq[1])*theta->mm_0[1];
		}
	}


	if ((theta->Prop_mosq[0] < 1.0e-10) && (theta->Prop_mosq[1] < 1.0e-10))
	{
		theta->mm_0[0] = 0.0;
		theta->mm_0[1] = 0.0;
		theta->mm_0[2] = theta->EIR_equil / (theta->aa[2] * I_M_star[2]);

		//theta->mm_0[2] = theta->EIR_equil / (theta->aa[2]*(theta->lam_M[2] / (theta->lam_M[2] + theta->mu_M[2]))*exp(-theta->mu_M[2]*theta->tau_M[2]));
	}



	for (int g = 0; g < N_spec; g++)
	{
		POP->yM[g][0] = 2.0*theta->omega_larvae[g] * theta->mu_M[g] * theta->d_L_larvae*(1.0 + theta->d_pupae*theta->mu_P)*theta->mm_0[g];
		POP->yM[g][1] = 2.0*theta->mu_M[g] * theta->d_L_larvae*(1.0 + theta->d_pupae*theta->mu_P)*theta->mm_0[g];
		POP->yM[g][2] = 2.0*theta->d_pupae*theta->mu_M[g] * theta->mm_0[g];
		POP->yM[g][3] = theta->mm_0[g] * (theta->mu_M[g] / (theta->lam_M[g] + theta->mu_M[g]));
		POP->yM[g][4] = theta->mm_0[g] * (theta->lam_M[g] / (theta->lam_M[g] + theta->mu_M[g]))*(1.0 - exp(-theta->mu_M[g] * theta->tau_M[g]));
		POP->yM[g][5] = theta->mm_0[g] * (theta->lam_M[g] / (theta->lam_M[g] + theta->mu_M[g]))*exp(-theta->mu_M[g] * theta->tau_M[g]);

		theta->Karry[g] = theta->mm_0[g] * 2.0*theta->d_L_larvae*theta->mu_M[g] * (1.0 + theta->d_pupae*theta->mu_P)*theta->gamma_larvae*(theta->omega_larvae[g] + 1.0) /
			(theta->omega_larvae[g] / (theta->mu_L0*theta->d_E_larvae) - 1.0 / (theta->mu_L0*theta->d_L_larvae) - 1.0);		            	 // Larval carry capacity


		if (theta->Karry[g] < 1.0e-10) { theta->Karry[g] = 1.0e-10; } //
	}



	for (int g = 0; g < N_spec; g++)
	{
		if (g == 0) { cout << "An. farauti:  " << 100.0 * theta->Prop_mosq[0] << "%" << endl; }
		if (g == 1) { cout << "An. punctulatus:  " << 100.0 * theta->Prop_mosq[1] << "%" << endl; }
		if (g == 2) { cout << "An. koliensis:  " << 100.0 * theta->Prop_mosq[2] << "%" << endl; }

		cout << "EL_M  " << POP->yM[g][0] << endl;
		cout << "LL_M  " << POP->yM[g][1] << endl;
		cout << "P_M  " << POP->yM[g][2] << endl;
		cout << "S_M  " << POP->yM[g][3] << endl;
		cout << "E_M  " << POP->yM[g][4] << endl;
		cout << "I_M  " << POP->yM[g][5] << endl;

		cout << "lam_M = " << theta->lam_M[g] << endl;

		cout << "I_M = " << POP->yM[g][5] << endl;

		cout << "mm = " << theta->mm_0[g] << endl;

		cout << endl;
	}
	cout << endl;

	cout << "lam_H = " << theta->bb*theta->EIR_equil << endl;

	double EIR_out = 0.0;
	for (int g = 0; g < N_spec; g++)
	{
		EIR_out = EIR_out + 365.0*theta->aa[g] * POP->yM[g][5];
	}

	cout << "EIR = " << EIR_out << endl;


	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	// 3.7.4. Proportion in each age and heterogeneity stratified category

	//////////////////////////////////////////
	// Fill out vector of lagged lam_M*S_M

	theta->lam_S_M_track.resize(N_spec);

	for (int g = 0; g<N_spec; g++)
	{
		for (int k = 0; k < theta->M_track; k++)
		{
			theta->lam_S_M_track[g].push_back(theta->lam_M[g] * POP->yM[g][3]);
		}
	}


	//////////////////////////////////////////////////////////////
	// 3.7.3. Calculate probability for each compartment 

	vector<vector<vector<double>>> yH_eq_cumsum;
	yH_eq_cumsum.resize(N_age);
	for (int i = 0; i<N_age; i++)
	{
		yH_eq_cumsum[i].resize(N_het);
		for (int j = 0; j<N_het; j++)
		{
			yH_eq_cumsum[i][j].resize(N_H_comp*(K_max + 1));
		}
	}


	for (int i = 0; i<N_age; i++)
	{
		for (int j = 0; j<N_het; j++)
		{
			yH_eq_cumsum[i][j][0] = yH_eq[i][j][0][0];

			for (int k = 1; k<(K_max + 1); k++)
			{
				yH_eq_cumsum[i][j][k] = yH_eq_cumsum[i][j][k - 1] + yH_eq[i][j][k][0];
			}

			for (int c = 1; c < N_H_comp; c++)
			{
				yH_eq_cumsum[i][j][c*(K_max + 1)] = yH_eq_cumsum[i][j][c*(K_max + 1) - 1] + yH_eq[i][j][0][c];

				for (int k = 1; k<(K_max + 1); k++)
				{
					yH_eq_cumsum[i][j][c*(K_max + 1) + k] = yH_eq_cumsum[i][j][c*(K_max + 1) + k - 1] + yH_eq[i][j][k][c];
				}
			}
		}
	}


	for (int i = 0; i < N_age; i++)
	{
		for (int j = 0; j < N_het; j++)
		{
			for (int k = 0; k < N_H_comp*(K_max + 1); k++)
			{
				yH_eq_cumsum[i][j][k] = yH_eq_cumsum[i][j][k] / yH_eq_cumsum[i][j][N_H_comp*(K_max + 1) - 1];
			}
		}
	}


	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	// 3.7.4. Initialise Population of individuals                           

	///////////////////////////////////////////////////////////////////////////
	// 3.7.4.1. Temporary objects for setting up individuals            

	double rand_comp;
	double age_start, zeta_start;

	int i_index, j_index;


	float GMN_parm[(N_int)*(N_int + 3) / 2 + 1];
	float GMN_work[N_int];
	float GMN_zero[N_int];
	float zz_GMN[N_int];

	for (int k = 0; k<N_int; k++)
	{
		GMN_zero[k] = 0.0;
	}


	///////////////////////////////////////////////////////////////////////////
	// 3.7.4.2 Loop through and create N_pop individuals            

	for (int n = 0; n<POP->N_pop; n++)
	{
		//////////////////////////////////////////////////////////////////
		// 3.7.4.2.1. Assign age and heterogeneity 

		age_start = genexp(theta->age_mean);

		while (age_start > theta->age_max)
		{
			age_start = genexp(theta->age_mean);
		}

		zeta_start = exp(gennor(-0.5*theta->sig_het*theta->sig_het, theta->sig_het));

		while (zeta_start > theta->het_max)
		{
			zeta_start = exp(gennor(-0.5*theta->sig_het*theta->sig_het, theta->sig_het));
		}


		//////////////////////////////////////////////////////////////////
		// 3.7.4.2.2. Find appropriate age and het compartment 

		i_index = 0;

		for (int i = 0; i<N_age; i++)
		{
			if ((age_start > age_bounds[i]) && (age_start <= age_bounds[i + 1]))
			{
				i_index = i;
			}
		}


		j_index = 0;

		for (int j = 0; j<N_het; j++)
		{
			if ((zeta_start > POP->x_het_bounds[j]) && (zeta_start < POP->x_het_bounds[j + 1]))
			{
				j_index = j;
			}
		}

		//////////////////////////////////////////////////////////////////
		// 3.7.4.2.3. Construct a new individual

		double q_rand;

		individual HH(age_start, zeta_start);

		HH.S = 0;
		HH.I_PCR = 0;
		HH.I_LM = 0;
		HH.I_D = 0;
		HH.T = 0;
		HH.P = 0;


		if (genunf(0.0, 1.0) < 0.5)
		{
			HH.gender = 0;
		}
		else {
			HH.gender = 1;
		}

		if (HH.gender == 0)
		{
			if (genunf(0.0, 1.0) < theta->G6PD_prev)
			{
				HH.G6PD_def = 1;
			}
			else {
				HH.G6PD_def = 0;
			}
		}
		else {
			q_rand = genunf(0.0, 1.0);

			if (q_rand <= theta->G6PD_prev*theta->G6PD_prev)
			{
				HH.G6PD_def = 2;
			}

			if ((q_rand > theta->G6PD_prev*theta->G6PD_prev) && (q_rand <= theta->G6PD_prev*theta->G6PD_prev + 2 * theta->G6PD_prev*(1.0 - theta->G6PD_prev)))
			{
				HH.G6PD_def = 1;
			}

			if (q_rand >  theta->G6PD_prev*theta->G6PD_prev + 2 * theta->G6PD_prev*(1.0 - theta->G6PD_prev))
			{
				HH.G6PD_def = 0;
			}
		}


		if (genunf(0.0, 1.0) < theta->CYP2D6_prev)
		{
			HH.CYP2D6 = 1;
		}
		else {
			HH.CYP2D6 = 0;
		}


		HH.T_last_BS = 1000000.0;

		///////////////////////////////////////////
		//  3.7.4.2.4. An indicator for pregnancy appropriate age
		//             Only women between ages 18 and 40 can be pregnant.

		if (HH.gender == 1)
		{
			if ((HH.age > 6570.0) && (HH.age < 14600.0))
			{
				HH.preg_age = 1;
			}
			else {
				HH.preg_age = 0;
			}
		}

		HH.pregnant = 0;
		HH.preg_timer = 0.0;


		///////////////////////////////////////////////////////////////////
		// Randomly assign a state according to equilibrium probabilities

		rand_comp = genunf(0.0, 1.0);

		if (rand_comp <= yH_eq_cumsum[i_index][j_index][0])
		{
			HH.S = 1;
			HH.Hyp = 0;
		}

		for (int k = 1; k < (K_max + 1); k++)
		{
			if ((rand_comp > yH_eq_cumsum[i_index][j_index][k - 1]) && (rand_comp <= yH_eq_cumsum[i_index][j_index][k]))
			{
				HH.S = 1;
				HH.Hyp = k;
			}
		}

		for (int c = 1; c < N_H_comp; c++)
		{
			for (int k = 0; k < (K_max + 1); k++)
			{
				if ((rand_comp > yH_eq_cumsum[i_index][j_index][c*(K_max + 1) + k - 1]) && (rand_comp <= yH_eq_cumsum[i_index][j_index][c*(K_max + 1) + k]))
				{
					if (c == 1) { HH.I_PCR = 1; }
					if (c == 2) { HH.I_LM = 1; }
					if (c == 3) { HH.I_D = 1; }
					if (c == 4) { HH.T = 1; }
					if (c == 5) { HH.P = 1; }

					HH.Hyp = k;
				}
			}
		}

		if (rand_comp > yH_eq_cumsum[i_index][j_index][N_H_comp*(K_max + 1) - 1])
		{
			HH.P = 1;
			HH.Hyp = K_max;
		}


		////////////////////////////////////////////
		//  3.7.4.2.5. Initialise immunity

		HH.A_par_mat = A_par_eq_mean[POP->index_age_20][j_index] * theta->P_mat*exp(-POP->age_mids[i_index] / theta->d_mat);
		HH.A_clin_mat = A_clin_eq_mean[POP->index_age_20][j_index] * theta->P_mat*exp(-POP->age_mids[i_index] / theta->d_mat);

		HH.A_par = A_par_eq[i_index][j_index][HH.Hyp] - HH.A_par_mat;
		HH.A_clin = A_clin_eq[i_index][j_index][HH.Hyp] - HH.A_clin_mat;


		HH.A_par_boost = 1;
		HH.A_clin_boost = 1;

		HH.A_par_timer = -1.0;
		HH.A_clin_timer = -1.0;

		HH.PQ_proph = 0;
		HH.PQ_proph_timer = -1.0;


		//////////////////////////////////////////////////
		// TO DO: set this up in equilibrium

		if ((HH.age > 6570.0) && (HH.age < 14600.0))
		{
			if (genunf(0.0, 1.0) < 0.05)
			{
				HH.pregnant = 1;
				HH.preg_timer = 0;
			}
		}


		////////////////////////////////////////////
		//  3.7.4.2.6. A vector for storing lagged force of infection

		for (int k = 0; k<theta->H_track; k++)
		{
			HH.lam_bite_track.push_back(lam_eq[i_index][j_index]);
		}

		for (int k = 0; k<theta->H_track; k++)
		{
			HH.lam_rel_track.push_back(HH.Hyp*theta->ff);
		}


		////////////////////////////////////////////////////////
		// 3.7.4.2.7. Give individuals their life-long intervention access score

		for (int p = 0; p<N_int; p++)
		{
			for (int q = 0; q<N_int; q++)
			{
				theta->V_int_dummy[p][q] = theta->V_int[p][q];
			}
		}

		setgmn(GMN_zero, *theta->V_int_dummy, N_int, GMN_parm);

		genmn(GMN_parm, zz_GMN, GMN_work);

		for (int k = 0; k<N_int; k++)
		{
			HH.zz_int[k] = zz_GMN[k];
		}


		////////////////////////////////////////////////////////
		// 3.7.4.2.8. Individuals begin without interventions

		HH.LLIN = 0;
		HH.IRS = 0;

		for (int g = 0; g < N_spec; g++)
		{
			HH.w_VC[g] = 1.0;
			HH.y_VC[g] = 1.0;
			HH.z_VC[g] = 0.0;
		}


		///////////////////////////////////////////////////
		////////////////////////////////////////////////////
		// 3.7.4.3.. Add this person to the vector of people

		POP->people.push_back(HH);
	}

	///////////////////////////////////////////////////////////
	// 3.7.5.1. Proportion of bites received by each person

	POP->pi_n.resize(POP->N_pop);
	for (int n = 0; n < POP->N_pop; n++)
	{
		POP->pi_n[n].resize(N_spec);
	}

	for (int n = 0; n<POP->N_pop; n++)
	{
		for (int g = 0; g < N_spec; g++)
		{
			POP->pi_n[n][g] = POP->people[n].zeta_het*(1.0 - theta->rho_age*exp(-POP->people[n].age*theta->age_0_inv));
		}
	}



	double SIGMA_PI[N_spec];

	for (int g = 0; g < N_spec; g++)
	{
		SIGMA_PI[g] = 0.0;
		for (int n = 0; n<POP->N_pop; n++)
		{
			SIGMA_PI[g] = SIGMA_PI[g] + POP->pi_n[n][g];
		}

		for (int n = 0; n<POP->N_pop; n++)
		{
			POP->pi_n[n][g] = POP->pi_n[n][g] / SIGMA_PI[g];
		}
	}

	///////////////////////////////////////////////////////////
	// 3.7.5. Initialise population level quantitites

	//////////////////////////////////////////////
	// 3.7.5.1. Vector control quantities

	for (int g = 0; g < N_spec; g++)
	{
		POP->SUM_pi_w[g] = 0;
		POP->SUM_pi_z[g] = 0;


		for (int n = 0; n < POP->N_pop; n++)
		{
			POP->SUM_pi_w[g] = POP->SUM_pi_w[g] + POP->pi_n[n][g] * POP->people[n].w_VC[g];
			POP->SUM_pi_z[g] = POP->SUM_pi_z[g] + POP->pi_n[n][g] * POP->people[n].z_VC[g];
		}
	}


	for (int g = 0; g < N_spec; g++)
	{
		POP->W_VC[g] = 1.0 - theta->Q_0[g] + theta->Q_0[g] * POP->SUM_pi_w[g];
		POP->Z_VC[g] = theta->Q_0[g] * POP->SUM_pi_z[g];

		POP->delta_1_VC[g] = theta->delta_1 / (1.0 - POP->Z_VC[g]);
		POP->delta_VC[g] = POP->delta_1_VC[g] + theta->delta_2;

		POP->p_1_VC[g] = exp(-theta->mu_M[g] * POP->delta_1_VC[g]);
		POP->mu_M_VC[g] = -log(POP->p_1_VC[g] * theta->p_2[g]) / POP->delta_VC[g];

		POP->Q_VC[g] = 1.0 - (1.0 - theta->Q_0[g]) / POP->W_VC[g];

		POP->aa_VC[g] = POP->Q_VC[g] / POP->delta_VC[g];

		POP->exp_muM_tauM_VC[g] = exp(-POP->mu_M_VC[g] * theta->tau_M[g]);
		POP->beta_VC[g] = theta->eps_max[g] * POP->mu_M_VC[g] / (exp(POP->delta_VC[g] * POP->mu_M_VC[g]) - 1.0);
	}


	/////////////////////////////////////////////////////////////////////////
	// 3.7.5.2. The rate at which person n is bitten by a single mosquito

	POP->lam_n.resize(POP->N_pop);
	for (int n = 0; n < POP->N_pop; n++)
	{
		POP->lam_n[n].resize(N_spec);
	}


	for (int n = 0; n < POP->N_pop; n++)
	{
		for (int g = 0; g < N_spec; g++)
		{
			POP->lam_n[n][g] = theta->aa[g] * POP->pi_n[n][g];
		}
	}


	/////////////////////////////////////////////////////////////////////////
	// 3.7.5.3. Output an overview of the initial set up

	cout << "Equilibrium set up......." << endl;

	double S_ind = 0.0, I_PCR_ind = 0.0, I_LM_ind = 0.0, I_D_ind = 0.0, T_ind = 0.0, P_ind = 0.0;
	double S_eqq = 0.0, I_PCR_eqq = 0.0, I_LM_eqq = 0.0, I_D_eqq = 0.0, T_eqq = 0.0, P_eqq = 0.0;

	for (int n = 0; n<POP->N_pop; n++)
	{
		S_ind = S_ind + POP->people[n].S;
		I_PCR_ind = I_PCR_ind + POP->people[n].I_PCR;
		I_LM_ind = I_LM_ind + POP->people[n].I_LM;
		I_D_ind = I_D_ind + POP->people[n].I_D;
		T_ind = T_ind + POP->people[n].T;
		P_ind = P_ind + POP->people[n].P;
	}


	for (int i = 0; i<N_age; i++)
	{
		for (int j = 0; j<N_het; j++)
		{
			for (int k = 0; k < (K_max + 1); k++)
			{
				S_eqq = S_eqq + yH_eq[i][j][k][0];
				I_PCR_eqq = I_PCR_eqq + yH_eq[i][j][k][1];
				I_LM_eqq = I_LM_eqq + yH_eq[i][j][k][2];
				I_D_eqq = I_D_eqq + yH_eq[i][j][k][3];
				T_eqq = T_eqq + yH_eq[i][j][k][4];
				P_eqq = P_eqq + yH_eq[i][j][k][5];
			}
		}
	}


	cout << "S = " << ((double)S_ind) / POP->N_pop << "\t" << S_eqq << endl;
	cout << "I_PCR = " << ((double)I_PCR_ind) / POP->N_pop << "\t" << I_PCR_eqq << endl;
	cout << "I_LM = " << ((double)I_LM_ind) / POP->N_pop << "\t" << I_LM_eqq << endl;
	cout << "I_D = " << ((double)I_D_ind) / POP->N_pop << "\t" << I_D_eqq << endl;
	cout << "T = " << ((double)T_ind) / POP->N_pop << "\t" << T_eqq << endl;
	cout << "P = " << ((double)P_ind) / POP->N_pop << "\t" << P_eqq << endl;

}


////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////// 
//          //                                                                    // 
//  ##      //  #### #    ## ####   #### ##   ## #### ####   ##  ##  ####  ##     //
//  ## ##   //   ##  ##   ## ## ##   ##  ##   ##  ##  ## ##  ##  ## ##  ## ##     //
//  ######  //   ##  ###  ## ##  ##  ##   ## ##   ##  ##  ## ##  ## ###### ##     //
//     ##   //   ##  ## #### ## ##   ##    ###    ##  ## ##  ##  ## ##  ## ##     //
//     ##   //  #### ##   ## ####   ####    #    #### ####    ####  ##  ## #####  //
//          //                                                                    //
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
//       //                                                                               //
// 4.1.  //  THE MODEL (within humans at least!)                                          //
//       //  Stochastic moves between compartments for each individual                    //
///////////  for a fixed time step                                                        //
///////////                                                                               //
///////////  TO DO: (i) lots of small stuff to speed things up                            //
///////////         (ii) might be able to save on some multiplications by not normalising //
///////////             vector of probabilities and doing it instead inside CH_sample     //
///////////         (iii) add new state for MDA prophylaxis                               //
///////////                                                                               //
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void individual::state_mover(params theta, double lam_bite)
{
	lam_bite_track.push_back(lam_bite);
	lam_bite_track.erase(lam_bite_track.begin());

	lam_bite_lag = lam_bite_track[0];


	lam_rel_track.push_back(((double)Hyp)*theta.ff);
	lam_rel_track.erase(lam_rel_track.begin());

	lam_rel_lag = lam_rel_track[0];


	lam_H_lag = lam_bite_lag + lam_rel_lag;


	///////////////////////////////////
	// Indicators for new events

	I_PCR_new = 0;
	I_LM_new = 0;
	I_D_new = 0;
	ACT_treat = 0;
	PQ_treat = 0;

	PQ_effective    = 0;
	PQ_overtreat    = 0;
	PQ_overtreat_9m = 0;


	//////////////////////////////////////////////////////////////////////
	//     //                                                           // 
	//  0  //  S: Susceptible                                           //
	//     //                                                           //
	//////////////////////////////////////////////////////////////////////

	if (S == 1)
	{
		theta.S_out = lam_H_lag;

		if (exp(-t_step*theta.S_out) < genunf(0, 1))
		{
			//theta.r_PCR   = 1.0/( theta.d_PCR_min + (theta.d_PCR_max-theta.d_PCR_min)/( 1.0 + pow((A_par+A_par_mat)*theta.A_PCR_50pc_inv,theta.K_PCR) )); 
			theta.phi_LM = theta.phi_LM_min + (theta.phi_LM_max - theta.phi_LM_min) / (1.0 + pow((A_par + A_par_mat)*theta.A_LM_50pc_inv, theta.K_LM));
			theta.phi_D = theta.phi_D_min + (theta.phi_D_max - theta.phi_D_min) / (1.0 + pow((A_clin + A_clin_mat)*theta.A_D_50pc_inv, theta.K_D));


			theta.S_move[0] = (1.0 - theta.phi_LM);                                  // Move to I_PCR  //  lam_H_lag*(1.0-theta.phi_LM)/theta.S_out; 
			theta.S_move[1] = theta.phi_LM*(1.0 - theta.phi_D);                      // Move to I_LM   //  lam_H_lag*theta.phi_LM*(1.0-theta.phi_D)/theta.S_out;
			theta.S_move[2] = theta.phi_LM*theta.phi_D*(1.0 - theta.treat_BScover);  // Move to I_D    //  lam_H_lag*theta.phi_LM*(1.0-theta.phi_D)*(1.0-theta.treat_cov)/theta.S_out;
			theta.S_move[3] = theta.phi_LM*theta.phi_D*theta.treat_BScover;          // Move to T      //  lam_H_lag*theta.phi_LM*(1.0-theta.phi_D)*theta.treat_cov/theta.S_out;

			CH_move = CH_sample(theta.S_move, 4);

			////////////////////////////////
			// S -> I_PCR

			if (CH_move == 0)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (PQ_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				I_PCR_new = 1;

				T_last_BS = 0.0;

				S = 0;
				I_PCR = 1;

				return;
			}

			////////////////////////////////
			// S -> I_LM

			if (CH_move == 1)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (PQ_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				I_LM_new = 1;
				I_PCR_new = 1;

				T_last_BS = 0.0;

				S = 0;
				I_LM = 1;

				return;
			}

			////////////////////////////////
			// S -> I_D

			if (CH_move == 2)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (PQ_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				I_PCR_new = 1;
				I_LM_new = 1;
				I_D_new = 1;

				T_last_BS = 0.0;

				S = 0;
				I_D = 1;

				return;
			}

			////////////////////////////////
			// S -> T

			if (CH_move == 3)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (PQ_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				/////////////////////////////////////////////////////////////////////
				// Blood-stage drug administered

				ACT_treat = 1;

				if( genunf(0.0, 1.0) < theta.treat_BSeff )
				{
					S = 0;
					T = 1;
				}


				I_PCR_new = 1;
				I_LM_new = 1;
				I_D_new = 1;

				T_last_BS = 0.0;


				/////////////////////////////////////////////////////////////////////
				// Is PQ administered?
				// For efficiency, the PQ related commands are only implemented if
				// PQ is available

				if( genunf(0.0, 1.0) < theta.PQ_treat_PQavail )
				{
					PQ_treat = 1;


					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of G6PD deficiency

					if ((theta.PQ_treat_G6PD_risk == 1) && (G6PD_def == 1))
					{
						PQ_treat = 0;
					}

					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of pregancy

					if ((theta.PQ_treat_preg_risk == 1) && (pregnant == 1))
					{
						PQ_treat = 0;
					}

					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of young age

					if (age < theta.PQ_treat_low_age)
					{
						PQ_treat = 0;
					}


					/////////////////////////////////////////////////////////////////////
					// Is PQ effective?

					PQ_effective = 0;

					if (genunf(0.0, 1.0) < theta.PQ_treat_PQeff)
					{
						PQ_effective = 0;
					}

					if ((theta.PQ_treat_CYP2D6_risk == 1) && (CYP2D6 == 1))
					{
						PQ_effective = 0;
					}

					if ((PQ_treat == 1) && (PQ_effective == 1))
					{
						Hyp = 0;
						PQ_proph = 1;
					}
				}


				return;
			}
		}
	}


	//////////////////////////////////////////////////////////////////////
	//     //                                                           // 
	//  1  //  I_PCR: PCR detectable BS infections                      //
	//     //                                                           //
	//////////////////////////////////////////////////////////////////////

	if (I_PCR == 1)
	{
		theta.r_PCR = 1.0 / (theta.d_PCR_min + (theta.d_PCR_max - theta.d_PCR_min) / (1.0 + pow((A_par + A_par_mat)*theta.A_PCR_50pc_inv, theta.K_PCR)));

		theta.I_PCR_out = lam_H_lag + theta.r_PCR;

		if (exp(-t_step*theta.I_PCR_out) < genunf(0, 1))
		{
			theta.phi_LM = theta.phi_LM_min + (theta.phi_LM_max - theta.phi_LM_min) / (1.0 + pow((A_par + A_par_mat)*theta.A_LM_50pc_inv, theta.K_LM));
			theta.phi_D = theta.phi_D_min + (theta.phi_D_max - theta.phi_D_min) / (1.0 + pow((A_clin + A_clin_mat)*theta.A_D_50pc_inv, theta.K_D));


			theta.I_PCR_move[0] = theta.r_PCR / theta.I_PCR_out;                                                       // Move to S 
			theta.I_PCR_move[1] = lam_H_lag*(1 - theta.phi_LM) / theta.I_PCR_out;                                      // Move to I_PCR
			theta.I_PCR_move[2] = lam_H_lag*theta.phi_LM*(1.0 - theta.phi_D) / theta.I_PCR_out;                        // Move to I_LM
			theta.I_PCR_move[3] = lam_H_lag*theta.phi_LM*theta.phi_D*(1.0 - theta.treat_BScover) / theta.I_PCR_out;        // Move to I_D
			theta.I_PCR_move[4] = lam_H_lag*theta.phi_LM*theta.phi_D*theta.treat_BScover / theta.I_PCR_out;                // Move to T

			CH_move = CH_sample(theta.I_PCR_move, 5);

			////////////////////////////////
			// I_PCR -> S

			if (CH_move == 0)
			{
				I_PCR = 0;
				S = 1;

				return;
			}

			////////////////////////////////
			// I_PCR -> I_PCR
			// 
			// This is a super-infection event and we assumes boosting of immunity

			if (CH_move == 1)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (PQ_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				return;
			}

			////////////////////////////////
			// I_PCR -> I_LM

			if (CH_move == 2)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (PQ_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				I_PCR_new = 1;
				I_LM_new = 1;

				I_PCR = 0;
				I_LM = 1;

				return;
			}

			////////////////////////////////
			// I_PCR -> I_D

			if (CH_move == 3)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (PQ_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				I_PCR_new = 1;
				I_LM_new = 1;
				I_D_new = 1;

				I_PCR = 0;
				I_D = 1;

				return;
			}


			////////////////////////////////
			// I_PCR -> T

			if (CH_move == 4)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (PQ_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				/////////////////////////////////////////////////////////////////////
				// Blood-stage drug administered

				ACT_treat = 1;

				if (genunf(0.0, 1.0) < theta.treat_BSeff)
				{
					I_PCR = 0;
					T = 1;
				}


				I_PCR_new = 1;
				I_LM_new = 1;
				I_D_new = 1;

				T_last_BS = 0.0;


				/////////////////////////////////////////////////////////////////////
				// Is PQ administered?
				// For efficiency, the PQ related commands are only implemented if
				// PQ is available

				if (genunf(0.0, 1.0) < theta.PQ_treat_PQavail)
				{
					PQ_treat = 1;


					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of G6PD deficiency

					if ((theta.PQ_treat_G6PD_risk == 1) && (G6PD_def == 1))
					{
						PQ_treat = 0;
					}

					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of pregancy

					if ((theta.PQ_treat_preg_risk == 1) && (pregnant == 1))
					{
						PQ_treat = 0;
					}

					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of young age

					if (age < theta.PQ_treat_low_age)
					{
						PQ_treat = 0;
					}


					/////////////////////////////////////////////////////////////////////
					// Is PQ effective?

					PQ_effective = 0;

					if (genunf(0.0, 1.0) < theta.PQ_treat_PQeff)
					{
						PQ_effective = 1;
					}

					if ((theta.PQ_treat_CYP2D6_risk == 1) && (CYP2D6 == 1))
					{
						PQ_effective = 0;
					}

					if ((PQ_treat == 1) && (PQ_effective == 1))
					{
						Hyp = 0;
						PQ_proph = 1;
					}
				}


				return;
			}

		}
	}


	//////////////////////////////////////////////////////////////////////
	//     //                                                           // 
	//  2  //  I_LM: LM detectable BS infection                         //
	//     //                                                           //
	//////////////////////////////////////////////////////////////////////

	if (I_LM == 1)
	{
		theta.I_LM_out = lam_H_lag + theta.r_LM;

		if (exp(-t_step*theta.I_LM_out) < genunf(0, 1))
		{
			//theta.r_PCR = 1.0/( theta.d_PCR_min + (theta.d_PCR_max-theta.d_PCR_min)/( 1.0 + pow((A_par+A_par_mat)*theta.A_PCR_50pc_inv,theta.K_PCR) ) ); 
			theta.phi_LM = theta.phi_LM_min + (theta.phi_LM_max - theta.phi_LM_min) / (1.0 + pow((A_par + A_par_mat)*theta.A_LM_50pc_inv, theta.K_LM));
			theta.phi_D = theta.phi_D_min + (theta.phi_D_max - theta.phi_D_min) / (1.0 + pow((A_clin + A_clin_mat)*theta.A_D_50pc_inv, theta.K_D));


			theta.I_LM_move[0] = theta.r_LM / theta.I_LM_out;                                            // Move to I_PCR
			theta.I_LM_move[1] = lam_H_lag*(1.0 - theta.phi_D) / theta.I_LM_out;                         // Move to I_LM
			theta.I_LM_move[2] = lam_H_lag*theta.phi_D*(1.0 - theta.treat_BScover) / theta.I_LM_out;     // Move to I_D
			theta.I_LM_move[3] = lam_H_lag*theta.phi_D*theta.treat_BScover / theta.I_LM_out;             // Move to T

			CH_move = CH_sample(theta.I_LM_move, 4);

			////////////////////////////////
			// I_LM -> I_PCR

			if (CH_move == 0)
			{
				I_LM = 0;
				I_PCR = 1;

				return;
			}


			////////////////////////////////
			// I_LM -> I_LM
			//
			// Super-infection boosts immunity

			if (CH_move == 1)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (PQ_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				return;
			}


			////////////////////////////////
			// I_LM -> I_D

			if (CH_move == 2)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (PQ_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				I_PCR_new = 1;
				I_LM_new = 1;
				I_D_new = 1;

				I_LM = 0;
				I_D = 1;

				return;
			}


			////////////////////////////////
			// I_LM -> T

			if (CH_move == 3)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (PQ_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				/////////////////////////////////////////////////////////////////////
				// Blood-stage drug administered

				ACT_treat = 1;

				if (genunf(0.0, 1.0) < theta.treat_BSeff)
				{
					I_LM = 0;
					T = 1;
				}



				I_PCR_new = 1;
				I_LM_new = 1;
				I_D_new = 1;

				T_last_BS = 0.0;


				/////////////////////////////////////////////////////////////////////
				// Is PQ administered?
				// For efficiency, the PQ related commands are only implemented if
				// PQ is available

				if (genunf(0.0, 1.0) < theta.PQ_treat_PQavail)
				{
					PQ_treat = 1;


					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of G6PD deficiency

					if ((theta.PQ_treat_G6PD_risk == 1) && (G6PD_def == 1))
					{
						PQ_treat = 0;
					}

					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of pregancy

					if ((theta.PQ_treat_preg_risk == 1) && (pregnant == 1))
					{
						PQ_treat = 0;
					}

					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of young age

					if (age < theta.PQ_treat_low_age)
					{
						PQ_treat = 0;
					}


					/////////////////////////////////////////////////////////////////////
					// Is PQ effective?

					PQ_effective = 0;

					if (genunf(0.0, 1.0) < theta.PQ_treat_PQeff)
					{
						PQ_effective = 1;
					}

					if ((theta.PQ_treat_CYP2D6_risk == 1) && (CYP2D6 == 1))
					{
						PQ_effective = 0;
					}

					if ((PQ_treat == 1) && (PQ_effective == 1))
					{
						Hyp = 0;
						PQ_proph = 1;
					}
				}


				return;
			}

		}
	}


	//////////////////////////////////////////////////////////////////////
	//     //                                                           // 
	//  3  //  I_D: Clinical disease                                    //
	//     //                                                           //
	//////////////////////////////////////////////////////////////////////

	if (I_D == 1)
	{
		theta.I_D_out = lam_H_lag + theta.r_D;

		if (exp(-t_step*theta.I_D_out) < genunf(0, 1))
		{
			theta.I_D_move[0] = theta.r_D / theta.I_D_out;              // Move to I_LM
			theta.I_D_move[1] = lam_H_lag / theta.I_D_out;              // Move to D

			CH_move = CH_sample(theta.I_D_move, 2);

			////////////////////////////////
			// I_D -> I_LM

			if (CH_move == 0)
			{
				I_D = 0;
				I_LM = 1;

				return;
			}


			////////////////////////////////
			// I_D -> I_D

			if (CH_move == 1)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (PQ_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				return;
			}

		}
	}


	//////////////////////////////////////////////////////////////////////
	//     //                                                           // 
	//  4  //  T : Clinical episode under treatment                     //
	//     //                                                           //
	//////////////////////////////////////////////////////////////////////

	if (T == 1)
	{
		theta.T_out = lam_H_lag + theta.r_T;

		if (exp(-t_step*theta.T_out) < genunf(0, 1))
		{
			theta.T_move[0] = theta.r_T / theta.T_out;              // Move to P
			theta.T_move[1] = lam_H_lag / theta.T_out;              // Move to T

			CH_move = CH_sample(theta.T_move, 2);

			////////////////////////////////
			// T -> P

			if (CH_move == 0)
			{
				T = 0;
				P = 1;

				return;
			}


			////////////////////////////////
			// T -> T
			//
			// Note that we don't allow immune boosting but we do
			// allow for the accumulation of new hypnozoites


			if (CH_move == 1)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (PQ_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				return;
			}

		}
	}




	//////////////////////////////////////////////////////////////////////
	//     //                                                           // 
	//  5  //  P: Treatment prophylaxis                                 //
	//     //                                                           //
	//////////////////////////////////////////////////////////////////////

	if (P == 1)
	{
		theta.P_out = lam_H_lag + theta.r_P;

		if (exp(-t_step*theta.P_out) < genunf(0, 1))
		{
			theta.P_move[0] = theta.r_P / theta.P_out;              // Move to S
			theta.P_move[1] = lam_H_lag / theta.P_out;              // Move to P

			CH_move = CH_sample(theta.P_move, 2);

			////////////////////////////////
			// P -> S

			if (CH_move == 0)
			{
				P = 0;
				S = 1;

				return;
			}


			////////////////////////////////
			// P -> P
			//
			// Note that we don't allow immune boosting but we do
			// allow for the accumulation of new hypnozoites


			if (CH_move == 1)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (PQ_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				return;
			}

		}
	}

}



////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
//       //                                                                               //
// 4.2.  //  Ageing and immune boosting                                                   //
//       //                                                                               //  
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////


void individual::ager(params theta)
{
	/////////////////////////
	// Ageing

	age = age + t_step;


	/////////////////////////
	// Loss of hypnozoites

	if (Hyp > K_max)
	{
		Hyp = K_max;
	}

	if (1.0 - exp(-t_step*theta.gamma_L*Hyp) > genunf(0, 1))
	{
		Hyp = Hyp - 1;
	}


	/////////////////////////
	// Immune decay

	A_par = A_par*theta.A_par_decay;
	A_clin = A_clin*theta.A_clin_decay;


	///////////////////////////////////////////////////
	// Maternal immunity decays exponentially for the first year
	// and is then set to zero afterwards

	if (age < 365.0)
	{
		A_par_mat = A_par_mat*theta.mat_decay;
		A_clin_mat = A_clin_mat*theta.mat_decay;
	}
	else {
		A_par_mat = 0.0;
		A_clin_mat = 0.0;
	}


	/////////////////////////
	// Of child-bearing age? ~ 20 years
	// age in (18, 22) years

	if (gender == 1)
	{
		if ((age>6570.0) && (age<14600.0))
		{
			preg_age = 1;
		}
		else {
			preg_age = 0;
		}

		if (pregnant == 1)
		{
			preg_timer = preg_timer + 1.0;
		}

		if (pregnant == 0)
		{
			if (preg_age == 1)
			{
				if (genunf(0.0, 1.0) < theta.P_preg)
				{
					pregnant = 1;
					preg_timer = 0.0;
				}
			}
		}

		if (pregnant == 1)
		{
			if (preg_timer > 270.0)
			{
				pregnant = 0;
				preg_timer = 0.0;
			}
		}
	}


	//////////////////////////////////////////////////////
	// Switches for refractory period of immune boosting

	if (A_par_boost == 0)
	{
		A_par_timer = A_par_timer - t_step;

		if (A_par_timer < 0.0)
		{
			A_par_boost = 1;
		}
	}

	if (A_clin_boost == 0)
	{
		A_clin_timer = A_clin_timer - t_step;

		if (A_clin_timer < 0.0)
		{
			A_clin_boost = 1;
		}
	}


	//////////////////////////////////////////////////////
	// Switches for primaquine prophylaxis

	if (PQ_proph == 1)
	{
		PQ_proph_timer = PQ_proph_timer - t_step;

		if (PQ_proph_timer < 0.0)
		{
			PQ_proph = 0;
		}
	}


	//////////////////////////////////////////////////////
	// Time since last blood-stage infection

	if (S == 1)
	{
		T_last_BS = T_last_BS + 1.0;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
//       //                                                                               //
// 4.3.  //  Individual-level vector control updater                                      //
//       //                                                                               //
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void individual::intervention_updater(params theta)
{

	///////////////////////////////////////////
	// Is LLIN lost

	if (LLIN == 1)
	{
		if (theta.P_LLIN_loss > genunf(0, 1))
		{
			LLIN = 0;
		}
	}

	///////////////////////////////////////////
	// IRS turned off if sufficiently long ago (3 years???)

	if (IRS_age > 1095.0)
	{
		IRS = 0;
	}


	if (LLIN == 1 && IRS == 0)
	{
		LLIN_age = LLIN_age + t_step;

		for (int g = 0; g < N_spec; g++)
		{
			r_LLIN[g] = theta.r_LLIN_net[g] + (r_LLIN[g] - theta.r_LLIN_net[g])*theta.P_PYR_decay;
			d_LLIN[g] = d_LLIN[g] * theta.P_PYR_decay;
			s_LLIN[g] = 1.0 - r_LLIN[g] - d_LLIN[g];

			w_VC[g] = 1 - theta.PSI_bed[g] + theta.PSI_bed[g] * s_LLIN[g];
			y_VC[g] = w_VC[g];
			z_VC[g] = theta.PSI_bed[g] * r_LLIN[g];
		}
	}

	if (LLIN == 0 && IRS == 1)
	{
		IRS_age = IRS_age + t_step;

		for (int g = 0; g < N_spec; g++)
		{
			r_IRS[g] = r_IRS[g] * theta.P_IRS_decay;
			d_IRS[g] = d_IRS[g] * theta.P_IRS_decay;
			s_IRS[g] = 1.0 - r_IRS[g] - d_IRS[g];

			w_VC[g] = 1 - theta.PSI_indoors[g] + theta.PSI_indoors[g] * (1.0 - r_IRS[g])*s_IRS[g];
			y_VC[g] = 1 - theta.PSI_indoors[g] + theta.PSI_indoors[g] * (1.0 - r_IRS[g]);
			z_VC[g] = theta.PSI_indoors[g] * r_IRS[g];
		}
	}

	if (LLIN == 1 && IRS == 1)
	{
		LLIN_age = LLIN_age + t_step;
		IRS_age = IRS_age + t_step;

		for (int g = 0; g < N_spec; g++)
		{
			r_LLIN[g] = theta.r_LLIN_net[g] + (r_LLIN[g] - theta.r_LLIN_net[g])*theta.P_PYR_decay;
			d_LLIN[g] = d_LLIN[g] * theta.P_PYR_decay;
			s_LLIN[g] = 1.0 - r_LLIN[g] - d_LLIN[g];

			r_IRS[g] = r_IRS[g] * theta.P_IRS_decay;
			d_IRS[g] = d_IRS[g] * theta.P_IRS_decay;
			s_IRS[g] = 1.0 - r_IRS[g] - d_IRS[g];

			w_VC[g] = 1.0 - theta.PSI_indoors[g] + theta.PSI_bed[g] * (1.0 - r_IRS[g])*s_LLIN[g] * s_IRS[g] + (theta.PSI_indoors[g] - theta.PSI_bed[g])*(1.0 - r_IRS[g])*s_IRS[g];
			y_VC[g] = 1.0 - theta.PSI_indoors[g] + theta.PSI_bed[g] * (1.0 - r_IRS[g])*s_LLIN[g] + (theta.PSI_indoors[g] - theta.PSI_bed[g])*(1.0 - r_IRS[g]);
			z_VC[g] = theta.PSI_bed[g] * (1.0 - r_IRS[g])*r_LLIN[g] + theta.PSI_indoors[g] * r_IRS[g];
		}

	}

	if (LLIN == 0 && IRS == 0)
	{
		for (int g = 0; g < N_spec; g++)
		{
			w_VC[g] = 1.0;
			y_VC[g] = 1.0;
			z_VC[g] = 0.0;
		}

	}

}

