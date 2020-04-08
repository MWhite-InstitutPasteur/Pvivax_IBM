/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///                                                                       ///
///  Individual-based Plasmodium vivax transmission model.                ///
///                                                                       ///
///  Dr Michael White                                                     ///
///  Imperial College London                                              ///
///  m.white08@imperial.ac.uk                                             ///
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
#include "rng.h"
#include "model.h"
#include <omp.h>
#include <vector>
#include <algorithm>
#include <Rcpp.h>

using namespace std;

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


//' Entrypoint for the simulation
//' 
//' @param model_param_path, the path of the model parameter file
//' @param fara_param_path, the path to the farauti mosquito parameter files.
//' @param punc_param_path, the path to the punctulatus mosquito parameter files.
//' @param koli_param_path, the path to the koliensis mosquito parameter files.
//' @param coverage_param_path, the path of the coverage parameter file
//' @param output_path, the path to write the results to
//' @export
// [[Rcpp::export]]
int run_simulation_from_path(
    std::string model_param_path,
    std::string fara_param_path,
    std::string punc_param_path,
    std::string koli_param_path,
    std::string coverage_param_path,
    std::string output_path
) {

  const char* mosquito_File[N_spec_max] = {
    fara_param_path.c_str(),
    punc_param_path.c_str(),
    koli_param_path.c_str()
  };

  return run_simulation(
    model_param_path.c_str(),
    mosquito_File,
    coverage_param_path.c_str(),
    output_path.c_str()
  );
}

int run_simulation(
  const char* parameter_File,
  const char** mosquito_File,
  const char* coverage_File,
  const char* output_File
) {

	clock_t clock_time;
	clock_time = clock();


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

	Rcpp::Rcout << "Reading in parameter file............." << endl;
	Rcpp::Rcout << endl;

	string discard;

	std::ifstream parameter_Stream(parameter_File);

	if (parameter_Stream.fail())
	{
		Rcpp::Rcout << "Failure reading in data." << endl;
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

	parameter_Stream >> discard >> Pv_mod_par.BS_treat_cov_base >> discard;      // proportion of episodes of symptomatic disease treated (baseline)
	parameter_Stream >> discard >> Pv_mod_par.BS_treat_eff_base >> discard;      // efficacy of front-line treatment (baseline)
	parameter_Stream >> discard >> Pv_mod_par.BS_treat_BSproph_base >> discard;  // duration of prophylaxis of front-line treatment (baseline)

	parameter_Stream >> discard >> Pv_mod_par.A_PCR_50pc >> discard;           // PCR_detectable infection scale parameter
	parameter_Stream >> discard >> Pv_mod_par.K_PCR >> discard;                // PCR_detectable infection shape parameter


	Pv_mod_par.H_track = int(Pv_mod_par.d_latent / t_step);                                 // Number of time steps for duration of latency

	Pv_mod_par.treat_cov = Pv_mod_par.BS_treat_cov_base;                       // Treatment coverage
	Pv_mod_par.treat_eff = Pv_mod_par.BS_treat_eff_base;                       // Efficacy of treatment
	Pv_mod_par.r_P = 1.0 / Pv_mod_par.BS_treat_BSproph_base;                   // rate of recovery from prophylaxis

	Pv_mod_par.PQ_treat_cover   = 0.0;
	Pv_mod_par.PQ_treat_PQcover = 0.0;
	Pv_mod_par.PQ_treat_BSeff   = 0.0;
	Pv_mod_par.PQ_treat_PQeff   = 0.0;
	Pv_mod_par.PQ_treat_BSproph = 10.0;

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


	Pv_mod_par.PQ_treat_CYP2D6 = 1;


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

	Rcpp::Rcout << "Parameter values read in from file!" << endl;
	Rcpp::Rcout << endl;


	////////////////////////////////////////////
	//                                        //
	// 1.4. Read in mosquito parameters       //
	//                                        //
	////////////////////////////////////////////

	Rcpp::Rcout << "Reading in mosquito files............." << endl;
	Rcpp::Rcout << endl;

	for (int g = 0; g < N_spec; g++)
	{

		std::ifstream mosquito_Stream(mosquito_File[g]);

		if (mosquito_Stream.fail())
		{
			Rcpp::Rcout << "Failure reading in mosquito parameters." << endl;
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
		Pv_mod_par.omega_larvae[g] = -0.5*Pv_mod_par.omega_larvae[g] + sqrt(0.25*Pv_mod_par.omega_larvae[g] * Pv_mod_par.omega_larvae[g] + 0.5*Pv_mod_par.gamma_larvae*Pv_mod_par.beta_larvae*Pv_mod_par.mu_L0*Pv_mod_par.d_E_larvae /
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

	Rcpp::Rcout << "Mosquito parameter values read in from file!" << endl;
	Rcpp::Rcout << endl;


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

	Pv_mod_par.A_par_decay = exp(-Pv_mod_par.r_par*t_step);
	Pv_mod_par.A_clin_decay = exp(-Pv_mod_par.r_clin*t_step);
	Pv_mod_par.mat_decay = exp(-Pv_mod_par.d_mat*t_step);

	Pv_mod_par.age_0_inv = 1.0 / Pv_mod_par.age_0;                 // Inverse of age-dependent biting parameter

	Pv_mod_par.A_PCR_50pc_inv = log2 / Pv_mod_par.A_PCR_50pc;      // Immune scalar for clearance of infection
	Pv_mod_par.A_LM_50pc_inv = 1.0 / Pv_mod_par.A_LM_50pc;        // Immune scalar for BS infection
	Pv_mod_par.A_D_50pc_inv = 1.0 / Pv_mod_par.A_D_50pc;         // Immune scalar for clinical disease

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
	//        Note that the matrix we read in may have variable size

	int N_cov_rounds = 0;

	Rcpp::Rcout << "Read in intervention coverage file............" << endl;

	std::ifstream coverage_Stream(coverage_File);

	if (coverage_Stream.fail())
	{
		Rcpp::Rcout << "Failure reading in data." << endl;
	}


	vector<vector<double>> coverage;
	coverage.resize(0);

	vector<double> cov_read(24);

	do {
		for (int j = 0; j<24; j++)
		{
			coverage_Stream >> cov_read[j];
		}
		if (cov_read[0] > -0.5)
		{
			N_cov_rounds = N_cov_rounds + 1;
			coverage.push_back(cov_read);
		}
	} while (cov_read[0] > -0.5);

	Rcpp::Rcout << "Intervention coverage file successfully read in!" << endl;


	/////////////////////////////////////////////////////////////////////////
	// 1.7.2. Fill out intervention structure

	intervention PNG_intven;

	for (int i = 0; i<N_cov_rounds; i++)
	{
		//////////////////////////////////////////////////////////////
		// LLINs

		if ((coverage[i][0] > -0.5) && (coverage[i][1] > -0.5))
		{
			PNG_intven.LLIN_year.push_back(coverage[i][0] * 365.0);
			PNG_intven.LLIN_cover.push_back(coverage[i][1]);
		}


		//////////////////////////////////////////////////////////////
		// IRS

		if ((coverage[i][0] > -0.5) && (coverage[i][2] > -0.5))
		{
			PNG_intven.IRS_year.push_back(coverage[i][0] * 365.0);
			PNG_intven.IRS_cover.push_back(coverage[i][2]);
		}


		//////////////////////////////////////////////////////////////
		// MDA - blood-stage drugs

		if ((coverage[i][0] > -0.5) && (coverage[i][3] > -0.5))
		{
			PNG_intven.MDA_BS_year.push_back(coverage[i][0] * 365.0);
			PNG_intven.MDA_BS_cover.push_back(coverage[i][3]);
			PNG_intven.MDA_BS_BSeff.push_back(coverage[i][4]);
			PNG_intven.MDA_BS_BSproph.push_back(coverage[i][5]);
		}


		//////////////////////////////////////////////////////////////
		// MDA - blood-stage drugs plus primaquine

		if ((coverage[i][0] > -0.5) && (coverage[i][6] > -0.5))
		{
			PNG_intven.MDA_PQ_year.push_back(coverage[i][0] * 365.0);
			PNG_intven.MDA_PQ_cover.push_back(coverage[i][6]);
			PNG_intven.MDA_PQ_BSeff.push_back(coverage[i][7]);
			PNG_intven.MDA_PQ_PQeff.push_back(coverage[i][8]);
			PNG_intven.MDA_PQ_BSproph.push_back(coverage[i][9]);
			PNG_intven.MDA_PQ_PQproph.push_back(coverage[i][10]);
			PNG_intven.MDA_PQ_CYP2D6.push_back((int) (coverage[i][11]));
		}


		//////////////////////////////////////////////////////////////
		// Front-line treatment - blood-stage drugs

		if ((coverage[i][0] > -0.5) && (coverage[i][12] > -0.5))
		{
			PNG_intven.BS_treat_year_on.push_back(coverage[i][0] * 365.0);
			PNG_intven.BS_treat_cover.push_back(coverage[i][12]);
			PNG_intven.BS_treat_BSeff.push_back(coverage[i][13]);
			PNG_intven.BS_treat_BSproph.push_back(coverage[i][14]);
			PNG_intven.BS_treat_year_off.push_back(coverage[i][15] * 365.0);
		}


		//////////////////////////////////////////////////////////////
		// Front-line treatment - primaquine

		if ((coverage[i][0] > -0.5) && (coverage[i][16] > -0.5))
		{
			PNG_intven.PQ_treat_year_on.push_back(coverage[i][0] * 365.0);
			PNG_intven.PQ_treat_cover.push_back(coverage[i][16]);
			PNG_intven.PQ_treat_PQcover.push_back(coverage[i][17]);
			PNG_intven.PQ_treat_BSeff.push_back(coverage[i][18]);
			PNG_intven.PQ_treat_PQeff.push_back(coverage[i][19]);
			PNG_intven.PQ_treat_BSproph.push_back(coverage[i][20]);
			PNG_intven.PQ_treat_PQproph.push_back(coverage[i][21]);
			PNG_intven.PQ_treat_CYP2D6.push_back((int) (coverage[i][22]));
			PNG_intven.PQ_treat_year_off.push_back(coverage[i][23] * 365.0);
		}
	}


	///////////////////////////////////////////////////////////////////////////
	//                                                                       //
	// 1.8. Initialise Population of individuals                             //
	//      Note that they begin with exponential age distribution           //
	//      and susceptible without immunity                                 //
	//                                                                       // 
	///////////////////////////////////////////////////////////////////////////

	Rcpp::Rcout << "Initialise population of individuals for simulation at equilbirium EIR of " << 365.0*Pv_mod_par.EIR_equil << endl;
	Rcpp::Rcout << endl;

	equi_pop_setup(&PNG_pop, &Pv_mod_par);

	Rcpp::Rcout << "Population of size " << PNG_pop.N_pop << " initialised!" << endl;
	Rcpp::Rcout << endl;


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

	PNG_sim.prev_2_10.resize(N_time);
	for (int i = 0; i<N_time; i++)
	{
		PNG_sim.prev_2_10[i].resize(11);
	}

	PNG_sim.EIR_t.resize(N_time);

	PNG_sim.LLIN_cov_t.resize(N_time);
	PNG_sim.IRS_cov_t.resize(N_time);
	PNG_sim.ACT_treat_t.resize(N_time);
	PNG_sim.PQ_treat_t.resize(N_time);
	PNG_sim.pregnant_t.resize(N_time);

	PNG_sim.A_par_mean_t.resize(N_time);
	PNG_sim.A_clin_mean_t.resize(N_time);


	//////////////////////////////////////////////////////
	//                                                  //
	// 1.10. Begin stochastic simulations               //
	//                                                  //
	////////////////////////////////////////////////////// 

	Rcpp::Rcout << "Starting model simulations......." << endl;

	model_simulator(&Pv_mod_par, &PNG_pop, &PNG_intven, &PNG_sim);

	Rcpp::Rcout << "Model simulations completed....." << endl;
	Rcpp::Rcout << endl;


	//////////////////////////////////////////////////////
	//                                                  //
	// 1.11. Output to file                             //
	//                                                  //
	////////////////////////////////////////////////////// 

	Rcpp::Rcout << "Start writing output to file......" << endl;
	Rcpp::Rcout << endl;

	ofstream output_Stream(output_File);

	for (int i = 0; i<N_time; i++)
	{
		output_Stream << PNG_sim.t_vec[i] << "\t";

		for (int k = 0; k<N_H_comp; k++)
		{
			output_Stream << PNG_sim.yH_t[i][k] << "\t";
		}

		for (int g = 0; g < N_spec; g++)
		{
			for (int k = 0; k < N_M_comp; k++)
			{
				output_Stream << PNG_sim.yM_t[i][g][k] << "\t";
			}
		}

		for (int k = 0; k<10; k++)
		{
			output_Stream << PNG_sim.prev_all[i][k] << "\t";
		}

		for (int k = 0; k<10; k++)
		{
			output_Stream << PNG_sim.prev_U5[i][k] << "\t";
		}

		for (int k = 0; k<10; k++)
		{
			output_Stream << PNG_sim.prev_2_10[i][k] << "\t";
		}

		output_Stream << PNG_sim.EIR_t[i] << "\t";
		output_Stream << PNG_sim.LLIN_cov_t[i] << "\t";
		output_Stream << PNG_sim.IRS_cov_t[i] << "\t";
		output_Stream << PNG_sim.ACT_treat_t[i] << "\t";
		output_Stream << PNG_sim.PQ_treat_t[i] << "\t";
		output_Stream << PNG_sim.pregnant_t[i] << "\t";

		output_Stream << PNG_sim.A_par_mean_t[i] << "\t";
		output_Stream << PNG_sim.A_clin_mean_t[i] << "\t";

		output_Stream << endl;
	}

	output_Stream.close();


	Rcpp::Rcout << "Output successfully written to file......" << endl;
	Rcpp::Rcout << endl;


	Rcpp::Rcout << "Time taken: " << ( (double) clock() - clock_time)/( (double) CLOCKS_PER_SEC ) << " seconds" << endl;


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

void mosq_derivs(const double t, double (&yM)[N_spec][N_M_comp], double (&dyMdt)[N_spec][N_M_comp], params* theta, population* POP)
{
	double Karry_seas_inv[N_spec];

	for (int g = 0; g < N_spec; g++)
	{
		Karry_seas_inv[g] = 1.0 / (theta->Karry[g] * ( theta->dry_seas[g] + (1 - theta->dry_seas[g])*pow(0.5 + 0.5*cos(0.01721421*(t - theta->t_peak_seas[g])), theta->kappa_seas[g])/ theta->denom_seas[g] ) );

		//Karry_seas_inv[g] = 1.0/theta->Karry[g];

		dyMdt[g][0] = POP->beta_VC[g] * (yM[g][3] + yM[g][4] + yM[g][5]) - yM[g][0] / theta->d_E_larvae - yM[g][0] * theta->mu_E0*(1.0 + (yM[g][0] + yM[g][1])*Karry_seas_inv[g]);
		dyMdt[g][1] = yM[g][0] / theta->d_E_larvae - yM[g][1] / theta->d_L_larvae - yM[g][1] * theta->mu_L0*(1.0 + theta->gamma_larvae*(yM[g][0] + yM[g][1])*Karry_seas_inv[g]);
		dyMdt[g][2] = yM[g][1] / theta->d_L_larvae - yM[g][2] / theta->d_pupae - yM[g][2] * theta->mu_P;
		dyMdt[g][3] = 0.5*yM[g][2] / theta->d_pupae - theta->lam_M[g] * yM[g][3] - POP->mu_M_VC[g] * yM[g][3];
		dyMdt[g][4] = + theta->lam_M[g] * yM[g][3] - theta->lam_S_M_track[g][0] * POP->exp_muM_tauM_VC[g] - POP->mu_M_VC[g] * yM[g][4];
		dyMdt[g][5] =                              + theta->lam_S_M_track[g][0] * POP->exp_muM_tauM_VC[g] - POP->mu_M_VC[g] * yM[g][5];
	}
}

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.2. Runge-Kutta 4 step updater for mosquito model                      //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void mosq_rk4(const double t, const double t_step_mosq, double (&yM)[N_spec][N_M_comp], params* theta, population* POP)
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
			theta->lam_S_M_track[g].push_back(theta->lam_M[g]* yM[g][3]);
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
			} else {
				HH.G6PD_def = 0;
			}
		} else {

			q_rand = genunf(0.0, 1.0);

			if(q_rand <= theta->G6PD_prev*theta->G6PD_prev)
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
		POP->prev_2_10[k] = 0.0;
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
		POP->prev_all[9]  = POP->prev_all[9]  + POP->people[n].ACT_new;
		POP->prev_all[10] = POP->prev_all[10] + POP->people[n].PQ_new;


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

			POP->prev_U5[6]  = POP->prev_U5[6]  + POP->people[n].I_PCR_new;
			POP->prev_U5[7]  = POP->prev_U5[7]  + POP->people[n].I_LM_new;
			POP->prev_U5[8]  = POP->prev_U5[8]  + POP->people[n].I_D_new;
			POP->prev_U5[9]  = POP->prev_U5[9]  + POP->people[n].ACT_new;
			POP->prev_U5[10] = POP->prev_U5[10] + POP->people[n].PQ_new;
		}

		//////////////////////////////////////////////
		//////////////////////////////////////////////
		// Summary - between 2 and 10's

		if (POP->people[n].age > 730.0 && POP->people[n].age < 3650.0)
		{
			////////////////////////////////////////
			// Prevalence

			POP->prev_2_10[0] = POP->prev_2_10[0] + 1;                                                            // Numbers - denominator
			POP->prev_2_10[1] = POP->prev_2_10[1] + POP->people[n].I_PCR + POP->people[n].I_LM
				                                + POP->people[n].I_D + POP->people[n].T;                          // PCR detectable infections
			POP->prev_2_10[2] = POP->prev_2_10[2] + POP->people[n].I_LM + POP->people[n].I_D + POP->people[n].T;    // LM detectable infections
			POP->prev_2_10[3] = POP->prev_2_10[3] + POP->people[n].I_D + POP->people[n].T;                          // Clinical episodes

			if (POP->people[n].Hyp > 0)
			{
				POP->prev_2_10[4] = POP->prev_2_10[4] + 1;                     // Hypnozoite positive

				POP->prev_2_10[5] = POP->prev_2_10[5] + POP->people[n].Hyp;    // Number of batches of hypnozoites
			}


			////////////////////////////////////////
			// Incidence

			POP->prev_2_10[6]  = POP->prev_2_10[6]  + POP->people[n].I_PCR_new;
			POP->prev_2_10[7]  = POP->prev_2_10[7]  + POP->people[n].I_LM_new;
			POP->prev_2_10[8]  = POP->prev_2_10[8]  + POP->people[n].I_D_new;
			POP->prev_2_10[9]  = POP->prev_2_10[9]  + POP->people[n].ACT_new;
			POP->prev_2_10[10] = POP->prev_2_10[10] + POP->people[n].PQ_new;
		}
	}


	//////////////////////////////
	// Intervention coverage

	POP->LLIN_cov_t = 0;
	POP->IRS_cov_t = 0;
	POP->ACT_treat_t = 0;
	POP->PQ_treat_t = 0;
	POP->pregnant_t = 0;

	for (int n = 0; n<POP->N_pop; n++)
	{
		POP->LLIN_cov_t  = POP->LLIN_cov_t  + POP->people[n].LLIN;
		POP->IRS_cov_t   = POP->IRS_cov_t   + POP->people[n].IRS;
		POP->ACT_treat_t = POP->ACT_treat_t + POP->people[n].ACT_new;
		POP->PQ_treat_t  = POP->PQ_treat_t + POP->people[n].PQ_new;
		POP->pregnant_t  = POP->pregnant_t  + POP->people[n].pregnant;
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
		//check for a user interupt every 100 timesteps
		if (i % 100 == 0) {
			Rcpp::checkUserInterrupt();
		}

		if (SIM->t_vec[i] / 365.0 - floor(SIM->t_vec[i] / 365.0) < 0.5*t_step / 365.0)
		{
			Rcpp::Rcout << "time = " << SIM->t_vec[i] / 365.0 << "\t" << 100.0*(SIM->t_vec[i] - SIM->t_vec[0]) / (double(t_step*SIM->N_time)) << "% complete" << endl;
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
			SIM->prev_U5[i][k]  = POP->prev_U5[k];
			SIM->prev_2_10[i][k] = POP->prev_2_10[k];
		}


		SIM->LLIN_cov_t[i]  = POP->LLIN_cov_t;
		SIM->IRS_cov_t[i]   = POP->IRS_cov_t;
		SIM->ACT_treat_t[i] = POP->ACT_treat_t;
		SIM->PQ_treat_t[i]  = POP->PQ_treat_t;
		SIM->pregnant_t[i]  = POP->pregnant_t;

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

	//////////////////////////////////////////////////////////
	// Intervention 1: LLINS

	for (int m = 0; m<INTVEN->LLIN_year.size(); m++)
	{
		if ((t > INTVEN->LLIN_year[m] - 0.5*t_step) &&
			(t < INTVEN->LLIN_year[m] + 0.51*t_step))
		{
			Rcpp::Rcout << "LLIN distribution" << endl;

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
			Rcpp::Rcout << "IRS distribution" << endl;

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
	// Intervention 3: MDA (blood-stage)

	for (int m = 0; m<INTVEN->MDA_BS_year.size(); m++)
	{
		if ((t > INTVEN->MDA_BS_year[m] - 0.5*t_step) &&
			(t < INTVEN->MDA_BS_year[m] + 0.51*t_step))
		{
			Rcpp::Rcout << "MDA (BS) distribution" << endl;

			theta->MDA_BS_cover   = INTVEN->MDA_BS_cover[m];
			theta->MDA_BS_BSeff   = INTVEN->MDA_BS_BSeff[m];
			theta->MDA_BS_BSproph = INTVEN->MDA_BS_BSproph[m];

			QQ = phi_inv(theta->MDA_BS_cover, 0.0, sqrt(1.0 + theta->sig_round_MDA*theta->sig_round_MDA));

			for (int n = 0; n<POP->N_pop; n++)
			{
				if (gennor(POP->people[n].zz_int[2], theta->sig_round_MDA) < QQ)
				{
					POP->people[n].ACT_new = 1;

					if (genunf(0.0, 1.0) < theta->MDA_BS_BSeff)
					{
						if (gennor(POP->people[n].zz_int[2], theta->sig_round_MDA) < QQ)
						{
							if (POP->people[n].S == 1    ) { POP->people[n].S = 0;     POP->people[n].P = 1; }
							if (POP->people[n].I_PCR == 1) { POP->people[n].I_PCR = 0; POP->people[n].P = 1;  }
							if (POP->people[n].I_LM == 1 ) { POP->people[n].I_LM = 0;  POP->people[n].P = 1;  }
							if (POP->people[n].I_D == 1  ) { POP->people[n].I_D = 0;   POP->people[n].T = 1;  }
						}
					}
				}
			}
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 4: MDA (blood-stage and liver-stage)

	for (int m = 0; m<INTVEN->MDA_PQ_year.size(); m++)
	{
		if ((t > INTVEN->MDA_PQ_year[m] - 0.5*t_step) &&
			(t < INTVEN->MDA_PQ_year[m] + 0.51*t_step))
		{
			Rcpp::Rcout << "MDA (BS+PQ) distribution" << endl;

			theta->MDA_PQ_cover   = INTVEN->MDA_PQ_cover[m];
			theta->MDA_PQ_BSeff   = INTVEN->MDA_PQ_BSeff[m];
			theta->MDA_PQ_PQeff   = INTVEN->MDA_PQ_PQeff[m];
			theta->MDA_PQ_BSproph = INTVEN->MDA_PQ_BSproph[m];
			theta->MDA_PQ_PQproph = INTVEN->MDA_PQ_PQproph[m];
			theta->MDA_PQ_CYP2D6  = INTVEN->MDA_PQ_CYP2D6[m];

			QQ = phi_inv(theta->MDA_PQ_cover, 0.0, sqrt(1.0 + theta->sig_round_MDA*theta->sig_round_MDA));

			for (int n = 0; n<POP->N_pop; n++)
			{
				if (gennor(POP->people[n].zz_int[3], theta->sig_round_MDA) < QQ)
				{
					POP->people[n].ACT_new = 1;

					if (genunf(0.0, 1.0) < theta->MDA_PQ_BSeff)
					{
						if (POP->people[n].S == 1    ) { POP->people[n].S = 0;     POP->people[n].P = 1; }
						if (POP->people[n].I_PCR == 1) { POP->people[n].I_PCR = 0; POP->people[n].P = 1; }
						if (POP->people[n].I_LM == 1 ) { POP->people[n].I_LM = 0;  POP->people[n].P = 1; }
						if (POP->people[n].I_D == 1  ) { POP->people[n].I_D = 0;   POP->people[n].T = 1; }
					}

					if( (POP->people[n].G6PD_def == 0) && (POP->people[n].pregnant == 0) && (POP->people[n].age > 180.0) )
					{
						POP->people[n].PQ_new = 1;

						if (theta->MDA_PQ_CYP2D6 == 0)    // Is CYP2D6 low metabolization a problem? No = 0, e.g. TQ; Otherwise Yes = 1, e.g. PQ
						{
							if (genunf(0.0, 1.0) < theta->MDA_PQ_PQeff)
							{
								POP->people[n].Hyp = 0;

								POP->people[n].PQ_proph = 1;
								POP->people[n].PQ_proph_timer = theta->PQ_treat_PQproph;
							}
						}else{
							if (POP->people[n].CYP2D6 == 0)          // If CYP2D6 low metabolization is a problem - it only effects the low metabolizers
							{
								if (genunf(0.0, 1.0) < theta->MDA_PQ_PQeff)
								{
									POP->people[n].Hyp = 0;

									POP->people[n].PQ_proph = 1;
									POP->people[n].PQ_proph_timer = theta->PQ_treat_PQproph;
								}
							}
						}

					}
				}
			}
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 7: first-line treatment (blood-stage)

	for (int m = 0; m<INTVEN->BS_treat_year_on.size(); m++)
	{
		if ((t > INTVEN->BS_treat_year_on[m] - 0.5*t_step) &&
			(t < INTVEN->BS_treat_year_on[m] + 0.51*t_step))
		{
			Rcpp::Rcout << "New front-line BS treatment" << endl;

			theta->BS_treat_cover   = INTVEN->BS_treat_cover[m];
			theta->BS_treat_BSeff   = INTVEN->BS_treat_BSeff[m];
			theta->BS_treat_BSproph = INTVEN->BS_treat_BSproph[m];

			theta->treat_cov = theta->BS_treat_cover;
			theta->treat_eff = theta->BS_treat_BSeff;
			theta->r_P = 1.0 / theta->BS_treat_BSproph;
		}
	}


	//////////////////////////////////////////////////////////
	// Switching back to baseline.

	for (int m = 0; m<INTVEN->BS_treat_year_on.size(); m++)
	{
		if ((t > INTVEN->BS_treat_year_off[m] - 0.5*t_step) &&
			(t < INTVEN->BS_treat_year_off[m] + 0.51*t_step))
		{
			Rcpp::Rcout << "End of changing front-line BS treatment" << endl;

			theta->treat_cov = theta->BS_treat_cov_base;
			theta->treat_eff = theta->BS_treat_eff_base;
			theta->r_P = 1.0 / theta->BS_treat_BSproph_base;
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 8: first-line treatment (primaquine)

	for (int m = 0; m<INTVEN->PQ_treat_year_on.size(); m++)
	{
		if ((t > INTVEN->PQ_treat_year_on[m] - 0.5*t_step) &&
			(t < INTVEN->PQ_treat_year_on[m] + 0.51*t_step))
		{
			Rcpp::Rcout << "New front-line PQ treatment" << endl;

			theta->PQ_treat_cover   = INTVEN->PQ_treat_cover[m];
			theta->PQ_treat_PQcover = INTVEN->PQ_treat_PQcover[m];
			theta->PQ_treat_BSeff   = INTVEN->PQ_treat_BSeff[m];
			theta->PQ_treat_PQeff   = INTVEN->PQ_treat_PQeff[m];
			theta->PQ_treat_BSproph = INTVEN->PQ_treat_BSproph[m];
			theta->PQ_treat_PQproph = INTVEN->PQ_treat_PQproph[m];
			theta->PQ_treat_CYP2D6  = INTVEN->PQ_treat_CYP2D6[m];

			theta->treat_cov = theta->PQ_treat_cover;
			theta->treat_eff = theta->PQ_treat_BSeff;
			theta->r_P = 1.0 / theta->PQ_treat_BSproph;
		}
	}


	//////////////////////////////////////////////////////////
	// Switching back to baseline.

	for (int m = 0; m<INTVEN->PQ_treat_year_on.size(); m++)
	{
		if ((t > INTVEN->PQ_treat_year_off[m] - 0.5*t_step) &&
			(t < INTVEN->PQ_treat_year_off[m] + 0.51*t_step))
		{
			Rcpp::Rcout << "End of changing front-line PQ treatment" << endl;

			theta->PQ_treat_cover = 0.0;
			theta->PQ_treat_PQeff = 0.0;
			theta->PQ_treat_BSproph = 10.0;
			theta->PQ_treat_PQproph = 10.0;
			theta->PQ_treat_CYP2D6 = 1;

			theta->treat_cov = theta->BS_treat_cov_base;
			theta->treat_eff = theta->BS_treat_eff_base;
			theta->r_P = 1.0 / theta->BS_treat_BSproph_base;
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
			MM[0 * (K_max + 1) + k1][0 * (K_max + 1) + k2] = - lam_eq[i][j] * theta->D_MAT[k1][k2] - theta->ff*theta->K_MAT[k1][k2]
				                                             + theta->gamma_L*theta->L_MAT[k1][k2] - theta->mu_H*theta->D_MAT[k1][k2] - POP->r_age[i] * theta->D_MAT[k1][k2];
			MM[0 * (K_max + 1) + k1][1 * (K_max + 1) + k2] = + r_PCR_eq[i][j][k2] * theta->D_MAT[k1][k2];
			MM[0 * (K_max + 1) + k1][5 * (K_max + 1) + k2] = + theta->r_P*theta->D_MAT[k1][k2];

			MM[1 * (K_max + 1) + k1][0 * (K_max + 1) + k2] = + lam_eq[i][j] * (1.0 - phi_LM_eq[i][j][k2])*theta->OD_MAT[k1][k2] + theta->ff*(1.0 - phi_LM_eq[i][j][k2])*theta->K_MAT[k1][k2];
			MM[1 * (K_max + 1) + k1][1 * (K_max + 1) + k2] = - lam_eq[i][j] * theta->D_MAT[k1][k2] - theta->ff*theta->K_MAT[k1][k2] - r_PCR_eq[i][j][k2] * theta->D_MAT[k1][k2]
				                                             + lam_eq[i][j] * (1.0 - phi_LM_eq[i][j][k2])*theta->OD_MAT[k1][k2] + theta->ff*(1.0 - phi_LM_eq[i][j][k2])*theta->K_MAT[k1][k2]
				                                             + theta->gamma_L*theta->L_MAT[k1][k2] - theta->mu_H*theta->D_MAT[k1][k2] - POP->r_age[i] * theta->D_MAT[k1][k2];
			MM[1 * (K_max + 1) + k1][2 * (K_max + 1) + k2] = + theta->r_LM*theta->D_MAT[k1][k2];

			MM[2 * (K_max + 1) + k1][0 * (K_max + 1) + k2] = + lam_eq[i][j] * phi_LM_eq[i][j][k2] * (1.0 - phi_D_eq[i][j][k2])*theta->OD_MAT[k1][k2] + theta->ff*phi_LM_eq[i][j][k2] * (1.0 - phi_D_eq[i][j][k2])*theta->K_MAT[k1][k2];
			MM[2 * (K_max + 1) + k1][1 * (K_max + 1) + k2] = + lam_eq[i][j] * phi_LM_eq[i][j][k2] * (1.0 - phi_D_eq[i][j][k2])*theta->OD_MAT[k1][k2] + theta->ff*phi_LM_eq[i][j][k2] * (1.0 - phi_D_eq[i][j][k2])*theta->K_MAT[k1][k2];
			MM[2 * (K_max + 1) + k1][2 * (K_max + 1) + k2] = - lam_eq[i][j] * theta->D_MAT[k1][k2] - theta->ff*theta->K_MAT[k1][k2] - theta->r_LM*theta->D_MAT[k1][k2]
				                                             + lam_eq[i][j] * (1.0 - phi_D_eq[i][j][k2])*theta->OD_MAT[k1][k2] + theta->ff*(1.0 - phi_D_eq[i][j][k2])*theta->K_MAT[k1][k2]
				                                             + theta->gamma_L*theta->L_MAT[k1][k2] - theta->mu_H*theta->D_MAT[k1][k2] - POP->r_age[i] * theta->D_MAT[k1][k2];
			MM[2 * (K_max + 1) + k1][3 * (K_max + 1) + k2] = + theta->r_D*theta->D_MAT[k1][k2];

			MM[3 * (K_max + 1) + k1][0 * (K_max + 1) + k2] = + lam_eq[i][j] * phi_LM_eq[i][j][k2] * phi_D_eq[i][j][k2] * (1.0 - theta->treat_cov*theta->treat_eff)*theta->OD_MAT[k1][k2] + theta->ff*phi_LM_eq[i][j][k2] * phi_D_eq[i][j][k2] * (1.0 - theta->treat_cov*theta->treat_eff)*theta->K_MAT[k1][k2];
			MM[3 * (K_max + 1) + k1][1 * (K_max + 1) + k2] = + lam_eq[i][j] * phi_LM_eq[i][j][k2] * phi_D_eq[i][j][k2] * (1.0 - theta->treat_cov*theta->treat_eff)*theta->OD_MAT[k1][k2] + theta->ff*phi_LM_eq[i][j][k2] * phi_D_eq[i][j][k2] * (1.0 - theta->treat_cov*theta->treat_eff)*theta->K_MAT[k1][k2];
			MM[3 * (K_max + 1) + k1][2 * (K_max + 1) + k2] = + lam_eq[i][j] * phi_D_eq[i][j][k2] * (1.0 - theta->treat_cov*theta->treat_eff)*theta->OD_MAT[k1][k2] + theta->ff*phi_D_eq[i][j][k2] * (1.0 - theta->treat_cov*theta->treat_eff)*theta->K_MAT[k1][k2];
			MM[3 * (K_max + 1) + k1][3 * (K_max + 1) + k2] = - lam_eq[i][j] * theta->D_MAT[k1][k2] - theta->r_D*theta->D_MAT[k1][k2] + lam_eq[i][j] * theta->OD_MAT[k1][k2]
				                                             + theta->gamma_L*theta->L_MAT[k1][k2] - theta->mu_H*theta->D_MAT[k1][k2] - POP->r_age[i] * theta->D_MAT[k1][k2];

			MM[4 * (K_max + 1) + k1][0 * (K_max + 1) + k2] = + lam_eq[i][j] * phi_LM_eq[i][j][k2] * phi_D_eq[i][j][k2] * theta->treat_cov*theta->treat_eff*theta->OD_MAT[k1][k2] + theta->ff*phi_LM_eq[i][j][k2] * phi_D_eq[i][j][k2] * theta->treat_cov*theta->treat_eff*theta->K_MAT[k1][k2];
			MM[4 * (K_max + 1) + k1][1 * (K_max + 1) + k2] = + lam_eq[i][j] * phi_LM_eq[i][j][k2] * phi_D_eq[i][j][k2] * theta->treat_cov*theta->treat_eff*theta->OD_MAT[k1][k2] + theta->ff*phi_LM_eq[i][j][k2] * phi_D_eq[i][j][k2] * theta->treat_cov*theta->treat_eff*theta->K_MAT[k1][k2];
			MM[4 * (K_max + 1) + k1][2 * (K_max + 1) + k2] = + lam_eq[i][j] * phi_D_eq[i][j][k2] * theta->treat_cov*theta->treat_eff*theta->OD_MAT[k1][k2] + theta->ff*phi_D_eq[i][j][k2] * theta->treat_cov*theta->treat_eff*theta->K_MAT[k1][k2];
			MM[4 * (K_max + 1) + k1][4 * (K_max + 1) + k2] = - lam_eq[i][j] * theta->D_MAT[k1][k2] - theta->r_T*theta->D_MAT[k1][k2] + lam_eq[i][j] * theta->OD_MAT[k1][k2]
				                                             + theta->gamma_L*theta->L_MAT[k1][k2] - theta->mu_H*theta->D_MAT[k1][k2] - POP->r_age[i] * theta->D_MAT[k1][k2];

			MM[5 * (K_max + 1) + k1][4 * (K_max + 1) + k2] = + theta->r_T*theta->D_MAT[k1][k2];
			MM[5 * (K_max + 1) + k1][5 * (K_max + 1) + k2] = - lam_eq[i][j] * theta->D_MAT[k1][k2] - theta->r_P*theta->D_MAT[k1][k2] + lam_eq[i][j] * theta->OD_MAT[k1][k2]
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

	Rcpp::Rcout << "Testing yH sum......" << endl;

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


	Rcpp::Rcout << "yH_sum = " << yH_sum << endl;
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


	if( (theta->Prop_mosq[0] < 1.0e-10) && (theta->Prop_mosq[1] < 1.0e-10) )
	{
		theta->mm_0[0] = 0.0;
		theta->mm_0[1] = 0.0;
		theta->mm_0[2] = theta->EIR_equil /( theta->aa[2] *I_M_star[2] );

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
		if (g == 0) { Rcpp::Rcout << "An. farauti:  " << 100.0 * theta->Prop_mosq[0] << "%" << endl; }
		if (g == 1) { Rcpp::Rcout << "An. punctulatus:  " << 100.0 * theta->Prop_mosq[1] << "%" << endl; }
		if (g == 2) { Rcpp::Rcout << "An. koliensis:  " << 100.0 * theta->Prop_mosq[2] << "%" << endl; }

		Rcpp::Rcout << "EL_M  " << POP->yM[g][0] <<  endl;
		Rcpp::Rcout << "LL_M  " << POP->yM[g][1] <<  endl;
		Rcpp::Rcout << "P_M  " << POP->yM[g][2]  <<  endl;
		Rcpp::Rcout << "S_M  " << POP->yM[g][3]  <<  endl;
		Rcpp::Rcout << "E_M  " << POP->yM[g][4]  <<  endl;
		Rcpp::Rcout << "I_M  " << POP->yM[g][5]  <<  endl;

		Rcpp::Rcout << "lam_M = " << theta->lam_M[g] << endl;

		Rcpp::Rcout << "I_M = " << POP->yM[g][5] << endl;
	
		Rcpp::Rcout << "mm = " << theta->mm_0[g] << endl;

		Rcpp::Rcout << endl;
	}
	Rcpp::Rcout << endl;

	Rcpp::Rcout << "lam_H = " << theta->bb*theta->EIR_equil << endl;

	double EIR_out = 0.0;
	for (int g = 0; g < N_spec; g++)
	{
		EIR_out = EIR_out + 365.0*theta->aa[g] * POP->yM[g][5];
	}

	Rcpp::Rcout << "EIR = " << EIR_out << endl;


	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	// 3.7.4. Proportion in each age and heterogeneity stratified category

	//////////////////////////////////////////
	// Fill out vector of lagged lam_M*S_M

	theta->lam_S_M_track.resize(N_spec);

	for( int g=0; g<N_spec; g++ )
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
		} else {
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

	Rcpp::Rcout << "Equilibrium set up......." << endl;

	double S_ind = 0.0, I_PCR_ind = 0.0, I_LM_ind = 0.0, I_D_ind = 0.0, T_ind = 0.0, P_ind = 0.0;
	double S_eqq = 0.0, I_PCR_eqq = 0.0, I_LM_eqq = 0.0, I_D_eqq = 0.0, T_eqq = 0.0, P_eqq = 0.0;

	for (int n = 0; n<POP->N_pop; n++)
	{
		S_ind     = S_ind     + POP->people[n].S;
		I_PCR_ind = I_PCR_ind + POP->people[n].I_PCR;
		I_LM_ind  = I_LM_ind  + POP->people[n].I_LM;
		I_D_ind   = I_D_ind   + POP->people[n].I_D;
		T_ind     = T_ind     + POP->people[n].T;
		P_ind     = P_ind     + POP->people[n].P;
	}


	for (int i = 0; i<N_age; i++)
	{
		for (int j = 0; j<N_het; j++)
		{
			for (int k = 0; k < (K_max + 1); k++)
			{
				S_eqq     = S_eqq     + yH_eq[i][j][k][0];
				I_PCR_eqq = I_PCR_eqq + yH_eq[i][j][k][1];
				I_LM_eqq  = I_LM_eqq  + yH_eq[i][j][k][2];
				I_D_eqq   = I_D_eqq   + yH_eq[i][j][k][3];
				T_eqq     = T_eqq     + yH_eq[i][j][k][4];
				P_eqq     = P_eqq     + yH_eq[i][j][k][5];
			}
		}
	}


	Rcpp::Rcout << "S = " << ((double)S_ind) / POP->N_pop << "\t" << S_eqq << endl;
	Rcpp::Rcout << "I_PCR = " << ((double)I_PCR_ind) / POP->N_pop << "\t" << I_PCR_eqq << endl;
	Rcpp::Rcout << "I_LM = " << ((double)I_LM_ind) / POP->N_pop << "\t" << I_LM_eqq << endl;
	Rcpp::Rcout << "I_D = " << ((double)I_D_ind) / POP->N_pop << "\t" << I_D_eqq << endl;
	Rcpp::Rcout << "T = " << ((double)T_ind) / POP->N_pop << "\t" << T_eqq << endl;
	Rcpp::Rcout << "P = " << ((double)P_ind) / POP->N_pop << "\t" << P_eqq << endl;

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
	ACT_new = 0;
	PQ_new = 0;


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
			theta.S_move[2] = theta.phi_LM*theta.phi_D*(1.0 - theta.treat_cov);      // Move to I_D    //  lam_H_lag*theta.phi_LM*(1.0-theta.phi_D)*(1.0-theta.treat_cov)/theta.S_out;
			theta.S_move[3] = theta.phi_LM*theta.phi_D*theta.treat_cov;              // Move to T      //  lam_H_lag*theta.phi_LM*(1.0-theta.phi_D)*theta.treat_cov/theta.S_out;

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

				//Rcpp::Rcout << theta.PQ_treat_cover << "\t" << theta.PQ_treat_eff << endl;
				//Rcpp::Rcout << G6PD_def << "\t" << CYP2D6 << "\t" << pregnant << "\t" << age << endl;


				if (theta.PQ_treat_PQcover > 0.0)
				{
					if ((G6PD_def == 0) && (pregnant == 0) && (age > 180.0))
					{
						if (genunf(0.0, 1.0) < theta.PQ_treat_PQcover)
						{
							PQ_new = 1;


							if (theta.PQ_treat_CYP2D6 == 0)   // Case where CYP2D6 low met is not be a problem (e.g. TQ) 
							{
								if (genunf(0.0, 1.0) < theta.PQ_treat_PQeff)
								{
									Hyp = 0;
									PQ_proph = 1;
								}
							}else{                            // Otherwise CYP2D6 low met might be a problem (e.g. PQ)
								if (CYP2D6 == 0)
								{
									if (genunf(0.0, 1.0) < theta.PQ_treat_PQeff)
									{
										Hyp = 0;
										PQ_proph = 1;
									}
								}
							}
						}

						
					}
				}

				I_PCR_new = 1;
				I_LM_new = 1;
				I_D_new = 1;
				ACT_new = 1;

				S = 0;
				if (genunf(0.0, 1.0) < theta.treat_eff)
				{
					T = 1;
				} else {
					I_D = 1;
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
			theta.I_PCR_move[3] = lam_H_lag*theta.phi_LM*theta.phi_D*(1.0 - theta.treat_cov) / theta.I_PCR_out;        // Move to I_D
			theta.I_PCR_move[4] = lam_H_lag*theta.phi_LM*theta.phi_D*theta.treat_cov / theta.I_PCR_out;                // Move to T

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


				if (theta.PQ_treat_PQcover > 0.0)
				{
					if ((G6PD_def == 0) && (pregnant == 0) && (age > 180.0))
					{
						if (genunf(0.0, 1.0) < theta.PQ_treat_PQcover)
						{
							PQ_new = 1;

							if (theta.PQ_treat_CYP2D6 == 0)   // Case where CYP2D6 low met is not be a problem (e.g. TQ) 
							{
								if (genunf(0.0, 1.0) < theta.PQ_treat_PQeff)
								{
									Hyp = 0;
									PQ_proph = 1;
								}
							}
							else {                            // Otherwise CYP2D6 low met might be a problem (e.g. PQ)
								if (CYP2D6 == 0)
								{
									if (genunf(0.0, 1.0) < theta.PQ_treat_PQeff)
									{
										Hyp = 0;
										PQ_proph = 1;
									}
								}
							}
						}

					}
				}


				I_PCR_new = 1;
				I_LM_new  = 1;
				I_D_new   = 1;
				ACT_new   = 1;

				I_PCR = 0;
				if (genunf(0.0, 1.0) < theta.treat_eff)
				{
					T = 1;
				}else {
					I_D = 1;
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
			theta.I_LM_move[2] = lam_H_lag*theta.phi_D*(1.0 - theta.treat_cov) / theta.I_LM_out;     // Move to I_D
			theta.I_LM_move[3] = lam_H_lag*theta.phi_D*theta.treat_cov / theta.I_LM_out;             // Move to T

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
				I_D_new  = 1;

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

				if (theta.PQ_treat_PQcover > 0.0)
				{
					if ((G6PD_def == 0) && (pregnant == 0) && (age > 180.0))
					{
						if (genunf(0.0, 1.0) < theta.PQ_treat_PQcover)
						{
							PQ_new = 1;

							if (theta.PQ_treat_CYP2D6 == 0)   // Case where CYP2D6 low met is not be a problem (e.g. TQ) 
							{
								if (genunf(0.0, 1.0) < theta.PQ_treat_PQeff)
								{
									Hyp = 0;
									PQ_proph = 1;
								}
							} else {                            // Otherwise CYP2D6 low met might be a problem (e.g. PQ)
								if (CYP2D6 == 0)
								{
									if (genunf(0.0, 1.0) < theta.PQ_treat_PQeff)
									{
										Hyp = 0;
										PQ_proph = 1;
									}
								}
							}
						}


					}
				}

				I_PCR_new = 1;
				I_LM_new = 1;
				I_D_new = 1;
				ACT_new = 1;

				I_LM = 0;
				if (genunf(0.0, 1.0) < theta.treat_eff)
				{
					T = 1;
				} else {
					I_D = 1;
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
		} else {
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

