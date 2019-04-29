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
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

#include "Params.hpp"

#include <iostream>
#include <fstream>
#include <cmath>


////////////////////////////////////////////////////////////
//                                                        //
//  Function declarations                                 //
//                                                        //
////////////////////////////////////////////////////////////

double gammln(const double xx);


/////////////////////////////////////
// Read parameters from input files
void Params::read(const char *parameter_File, const char *mosquito_File[N_spec_max])
{
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

    parameter_Stream >> discard >> N_pop >> discard;               // Number of participants

    parameter_Stream >> discard >> EIR_equil >> discard;        // EIR at equilibrium

    for (int g = 0; g < N_spec_max; g++)
    {
        parameter_Stream >> discard >> Prop_mosq[g] >> discard; // Proportions of each mosquito species
    }

    parameter_Stream >> discard >> time_start >> discard;                  // Start time for simulation
    parameter_Stream >> discard >> time_end >> discard;                    // End time for simulation
    parameter_Stream >> discard >> burnin_time >> discard;                 // End time for simulation


    /////////////////////////////////
    // Human demography
    // Note that the mean age isn't technically the mean 
    // due to the effects of truncation at maximum age

    parameter_Stream >> discard >> age_mean >> discard;    // mean age of human population
    parameter_Stream >> discard >> age_max >> discard;     // maximum age in human population

    mu_H = 1.0 / age_mean;                      // human death rate

    het_max = 100.0;                                       // Maximum relative heterogeneity value


    ////////////////////////////////////////////////////////
    // Heterogeneity in exposure and age-dependent biting

    parameter_Stream >> discard >> age_0 >> discard;       // age-dependent biting parameter
    parameter_Stream >> discard >> rho_age >> discard;     // age-dependent biting parameter
    parameter_Stream >> discard >> sig_het >> discard;     // heterogeneity in exposure - standard deviation on a log scale


    /////////////////////////////////
    // Transmission probabilities

    parameter_Stream >> discard >> bb >> discard;          // mosquito -> human transmission probability
    parameter_Stream >> discard >> c_PCR >> discard;       // human -> mosquito transmission probability (PCR detectable)
    parameter_Stream >> discard >> c_LM >> discard;        // human -> mosquito transmission probability (LM detectable)
    parameter_Stream >> discard >> c_D >> discard;         // human -> mosquito transmission probability (clinical disease)
    parameter_Stream >> discard >> c_T >> discard;         // human -> mosquito transmission probability (during treatment)


    ////////////////////////////////
    // Human recovery paramters

    parameter_Stream >> discard >> d_latent >> discard;             // latent period in liver
    parameter_Stream >> discard >> r_LM >> discard;                 // rate of recovery from LM detectable infection
    parameter_Stream >> discard >> r_D >> discard;                  // rate of recovery from symptomatic disease
    parameter_Stream >> discard >> r_T >> discard;                  // rate of progression through treatment 
    parameter_Stream >> discard >> d_PCR_min >> discard;            // minimum duration of PCR-detectable infection - full immunity
    parameter_Stream >> discard >> d_PCR_max >> discard;            // maximum duration of PCR-detectable infection - no immunity

    parameter_Stream >> discard >> BS_treat_BScover_base >> discard;      // proportion of episodes of symptomatic disease treated (baseline)
    parameter_Stream >> discard >> BS_treat_BSeff_base >> discard;      // efficacy of front-line treatment (baseline)
    parameter_Stream >> discard >> BS_treat_BSproph_base >> discard;  // duration of prophylaxis of front-line treatment (baseline)

    parameter_Stream >> discard >> A_PCR_50pc >> discard;           // PCR_detectable infection scale parameter
    parameter_Stream >> discard >> K_PCR >> discard;                // PCR_detectable infection shape parameter


    H_track = int(d_latent / t_step);                       // Number of time steps for duration of latency

    treat_BScover = BS_treat_BScover_base;                  // Treatment coverage
    treat_BSeff   = BS_treat_BSeff_base;                    // Efficacy of treatment
    r_P           = 1.0 / BS_treat_BSproph_base;            // rate of recovery from prophylaxis
    treat_PQavail = MDA_PQ_PQavail;                         // PQ availability

    PQ_treat_BScover     = 0.0;
    PQ_treat_BSeff       = 0.0;
    PQ_treat_BSproph     = 10.0;
    PQ_treat_PQavail     = 0.0;
    PQ_treat_PQeff       = 0.0;
    PQ_treat_PQproph     = 0.0;
    PQ_treat_G6PD_risk   = 1;
    PQ_treat_CYP2D6_risk = 1;
    PQ_treat_preg_risk   = 1;
    PQ_treat_low_age     = 180.0;


    //////////////////
    // temporary setting treatment cov to zero
    //
    //treat_cov_eff = 0.0;


    /////////////////////////////////
    // Blood-stage immunity paramters

    parameter_Stream >> discard >> u_par >> discard;         // scale paramter for acquisition of blood-stage immunity
    parameter_Stream >> discard >> r_par >> discard;         // rate of decay of blood-stage immunity

    parameter_Stream >> discard >> phi_LM_max >> discard;    // probability of blood-stage infection with no immunity 
    parameter_Stream >> discard >> phi_LM_min >> discard;    // probability of blood-stage infection with maximum immunity
    parameter_Stream >> discard >> A_LM_50pc >> discard;     // blood-stage immunity scale parameter
    parameter_Stream >> discard >> K_LM >> discard;          // blood-stage immunity shape parameter


    /////////////////////////////////
    // Clinical immunity paramters

    parameter_Stream >> discard >> u_clin >> discard;       // scale paramter for acquisition of blood-stage immunity
    parameter_Stream >> discard >> r_clin >> discard;       // rate of decay of clinical immunity

    parameter_Stream >> discard >> phi_D_max >> discard;    // probability of clinical episode with no immunity
    parameter_Stream >> discard >> phi_D_min >> discard;    // probability of clinical episode with maximum immunity
    parameter_Stream >> discard >> A_D_50pc >> discard;     // clinical immunity scale parameter
    parameter_Stream >> discard >> K_D >> discard;          // clinical immunity shape parameter


    /////////////////////////////////////////////
    // maternal immunity

    parameter_Stream >> discard >> P_mat >> discard;        // New-born immunity relative to mother's
    parameter_Stream >> discard >> d_mat >> discard;        // Inverse of decay rate of maternal immunity


    /////////////////////////////////
    // Relapse paramters

    parameter_Stream >> discard >> ff >> discard;           // relapse rate
    parameter_Stream >> discard >> gamma_L >> discard;      // liver clearance rate


    ////////////////////////////////
    // Human genotype prevalences

    parameter_Stream >> discard >> G6PD_prev >> discard;    // prevalence of G6PD deficiency
    parameter_Stream >> discard >> CYP2D6_prev >> discard;  // prevalence of CYP2D6 phenotype 


    ////////////////////////////////////////////////////////
    // Intervention distribution parameters

    parameter_Stream >> discard >> rho_round_LLIN >> discard;
    parameter_Stream >> discard >> rho_round_IRS >> discard;
    parameter_Stream >> discard >> rho_round_MDA >> discard;

    parameter_Stream >> discard >> rho_LLIN_IRS >> discard;
    parameter_Stream >> discard >> rho_MDA_IRS >> discard;
    parameter_Stream >> discard >> rho_MDA_LLIN >> discard;


    sig_round_LLIN = sqrt((1.0 - rho_round_LLIN) / rho_round_LLIN);
    sig_round_IRS = sqrt((1.0 - rho_round_IRS) / rho_round_IRS);
    sig_round_MDA = sqrt((1.0 - rho_round_MDA) / rho_round_MDA);


    //////////////////////////////////////////////
    // TO DO - no need to have different possible correlations
    //         for the various MDA interventions - can assume access
    //         to individuals is all the same

    V_int[0][0] = 1.0;
    V_int[0][1] = rho_LLIN_IRS;
    V_int[0][2] = rho_MDA_LLIN;
    V_int[0][3] = rho_MDA_LLIN;
    V_int[0][4] = rho_MDA_LLIN;
    V_int[0][5] = rho_MDA_LLIN;

    V_int[1][0] = rho_LLIN_IRS;
    V_int[1][1] = 1.0;
    V_int[1][2] = rho_MDA_IRS;
    V_int[1][3] = rho_MDA_IRS;
    V_int[1][4] = rho_MDA_IRS;
    V_int[1][5] = rho_MDA_IRS;

    V_int[2][0] = rho_MDA_LLIN;
    V_int[2][1] = rho_MDA_IRS;
    V_int[2][2] = 1.0;
    V_int[2][3] = rho_round_MDA;
    V_int[2][4] = rho_round_MDA;
    V_int[2][5] = rho_round_MDA;

    V_int[3][0] = rho_MDA_LLIN;
    V_int[3][1] = rho_MDA_IRS;
    V_int[3][2] = rho_round_MDA;
    V_int[3][3] = 1.0;
    V_int[3][4] = rho_round_MDA;
    V_int[3][5] = rho_round_MDA;

    V_int[4][0] = rho_MDA_LLIN;
    V_int[4][1] = rho_MDA_IRS;
    V_int[4][2] = rho_round_MDA;
    V_int[4][3] = rho_round_MDA;
    V_int[4][4] = 1.0;
    V_int[4][5] = rho_round_MDA;

    V_int[5][0] = rho_MDA_LLIN;
    V_int[5][1] = rho_MDA_IRS;
    V_int[5][2] = rho_round_MDA;
    V_int[5][3] = rho_round_MDA;
    V_int[5][4] = rho_round_MDA;
    V_int[5][5] = 1.0;


    /////////////////////////////////////////////////////////
    // We need to make a dummy covariance matrix for intervention
    // distribution because of the way genmn works

    for (int p = 0; p<N_int; p++)
    {
        for (int q = 0; q<N_int; q++)
        {
            V_int_dummy[p][q] = V_int[p][q];
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

        mosquito_Stream >> discard >> mu_M[g] >> discard;      // mosquito death rate
        mosquito_Stream >> discard >> tau_M[g] >> discard;     // duration of sporogony

        M_track = int(tau_M[g] * mosq_steps / t_step);


        /////////////////////////////////
        // Larval paramters
        // From White et al (2011) P&V

        mosquito_Stream >> discard >> d_E_larvae >> discard;      // Development time of early larval instars
        mosquito_Stream >> discard >> d_L_larvae >> discard;      // Development time of late larval instars
        mosquito_Stream >> discard >> d_pupae >> discard;         // Development time of pupae
        mosquito_Stream >> discard >> mu_E0 >> discard;           // Mortality rate of early larval instars (low density)
        mosquito_Stream >> discard >> mu_L0 >> discard;           // Mortality rate of late larval instars (low density)
        mosquito_Stream >> discard >> mu_P >> discard;            // Mortality rate of pupae
        mosquito_Stream >> discard >> beta_larvae >> discard;     // Number of eggs laid per day per mosquito
        mosquito_Stream >> discard >> gamma_larvae >> discard;    // Effect of density dependence on late instars relative to early instars

        omega_larvae[g] = gamma_larvae*mu_L0 / mu_E0 - d_E_larvae / d_L_larvae + (gamma_larvae - 1.0)*mu_L0*d_E_larvae;
        omega_larvae[g] = - 0.5*omega_larvae[g] + sqrt(0.25*omega_larvae[g] * omega_larvae[g] + 0.5*gamma_larvae*beta_larvae*mu_L0*d_E_larvae /
                                       (mu_E0*mu_M[g] * d_L_larvae*(1.0 + d_pupae*mu_P)));


        /////////////////////////////////
        // Seasonality paramters
        // Denominator for seasonality - see Griffin (2015) PLoS Comp Biol

        mosquito_Stream >> discard >> dry_seas[g] >> discard;      // Proportion of dry season transmission compared to mean 
        mosquito_Stream >> discard >> kappa_seas[g] >> discard;    // Shape parameter for seasonality
        mosquito_Stream >> discard >> t_peak_seas[g] >> discard;   // Timing of peak for seasonal transmission

        denom_seas[g] = exp(gammln(0.5) + gammln(kappa_seas[g] + 0.5) - gammln(kappa_seas[g] + 1.0)) / 3.14159265359;


        //////////////////////////////////////////////
        // Entomology paramters

        mosquito_Stream >> discard >> Q_0[g] >> discard;           // Human Blood Index (proportion of blood meals taken on humans)
        mosquito_Stream >> discard >> CHI_endo[g] >> discard;      // Endophily - proportion of mosquitoes resting indoors after feeding (no intervention)
        mosquito_Stream >> discard >> PSI_indoors[g] >> discard;   // Proportion of bites taken on humans indoors
        mosquito_Stream >> discard >> PSI_bed[g] >> discard;       // Proportion of bites taken on humans in bed

        mosquito_Stream >> discard >> delta_1 >> discard;          // Time spent foraging for a blood meal
        mosquito_Stream >> discard >> delta >> discard;            // Duration of gonotrophic cycle

        delta_2 = delta - delta_1;

        p_1[g] = exp(-mu_M[g] * delta_1);
        p_2[g] = exp(-mu_M[g] * delta_2);

        aa[g] = Q_0[g] / (delta_1 + delta_2);

        eps_max[g] = beta_larvae*(exp(delta*mu_M[g]) - 1.0) / mu_M[g];


        //////////////////////////////////////////////
        // LLIN paramters

        mosquito_Stream >> discard >> LLIN_half_life >> discard;

        mosquito_Stream >> discard >> PYR_half_life >> discard;

        mosquito_Stream >> discard >> r_LLIN_0[g] >> discard;       // Probability mosquito repelled (with full insecticide activity)
        mosquito_Stream >> discard >> r_LLIN_net[g] >> discard;     // Probability mosquito repelled due to barrier effect of net (no insecticide)
        mosquito_Stream >> discard >> d_LLIN_0[g] >> discard;       // Probability mosquito dies during feeding attempt


        P_LLIN_loss = 1.0 - exp(-t_step*log(2.0) / LLIN_half_life);   // Probability of losing LLIN in a time step
        PYR_decay = log(2.0) / PYR_half_life;                            // Rate of pyrethroid decay

        s_LLIN_0[g] = 1.0 - r_LLIN_0[g] - d_LLIN_0[g];     // Probability mosquito feeds successfully


        ////////////////////////////////////////////////////////
        // IRS parameters

        mosquito_Stream >> discard >> IRS_half_life >> discard;    // IRS insecticide half-life

        mosquito_Stream >> discard >> r_IRS_0[g] >> discard;          // IRS repellency
        mosquito_Stream >> discard >> d_IRS_0[g] >> discard;          // IRS death

        IRS_decay = log(2.0) / IRS_half_life;         // IRS decay rate

        s_IRS_0[g] = 1.0 - d_IRS_0[g] - s_IRS_0[g];   // Feeding success of mosquito on IRS protected person
    }

    cout << "Mosquito parameter values read in from file!" << endl;
    cout << endl;


    //////////////////////////////////////////////
    // Normalise relative proprotions of 
    // different mosquito species

    double Prop_mosq_denom = 0.0;

    for (int g = 0; g < N_spec; g++)
    {
        Prop_mosq_denom = Prop_mosq_denom + Prop_mosq[g];
    }

    for (int g = 0; g < N_spec; g++)
    {
        Prop_mosq[g] = Prop_mosq[g] / Prop_mosq_denom;
    }


    ///////////////////////////////////////////////////////////
    //                                                       //
    // 1.5. Pre-multiplication of quantities for efficiency  //
    //                                                       //
    ///////////////////////////////////////////////////////////

    A_par_decay  = exp(-r_par*t_step);
    A_clin_decay = exp(-r_clin*t_step);
    mat_decay    = exp(-d_mat*t_step);

    age_0_inv = 1.0 / age_0;                 // Inverse of age-dependent biting parameter

    A_PCR_50pc_inv = CONST_LOG_2 / A_PCR_50pc;      // Immune scalar for clearance of infection
    A_LM_50pc_inv  = 1.0 / A_LM_50pc;        // Immune scalar for BS infection
    A_D_50pc_inv   = 1.0 / A_D_50pc;         // Immune scalar for clinical disease

    P_dead = 1.0 - exp(-t_step*mu_H);
    P_preg = 0.0014189;


    P_PYR_decay = exp(-PYR_decay*t_step);
    P_IRS_decay = exp(-IRS_decay*t_step);


    ///////////////////////////////////////////////////////////
    //                                                       //
    // 1.6. Fill out hypnozoite transition matrices          //
    //                                                       //
    ///////////////////////////////////////////////////////////

    for (int k1 = 0; k1<(K_max + 1); k1++)
    {
        for (int k2 = 0; k2 < (K_max + 1); k2++)
        {
            D_MAT[k1][k2] = 0.0;
            OD_MAT[k1][k2] = 0.0;
            K_MAT[k1][k2] = 0.0;
            L_MAT[k1][k2] = 0.0;
            H_MAT[k1][k2] = 0.0;
        }
    }

    for (int k = 0; k < (K_max + 1); k++)
    {
        D_MAT[k][k] = 1.0;
    }

    for (int k = 0; k < K_max; k++)
    {
        OD_MAT[k + 1][k] = 1.0;
    }
    OD_MAT[K_max][K_max] = 1.0;

    for (int k = 0; k < (K_max + 1); k++)
    {
        K_MAT[k][k] = (double)(k);
    }

    for (int k = 0; k < K_max; k++)
    {
        L_MAT[k][k + 1] = +(double)(k + 1);
        L_MAT[k + 1][k + 1] = -(double)(k + 1);
    }

    for (int k = 0; k < K_max; k++)
    {
        H_MAT[k][k] = -1.0;
        H_MAT[k + 1][k] = +1.0;
    }
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
