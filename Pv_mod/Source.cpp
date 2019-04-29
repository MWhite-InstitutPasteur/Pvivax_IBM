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
///  Model code is split up into multiple files as follows:               ///
///                                                                       ///
///  -  Params.hpp, Params.cpp                                            ///
///     The Params structure stores input parameters and has              ///
///     associated code for reading parameters from input files.          ///
///                                                                       ///
///  -  Individual.hpp, Individual.cpp                                    ///
///     A class is created which stores all the information of a          ///
///     single individual.                                                ///
///     Details of the stochastic individual-based model for each         ///
///     person. Transitions occur with a fixed time step according to     ///
///     compting hazards                                                  ///
///                                                                       ///
///  -  Population.hpp, Population.cpp                                    ///
///     A structure called Population stores all individuals.             ///
///     This set of functions calculates the equilibrium set up of the    ///
///     population. It is only called once while the population is        ///
///     being initialised.                                                ///
///                                                                       ///
///  -  Intervention.hpp, Intervention.cpp                                ///
///     The Intervention struct stores intervention parameters.           ///
///                                                                       ///
///                                                                       ///
///  The code below is structured as follows:                             ///
///                                                                       ///
///  0. SETTING UP STRUCTURES AND CLASSES                                 ///
///     The time-dependent output of the model is stored in a             ///
///     structure called Simulation.                                      ///
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
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

#include "Simulation.hpp"

#include <iostream>
#include <cmath>
#include "randlib.h"


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


////////////////////////////////////////////////////////////
//                                                        //
// 0.6. Function declarations                             //
//                                                        //
////////////////////////////////////////////////////////////

void mosq_derivs(const double t, double(&yM)[N_spec][N_M_comp], double(&dyMdt)[N_spec][N_M_comp], Params& theta, Population& POP);
void mosq_rk4(const double t, const double t_step_mosq, double(&yM)[N_spec][N_M_comp], Params& theta, Population& POP);


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
    if (argc != 4 + N_spec_max)
    {
        std::cout << "Incorrect command line.\n";
        return 0;
    }

    const char* parameter_File = argv[1];

    const char* mosquito_File[N_spec_max];
    for (int g = 0; g < N_spec_max; g++)
    {
        mosquito_File[g] = argv[2 + g];
    }

    const char* coverage_File = argv[5];
    const char* output_File = argv[6];


    ////////////////////////////////////////////
    //                                        //
    // 1.2. Initialise objects                //
    //                                        //
    ////////////////////////////////////////////

    Population PNG_pop;
    Params Pv_mod_par;


    ////////////////////////////////////////////
    //                                        //
    // 1.3. Read in model parameters          //
    //                                        //
    ////////////////////////////////////////////

    SimTimes times = Pv_mod_par.read(parameter_File, mosquito_File);
    PNG_pop.N_pop = Pv_mod_par.N_pop;

    Intervention PNG_intven(coverage_File);

    ///////////////////////////////////////////////////////////////////////////
    //                                                                       //
    // 1.8. Initialise Population of individuals                             //
    //      Note that they begin with exponential age distribution           //
    //      and susceptible without immunity                                 //
    //                                                                       // 
    ///////////////////////////////////////////////////////////////////////////

    cout << "Initialise population of individuals for simulation at equilbirium EIR of " << 365.0*Pv_mod_par.EIR_equil << endl;
    cout << endl;

    equi_pop_setup(PNG_pop, Pv_mod_par);

    cout << "Population of size " << PNG_pop.N_pop << " initialised!" << endl;
    cout << endl;


    /////////////////////////////////////////////////////////////////////////
    //                                                                     //
    // 1.9. Create Simulation object                                       //
    //                                                                     //
    /////////////////////////////////////////////////////////////////////////

    Simulation PNG_sim(times);


    //////////////////////////////////////////////////////
    //                                                  //
    // 1.10. Begin stochastic simulations               //
    //                                                  //
    ////////////////////////////////////////////////////// 

    cout << "Starting model simulations......." << endl;

    PNG_sim.run(Pv_mod_par, PNG_pop, PNG_intven);

    cout << "Model simulations completed....." << endl;
    cout << endl;


    //////////////////////////////////////////////////////
    //                                                  //
    // 1.11. Output to file                             //
    //                                                  //
    ////////////////////////////////////////////////////// 

    PNG_sim.write_output(output_File);


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

void mosq_derivs(const double t, double(&yM)[N_spec][N_M_comp], double(&dyMdt)[N_spec][N_M_comp], Params& theta, Population& POP)
{
    double Karry_seas_inv[N_spec];

    for (int g = 0; g < N_spec; g++)
    {
        Karry_seas_inv[g] = 1.0 / (theta.Karry[g] * (theta.dry_seas[g] + (1 - theta.dry_seas[g])*pow(0.5 + 0.5*cos(0.01721421*(t - theta.t_peak_seas[g])), theta.kappa_seas[g]) / theta.denom_seas[g]));

        //Karry_seas_inv[g] = 1.0/theta.Karry[g];

        dyMdt[g][0] = POP.beta_VC[g] * (yM[g][3] + yM[g][4] + yM[g][5]) - yM[g][0] / theta.d_E_larvae - yM[g][0] * theta.mu_E0*(1.0 + (yM[g][0] + yM[g][1])*Karry_seas_inv[g]);
        dyMdt[g][1] = yM[g][0] / theta.d_E_larvae - yM[g][1] / theta.d_L_larvae - yM[g][1] * theta.mu_L0*(1.0 + theta.gamma_larvae*(yM[g][0] + yM[g][1])*Karry_seas_inv[g]);
        dyMdt[g][2] = yM[g][1] / theta.d_L_larvae - yM[g][2] / theta.d_pupae - yM[g][2] * theta.mu_P;
        dyMdt[g][3] = 0.5*yM[g][2] / theta.d_pupae - theta.lam_M[g] * yM[g][3] - POP.mu_M_VC[g] * yM[g][3];
        dyMdt[g][4] = +theta.lam_M[g] * yM[g][3] - theta.lam_S_M_track[g][0] * POP.exp_muM_tauM_VC[g] - POP.mu_M_VC[g] * yM[g][4];
        dyMdt[g][5] = +theta.lam_S_M_track[g][0] * POP.exp_muM_tauM_VC[g] - POP.mu_M_VC[g] * yM[g][5];
    }
}

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.2. Runge-Kutta 4 step updater for mosquito model                      //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void mosq_rk4(const double t, const double t_step_mosq, double(&yM)[N_spec][N_M_comp], Params& theta, Population& POP)
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

void mosquito_step(double t, Params& theta, Population& POP)
{
    //////////////////////////////////
    // Set up mosquito state vector

    double yM[N_spec][N_M_comp];

    for (int g = 0; g < N_spec; g++)
    {
        for (int k = 0; k<N_M_comp; k++)
        {
            yM[g][k] = POP.yM[g][k];
        }
    }


    double t_step_mosq = (double(t_step)) / (double(mosq_steps));


    //////////////////////////////////
    // Force of infection on mosquitoes

    for (int g = 0; g < N_spec; g++)
    {
        theta.lam_M[g] = 0.0;
    }

    for (int n = 0; n < POP.N_pop; n++)
    {
        for (int g = 0; g < N_spec; g++)
        {
            theta.lam_M[g] = theta.lam_M[g] + POP.lam_n[n][g] * (theta.c_PCR*POP.people[n].I_PCR + theta.c_LM*POP.people[n].I_LM +
                              theta.c_D*POP.people[n].I_D + theta.c_T*POP.people[n].T);
        }
    }


    //////////////////////////////////////
    // Carry out the mosq_steps

    for (int j = 0; j<mosq_steps; j++)
    {
        mosq_rk4(t, t_step_mosq, yM, theta, POP);

        for (int g = 0; g < N_spec; g++)
        {
            theta.lam_S_M_track[g].push_back(theta.lam_M[g] * yM[g][3]);
            theta.lam_S_M_track[g].erase(theta.lam_S_M_track[g].begin());
        }
    }

    for (int g = 0; g < N_spec; g++)
    {
        for (int k = 0; k < N_M_comp; k++)
        {
            POP.yM[g][k] = yM[g][k];
        }
    }

}
