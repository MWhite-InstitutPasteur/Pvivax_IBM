/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///                                                                       ///
///  Individual-based Plasmodium vivax transmission model.                ///
///                                                                       ///
///  Dr Michael White                                                     ///
///  Institut Pasteur                                                     ///
///  michael.white@pasteur.fr                                             ///
///                                                                       ///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

#include "Population.hpp"

#include <cmath>


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
