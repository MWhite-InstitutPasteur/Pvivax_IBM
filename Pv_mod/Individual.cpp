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

#include "Individual.hpp"

#include <cmath>
#include "randlib.h"


////////////////////////////////////////////////////////////
//                                                        //
//  Function declarations                                 //
//                                                        //
////////////////////////////////////////////////////////////

int CH_sample(double *xx, int nn);


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

void Individual::state_mover(Params& theta, double lam_bite)
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


    //////////////////////////////////////////////////////////////////////
    //     //                                                           // 
    //  0  //  S: Susceptible                                           //
    //     //                                                           //
    //////////////////////////////////////////////////////////////////////

    if (S == 1)
    {
        S_out = lam_H_lag;

        if (exp(-t_step*S_out) < genunf(0, 1))
        {
            //theta.r_PCR   = 1.0/( theta.d_PCR_min + (theta.d_PCR_max-theta.d_PCR_min)/( 1.0 + pow((A_par+A_par_mat)*theta.A_PCR_50pc_inv,theta.K_PCR) )); 
            theta.phi_LM = theta.phi_LM_min + (theta.phi_LM_max - theta.phi_LM_min) / (1.0 + pow((A_par + A_par_mat)*theta.A_LM_50pc_inv, theta.K_LM));
            theta.phi_D = theta.phi_D_min + (theta.phi_D_max - theta.phi_D_min) / (1.0 + pow((A_clin + A_clin_mat)*theta.A_D_50pc_inv, theta.K_D));


            S_move[0] = (1.0 - theta.phi_LM);                                  // Move to I_PCR  //  lam_H_lag*(1.0-theta.phi_LM)/S_out; 
            S_move[1] = theta.phi_LM*(1.0 - theta.phi_D);                      // Move to I_LM   //  lam_H_lag*theta.phi_LM*(1.0-theta.phi_D)/S_out;
            S_move[2] = theta.phi_LM*theta.phi_D*(1.0 - theta.treat_BScover);  // Move to I_D    //  lam_H_lag*theta.phi_LM*(1.0-theta.phi_D)*(1.0-theta.treat_cov)/S_out;
            S_move[3] = theta.phi_LM*theta.phi_D*theta.treat_BScover;          // Move to T      //  lam_H_lag*theta.phi_LM*(1.0-theta.phi_D)*theta.treat_cov/S_out;

            CH_move = CH_sample(S_move, 4);

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

        I_PCR_out = lam_H_lag + theta.r_PCR;

        if (exp(-t_step*I_PCR_out) < genunf(0, 1))
        {
            theta.phi_LM = theta.phi_LM_min + (theta.phi_LM_max - theta.phi_LM_min) / (1.0 + pow((A_par + A_par_mat)*theta.A_LM_50pc_inv, theta.K_LM));
            theta.phi_D = theta.phi_D_min + (theta.phi_D_max - theta.phi_D_min) / (1.0 + pow((A_clin + A_clin_mat)*theta.A_D_50pc_inv, theta.K_D));


            I_PCR_move[0] = theta.r_PCR / I_PCR_out;                                                       // Move to S 
            I_PCR_move[1] = lam_H_lag*(1 - theta.phi_LM) / I_PCR_out;                                      // Move to I_PCR
            I_PCR_move[2] = lam_H_lag*theta.phi_LM*(1.0 - theta.phi_D) / I_PCR_out;                        // Move to I_LM
            I_PCR_move[3] = lam_H_lag*theta.phi_LM*theta.phi_D*(1.0 - theta.treat_BScover) / I_PCR_out;        // Move to I_D
            I_PCR_move[4] = lam_H_lag*theta.phi_LM*theta.phi_D*theta.treat_BScover / I_PCR_out;                // Move to T

            CH_move = CH_sample(I_PCR_move, 5);

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
        I_LM_out = lam_H_lag + theta.r_LM;

        if (exp(-t_step*I_LM_out) < genunf(0, 1))
        {
            //theta.r_PCR = 1.0/( theta.d_PCR_min + (theta.d_PCR_max-theta.d_PCR_min)/( 1.0 + pow((A_par+A_par_mat)*theta.A_PCR_50pc_inv,theta.K_PCR) ) ); 
            theta.phi_LM = theta.phi_LM_min + (theta.phi_LM_max - theta.phi_LM_min) / (1.0 + pow((A_par + A_par_mat)*theta.A_LM_50pc_inv, theta.K_LM));
            theta.phi_D = theta.phi_D_min + (theta.phi_D_max - theta.phi_D_min) / (1.0 + pow((A_clin + A_clin_mat)*theta.A_D_50pc_inv, theta.K_D));


            I_LM_move[0] = theta.r_LM / I_LM_out;                                            // Move to I_PCR
            I_LM_move[1] = lam_H_lag*(1.0 - theta.phi_D) / I_LM_out;                         // Move to I_LM
            I_LM_move[2] = lam_H_lag*theta.phi_D*(1.0 - theta.treat_BScover) / I_LM_out;     // Move to I_D
            I_LM_move[3] = lam_H_lag*theta.phi_D*theta.treat_BScover / I_LM_out;             // Move to T

            CH_move = CH_sample(I_LM_move, 4);

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
        I_D_out = lam_H_lag + theta.r_D;

        if (exp(-t_step*I_D_out) < genunf(0, 1))
        {
            I_D_move[0] = theta.r_D / I_D_out;              // Move to I_LM
            I_D_move[1] = lam_H_lag / I_D_out;              // Move to D

            CH_move = CH_sample(I_D_move, 2);

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
        T_out = lam_H_lag + theta.r_T;

        if (exp(-t_step*T_out) < genunf(0, 1))
        {
            T_move[0] = theta.r_T / T_out;              // Move to P
            T_move[1] = lam_H_lag / T_out;              // Move to T

            CH_move = CH_sample(T_move, 2);

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
        P_out = lam_H_lag + theta.r_P;

        if (exp(-t_step*P_out) < genunf(0, 1))
        {
            P_move[0] = theta.r_P / P_out;              // Move to S
            P_move[1] = lam_H_lag / P_out;              // Move to P

            CH_move = CH_sample(P_move, 2);

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


void Individual::ager(Params& theta)
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

void Individual::intervention_updater(Params& theta)
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
