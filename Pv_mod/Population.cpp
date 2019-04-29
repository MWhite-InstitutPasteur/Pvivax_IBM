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

#include <iostream>
#include <cmath>
#include "randlib.h"


////////////////////////////////////////////////////////////
//                                                        //
//  Function declarations                                 //
//                                                        //
////////////////////////////////////////////////////////////

void ludcmp(vector<vector<double>> &a, int n_dim, vector<int> &indx, double &d);
void lubksb(vector<vector<double>> &a, int n_dim, vector<int> &indx, vector<double> &b);
void matrix_inv(vector<vector<double>> &a, int n, vector<vector<double>> &a_inv);
void inv_MM_bb(vector<vector<double>> &MM, vector<double> &bb, vector<double> &xx, int n_dim);
void MM_ij(int i, int j, Params& theta, Population& POP, vector<vector<double>> &MM,
           vector<vector<double>> lam_eq, vector<vector<vector<double>>> phi_LM_eq,
           vector<vector<vector<double>>> phi_D_eq, vector<vector<vector<double>>> r_PCR_eq);
void gauher(Population& POP, Params& theta);


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.4. Update the vector of human classes                                 //
//                                                                          // 
//       THINK CAREFULLY ABOUT THE ORDERING OF EVENTS                       //
//////////////////////////////////////////////////////////////////////////////

void Population::human_step(Params& theta)
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

    for (int n = 0; n<N_pop; n++)
    {
        people[n].ager(theta);
    }


    ///////////////////////////////////////////////
    // 2.4.3. Deaths
    //
    // Look again at how things are erased from vectors.

    int N_dead = 0;

    for (size_t n = 0; n<people.size(); n++)
    {
        /////////////////////////////////////////////
        // Everyone has an equal probability of dying

        if (theta.P_dead > genunf(0, 1))
        {
            people.erase(people.begin() + n);

            pi_n.erase(pi_n.begin() + n);
            lam_n.erase(lam_n.begin() + n);

            N_dead = N_dead + 1;
            n = n - 1;      // If we erase something, the next one moves into it's place so we don't want to step forward.
        }
        else {

            ///////////////////////////////////////////
            // People die once they reach the maximum age

            if (people[n].age > theta.age_max)
            {
                people.erase(people.begin() + n);

                pi_n.erase(pi_n.begin() + n);
                lam_n.erase(lam_n.begin() + n);

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
        zeta_start = exp(gennor(-0.5*theta.sig_het*theta.sig_het, theta.sig_het));

        while (zeta_start > theta.het_max)
        {
            zeta_start = exp(gennor(-0.5*theta.sig_het*theta.sig_het, theta.sig_het));
        }

        Individual HH(0.0, zeta_start);

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
            if (genunf(0.0, 1.0) < theta.G6PD_prev)
            {
                HH.G6PD_def = 1;
            }
            else {
                HH.G6PD_def = 0;
            }
        }
        else {

            q_rand = genunf(0.0, 1.0);

            if (q_rand <= theta.G6PD_prev*theta.G6PD_prev)
            {
                HH.G6PD_def = 2;
            }

            if ((q_rand > theta.G6PD_prev*theta.G6PD_prev) && (q_rand <= theta.G6PD_prev*theta.G6PD_prev + 2 * theta.G6PD_prev*(1.0 - theta.G6PD_prev)))
            {
                HH.G6PD_def = 1;
            }

            if (q_rand > (theta.G6PD_prev*theta.G6PD_prev + 2 * theta.G6PD_prev*(1.0 - theta.G6PD_prev)))
            {
                HH.G6PD_def = 0;
            }
        }

        if (genunf(0.0, 1.0) < theta.CYP2D6_prev)
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

        for (size_t j = 0; j<people.size(); j++)
        {
            if (people[j].preg_age == 1)
            {
                if (abs(HH.zeta_het - people[j].zeta_het) < het_dif_track)
                {
                    HH.A_par_mat = theta.P_mat*people[j].A_par_mat;
                    HH.A_clin_mat = theta.P_mat*people[j].A_clin_mat;

                    het_dif_track = (HH.zeta_het - people[j].zeta_het)*(HH.zeta_het - people[j].zeta_het);
                }
            }
        }


        ///////////////////////////////////////////////////
        // Lagged exposure equals zero - they're not born yet!

        for (int k = 0; k<theta.H_track; k++)
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
                theta.V_int_dummy[p][q] = theta.V_int[p][q];
            }
        }

        setgmn(GMN_zero, *theta.V_int_dummy, N_int, GMN_parm);

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

        people.push_back(move(HH));

        pi_n.push_back(zero_push);
        lam_n.push_back(zero_push);
    }



    ///////////////////////////////////////////////////
    // 2.4.6. Update individual-level vector control

    for (int n = 0; n<N_pop; n++)
    {
        people[n].intervention_updater(theta);
    }


    ///////////////////////////////////////////////////
    // 2.4.7. Update proportion of bites
    //
    //        Note the ordering of n and g loops. Need to 
    //        check if this makes a difference for speed.
    //
    //        Should be able to make this quicker


    for (int n = 0; n<N_pop; n++)
    {
        for (int g = 0; g < N_spec; g++)
        {
            pi_n[n][g] = people[n].zeta_het*(1.0 - theta.rho_age*exp(-people[n].age*theta.age_0_inv));

            //pi_n[n][g] = people[n].zeta_het - (people[n].zeta_het - people[n].zeta_het)*P_age_bite;   // Slightly quicker - no calling of exponentials
        }
    }

    double SIGMA_PI[N_spec];
    for (int g = 0; g < N_spec; g++)
    {
        SIGMA_PI[g] = 0.0;
    }

    for (int n = 0; n < N_pop; n++)
    {
        for (int g = 0; g < N_spec; g++)
        {
            SIGMA_PI[g] = SIGMA_PI[g] + pi_n[n][g];
        }
    }

    for (int g = 0; g < N_spec; g++)
    {
        SIGMA_PI[g] = 1.0 / SIGMA_PI[g];
    }

    for (int n = 0; n < N_pop; n++)
    {
        for (int g = 0; g < N_spec; g++)
        {
            pi_n[n][g] = pi_n[n][g] * SIGMA_PI[g];
        }
    }


    ///////////////////////////////////////////////////
    // 2.4.8 Update population-level vector control quantities

    for (int g = 0; g < N_spec; g++)
    {
        SUM_pi_w[g] = 0;
    }

    for (int n = 0; n < N_pop; n++)
    {
        for (int g = 0; g < N_spec; g++)
        {
            SUM_pi_w[g] = SUM_pi_w[g] + pi_n[n][g] * people[n].w_VC[g];
        }
    }


    for (int g = 0; g < N_spec; g++)
    {
        W_VC[g] = 1.0 - theta.Q_0[g] + theta.Q_0[g] * SUM_pi_w[g];
        Z_VC[g] = theta.Q_0[g] * SUM_pi_z[g];

        delta_1_VC[g] = theta.delta_1 / (1.0 - Z_VC[g]);
        delta_VC[g] = delta_1_VC[g] + theta.delta_2;

        p_1_VC[g] = theta.p_1[g] * W_VC[g] / (1.0 - Z_VC[g] * theta.p_1[g]);

        mu_M_VC[g] = -log(p_1_VC[g] * theta.p_2[g]) / delta_VC[g];

        Q_VC[g] = 1.0 - (1.0 - theta.Q_0[g]) / W_VC[g];

        aa_VC[g] = Q_VC[g] / delta_VC[g];

        exp_muM_tauM_VC[g] = exp(-mu_M_VC[g] * theta.tau_M[g]);
        beta_VC[g] = theta.eps_max[g] * mu_M_VC[g] / (exp(delta_VC[g] * mu_M_VC[g]) - 1.0);
    }


    ///////////////////////////////////////////////////
    // 2.4.9. Update individual-level force of infection on humans

    for (int n = 0; n < N_pop; n++)
    {
        for (int g = 0; g < N_spec; g++)
        {
            lam_n[n][g] = aa_VC[g] * pi_n[n][g] * people[n].w_VC[g] / SUM_pi_w[g];
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
        lam_bite_base[g] = (double(N_pop))*theta.bb*yM[g][5];
    }

    for (int n = 0; n<N_pop; n++)
    {
        lam_bite_n = 0.0;

        for (int g = 0; g < N_spec; g++)
        {
            lam_bite_n = lam_bite_n + lam_n[n][g] * lam_bite_base[g];
        }

        people[n].state_mover(theta, lam_bite_n);
    }

}


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.5. Summarise the output from the population                           //
//                                                                          // 
//////////////////////////////////////////////////////////////////////////////

void Population::summary()
{
    for (int k = 0; k<N_H_comp; k++)
    {
        yH[k] = 0.0;
    }

    for (int k = 0; k<10; k++)
    {
        prev_all[k] = 0.0;
        prev_U5[k] = 0.0;
        prev_U10[k] = 0.0;
    }


    for (int n = 0; n<N_pop; n++)
    {
        ////////////////////////////////////////
        // Numbers in each compartment

        yH[0] = yH[0] + people[n].S;
        yH[1] = yH[1] + people[n].I_PCR;
        yH[2] = yH[2] + people[n].I_LM;
        yH[3] = yH[3] + people[n].I_D;
        yH[4] = yH[4] + people[n].T;
        yH[5] = yH[5] + people[n].P;


        //////////////////////////////////////////////
        //////////////////////////////////////////////
        // Summary - full population

        ////////////////////////////////////////
        // Prevalence

        prev_all[0] = prev_all[0] + 1;                                                                        // Numbers - denominator
        prev_all[1] = prev_all[1] + people[n].I_PCR + people[n].I_LM +
                                            + people[n].I_D + people[n].T;                                      // PCR detectable infections
        prev_all[2] = prev_all[2] + people[n].I_LM + people[n].I_D + people[n].T;                // LM detectable infections
        prev_all[3] = prev_all[3] + people[n].I_D + people[n].T;                                      // Clinical episodes

        if (people[n].Hyp > 0)
        {
            prev_all[4] = prev_all[4] + 1;                     // Hypnozoite positive

            prev_all[5] = prev_all[5] + people[n].Hyp;    // Number of batches of hypnozoites
        }


        ////////////////////////////////////////
        // Incidence

        prev_all[6]  = prev_all[6]  + people[n].I_PCR_new;
        prev_all[7]  = prev_all[7]  + people[n].I_LM_new;
        prev_all[8]  = prev_all[8]  + people[n].I_D_new;
        prev_all[9]  = prev_all[9]  + people[n].ACT_treat;
        prev_all[10] = prev_all[10] + people[n].PQ_treat;


        //////////////////////////////////////////////
        //////////////////////////////////////////////
        // Summary - under 5's

        if (people[n].age < 1825.0)
        {
            ////////////////////////////////////////
            // Prevalence

            prev_U5[0] = prev_U5[0] + 1;                                                                // Numbers - denominator
            prev_U5[1] = prev_U5[1] + people[n].I_PCR + people[n].I_LM
                                              + people[n].I_D + people[n].T;                              // PCR detectable infections
            prev_U5[2] = prev_U5[2] + people[n].I_LM + people[n].I_D + people[n].T;        // LM detectable infections
            prev_U5[3] = prev_U5[3] + people[n].I_D + people[n].T;                              // Clinical episodes

            if (people[n].Hyp > 0)
            {
                prev_U5[4] = prev_U5[4] + 1;                     // Hypnozoite positive

                prev_U5[5] = prev_U5[5] + people[n].Hyp;    // Number of batches of hypnozoites
            }


            ////////////////////////////////////////
            // Incidence

            prev_U5[6] = prev_U5[6] + people[n].I_PCR_new;
            prev_U5[7] = prev_U5[7] + people[n].I_LM_new;
            prev_U5[8] = prev_U5[8] + people[n].I_D_new;
            prev_U5[9] = prev_U5[9] + people[n].ACT_treat;
            prev_U5[10] = prev_U5[10] + people[n].PQ_treat;
        }

        //////////////////////////////////////////////
        //////////////////////////////////////////////
        // Summary - under 10's

        if (people[n].age < 3650.0)
        {
            ////////////////////////////////////////
            // Prevalence

            prev_U10[0] = prev_U10[0] + 1;                                                            // Numbers - denominator
            prev_U10[1] = prev_U10[1] + people[n].I_PCR + people[n].I_LM
                                                + people[n].I_D + people[n].T;                          // PCR detectable infections
            prev_U10[2] = prev_U10[2] + people[n].I_LM + people[n].I_D + people[n].T;    // LM detectable infections
            prev_U10[3] = prev_U10[3] + people[n].I_D + people[n].T;                          // Clinical episodes

            if (people[n].Hyp > 0)
            {
                prev_U10[4] = prev_U10[4] + 1;                     // Hypnozoite positive

                prev_U10[5] = prev_U10[5] + people[n].Hyp;    // Number of batches of hypnozoites
            }


            ////////////////////////////////////////
            // Incidence

            prev_U10[6]  = prev_U10[6]  + people[n].I_PCR_new;
            prev_U10[7]  = prev_U10[7]  + people[n].I_LM_new;
            prev_U10[8]  = prev_U10[8]  + people[n].I_D_new;
            prev_U10[9]  = prev_U10[9]  + people[n].ACT_treat;
            prev_U10[10] = prev_U10[10] + people[n].PQ_treat;
        }
    }


    //////////////////////////////
    // Intervention coverage

    LLIN_cov_t = 0;
    IRS_cov_t = 0;
    ACT_treat_t = 0;
    PQ_treat_t = 0;
    pregnant_t = 0;

    PQ_overtreat_t = 0;
    PQ_overtreat_9m_t = 0;


    for (int n = 0; n<N_pop; n++)
    {
        LLIN_cov_t  = LLIN_cov_t  + people[n].LLIN;
        IRS_cov_t   = IRS_cov_t   + people[n].IRS;
        ACT_treat_t = ACT_treat_t + people[n].ACT_treat;
        PQ_treat_t  = PQ_treat_t  + people[n].PQ_treat;
        pregnant_t  = pregnant_t  + people[n].pregnant;

        PQ_overtreat_t    = PQ_overtreat_t    + people[n].PQ_overtreat;
        PQ_overtreat_9m_t = PQ_overtreat_9m_t + people[n].PQ_overtreat_9m;
    }


    //////////////////////////////
    // Immunity

    double A_par_mean = 0.0, A_clin_mean = 0.0;

    for (int n = 0; n<N_pop; n++)
    {
        A_par_mean = A_par_mean + people[n].A_par;
        A_clin_mean = A_clin_mean + people[n].A_clin;
    }

    A_par_mean_t = A_par_mean / ((double)N_pop);
    A_clin_mean_t = A_clin_mean / ((double)N_pop);
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

void MM_ij(int i, int j, Params& theta, Population& POP, vector<vector<double>> &MM,
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
            MM[0 * (K_max + 1) + k1][0 * (K_max + 1) + k2] = - lam_eq[i][j]*theta.D_MAT[k1][k2] - theta.ff*theta.K_MAT[k1][k2]
                                                             + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - POP.r_age[i] * theta.D_MAT[k1][k2];
            MM[0 * (K_max + 1) + k1][1 * (K_max + 1) + k2] = + r_PCR_eq[i][j][k2]*theta.D_MAT[k1][k2];
            MM[0 * (K_max + 1) + k1][5 * (K_max + 1) + k2] = +theta.r_P*theta.D_MAT[k1][k2];

            MM[1 * (K_max + 1) + k1][0 * (K_max + 1) + k2] = + lam_eq[i][j]*(1.0 - phi_LM_eq[i][j][k2])*theta.OD_MAT[k1][k2] + theta.ff*(1.0 - phi_LM_eq[i][j][k2])*theta.K_MAT[k1][k2];
            MM[1 * (K_max + 1) + k1][1 * (K_max + 1) + k2] = - lam_eq[i][j]*theta.D_MAT[k1][k2] - theta.ff*theta.K_MAT[k1][k2] - r_PCR_eq[i][j][k2]*theta.D_MAT[k1][k2]
                                                             + lam_eq[i][j]*(1.0 - phi_LM_eq[i][j][k2])*theta.OD_MAT[k1][k2] + theta.ff*(1.0 - phi_LM_eq[i][j][k2])*theta.K_MAT[k1][k2]
                                                             + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - POP.r_age[i] * theta.D_MAT[k1][k2];
            MM[1 * (K_max + 1) + k1][2 * (K_max + 1) + k2] = + theta.r_LM*theta.D_MAT[k1][k2];

            MM[2 * (K_max + 1) + k1][0 * (K_max + 1) + k2] = + lam_eq[i][j]*phi_LM_eq[i][j][k2]*(1.0 - phi_D_eq[i][j][k2])*theta.OD_MAT[k1][k2] + theta.ff*phi_LM_eq[i][j][k2]*(1.0 - phi_D_eq[i][j][k2])*theta.K_MAT[k1][k2];
            MM[2 * (K_max + 1) + k1][1 * (K_max + 1) + k2] = + lam_eq[i][j]*phi_LM_eq[i][j][k2]*(1.0 - phi_D_eq[i][j][k2])*theta.OD_MAT[k1][k2] + theta.ff*phi_LM_eq[i][j][k2] * (1.0 - phi_D_eq[i][j][k2])*theta.K_MAT[k1][k2];
            MM[2 * (K_max + 1) + k1][2 * (K_max + 1) + k2] = - lam_eq[i][j]*theta.D_MAT[k1][k2] - theta.ff*theta.K_MAT[k1][k2] - theta.r_LM*theta.D_MAT[k1][k2]
                                                             + lam_eq[i][j]*(1.0 - phi_D_eq[i][j][k2])*theta.OD_MAT[k1][k2] + theta.ff*(1.0 - phi_D_eq[i][j][k2])*theta.K_MAT[k1][k2]
                                                             + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - POP.r_age[i] * theta.D_MAT[k1][k2];
            MM[2 * (K_max + 1) + k1][3 * (K_max + 1) + k2] = + theta.r_D*theta.D_MAT[k1][k2];

            MM[3 * (K_max + 1) + k1][0 * (K_max + 1) + k2] = + lam_eq[i][j]*phi_LM_eq[i][j][k2]*phi_D_eq[i][j][k2]*(1.0 - theta.treat_BScover*theta.treat_BSeff)*theta.OD_MAT[k1][k2] 
                                                             + theta.ff*phi_LM_eq[i][j][k2]*phi_D_eq[i][j][k2]*(1.0 - theta.treat_BScover*theta.treat_BSeff)*theta.K_MAT[k1][k2];
            MM[3 * (K_max + 1) + k1][1 * (K_max + 1) + k2] = + lam_eq[i][j]*phi_LM_eq[i][j][k2]*phi_D_eq[i][j][k2] * (1.0 - theta.treat_BScover*theta.treat_BSeff)*theta.OD_MAT[k1][k2] 
                                                             + theta.ff*phi_LM_eq[i][j][k2]*phi_D_eq[i][j][k2]*(1.0 - theta.treat_BScover*theta.treat_BSeff)*theta.K_MAT[k1][k2];
            MM[3 * (K_max + 1) + k1][2 * (K_max + 1) + k2] = + lam_eq[i][j]*phi_D_eq[i][j][k2]*(1.0 - theta.treat_BScover*theta.treat_BSeff)*theta.OD_MAT[k1][k2] + theta.ff*phi_D_eq[i][j][k2]*(1.0 - theta.treat_BScover*theta.treat_BSeff)*theta.K_MAT[k1][k2];
            MM[3 * (K_max + 1) + k1][3 * (K_max + 1) + k2] = - lam_eq[i][j]*theta.D_MAT[k1][k2] - theta.r_D*theta.D_MAT[k1][k2] + lam_eq[i][j]*theta.OD_MAT[k1][k2]
                                                             + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - POP.r_age[i]*theta.D_MAT[k1][k2];

            MM[4 * (K_max + 1) + k1][0 * (K_max + 1) + k2] = + lam_eq[i][j]*phi_LM_eq[i][j][k2]*phi_D_eq[i][j][k2]*theta.treat_BScover*theta.treat_BSeff*theta.OD_MAT[k1][k2] 
                                                             + theta.ff*phi_LM_eq[i][j][k2]*phi_D_eq[i][j][k2]*theta.treat_BScover*theta.treat_BSeff*theta.K_MAT[k1][k2];
            MM[4 * (K_max + 1) + k1][1 * (K_max + 1) + k2] = + lam_eq[i][j]*phi_LM_eq[i][j][k2]*phi_D_eq[i][j][k2]*theta.treat_BScover*theta.treat_BSeff*theta.OD_MAT[k1][k2] 
                                                             + theta.ff*phi_LM_eq[i][j][k2]*phi_D_eq[i][j][k2]*theta.treat_BScover*theta.treat_BSeff*theta.K_MAT[k1][k2];
            MM[4 * (K_max + 1) + k1][2 * (K_max + 1) + k2] = + lam_eq[i][j]*phi_D_eq[i][j][k2]*theta.treat_BScover*theta.treat_BSeff*theta.OD_MAT[k1][k2] + theta.ff*phi_D_eq[i][j][k2]*theta.treat_BScover*theta.treat_BSeff*theta.K_MAT[k1][k2];
            MM[4 * (K_max + 1) + k1][4 * (K_max + 1) + k2] = - lam_eq[i][j]*theta.D_MAT[k1][k2] - theta.r_T*theta.D_MAT[k1][k2] + lam_eq[i][j] * theta.OD_MAT[k1][k2]
                                                             + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - POP.r_age[i] * theta.D_MAT[k1][k2];

            MM[5 * (K_max + 1) + k1][4 * (K_max + 1) + k2] = + theta.r_T*theta.D_MAT[k1][k2];
            MM[5 * (K_max + 1) + k1][5 * (K_max + 1) + k2] = - lam_eq[i][j]*theta.D_MAT[k1][k2] - theta.r_P*theta.D_MAT[k1][k2] + lam_eq[i][j] * theta.OD_MAT[k1][k2]
                                                             + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - POP.r_age[i] * theta.D_MAT[k1][k2];
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

void gauher(Population& POP, Params& theta)
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
        POP.x_het[j] = exp(theta.sig_het*x[N_het - 1 - j] * sqrt(2.0) - 0.5*theta.sig_het*theta.sig_het);
        POP.w_het[j] = w[j] / w_sum;
    }


    ////////////////////////////
    // temporary for N_het = 1

    if (N_het == 1)
    {
        POP.x_het[0] = 1.0;
        POP.w_het[0] = 1.0;
    }

    for (int i = 0; i<N_age; i++)
    {
        for (int j = 0; j<N_het; j++)
        {
            POP.x_age_het[i][j] = POP.age_bite[i] * POP.x_het[j];
            POP.w_age_het[i][j] = POP.age_demog[i] * POP.w_het[j];
        }
    }

    ////////////////////////////////
    // Boundaries of heterogeneity compartments

    POP.x_het_bounds[0] = 0.0;

    for (int j = 1; j<N_het; j++)
    {
        POP.x_het_bounds[j] = exp(0.5*(log(POP.x_het[j - 1]) + log(POP.x_het[j])));
    }

    POP.x_het_bounds[N_het] = theta.het_max;
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

void equi_pop_setup(Population& POP, Params& theta)
{
    //////////////////////////////////////////////////////
    // 3.7.1. Set up age and heterogeneity compartments

    ////////////////////////////////////////
    // 3.7.1.1. Bounds of age bins

    // double age_bounds[N_age+1] = {0.0*365.0, 20.0*365.0, 40.0*365.0, 60.0*365.0, 80.0*365.0};

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
        POP.age_demog[i] = exp(-theta.mu_H*age_bounds[i]) - exp(-theta.mu_H*age_bounds[i + 1]);
    }

    POP.age_demog[N_age - 1] = 1.0;

    for (int i = 0; i<(N_age - 1); i++)
    {
        POP.age_demog[N_age - 1] = POP.age_demog[N_age - 1] - POP.age_demog[i];
    }


    ////////////////////////////////////////
    // 3.7.1.3. Ageing rates - formula below ensures
    //          balanced demography

    POP.r_age[0] = theta.mu_H*(1.0 - POP.age_demog[0]) / POP.age_demog[0];

    for (int i = 1; i<(N_age - 1); i++)
    {
        POP.r_age[i] = (POP.r_age[i - 1] * POP.age_demog[i - 1] - theta.mu_H*POP.age_demog[i]) / POP.age_demog[i];
    }

    POP.r_age[N_age - 1] = 0.0;


    ////////////////////////////////////////
    // 3.7.1.4. Age-dependent mosquito biting rates

    for (int i = 0; i<N_age; i++)
    {
        POP.age_mids[i] = 0.5*(age_bounds[i] + age_bounds[i + 1]);
    }

    for (int i = 0; i<N_age; i++)
    {
        POP.age_bite[i] = 1.0 - theta.rho_age*exp(-POP.age_mids[i] / theta.age_0);
    }

    POP.P_age_bite = exp(-t_step / theta.age_0);


    ///////////////////////////////////////
    // Ensure total bites are normalised

    double omega_age = 0.0;

    for (int i = 0; i<N_age; i++)
    {
        omega_age = omega_age + POP.age_demog[i] * POP.age_bite[i];
    }

    omega_age = 1 / omega_age;


    for (int i = 0; i<N_age; i++)
    {
        POP.age_bite[i] = omega_age*POP.age_bite[i];
    }


    //////////////////////////////////////////////////////
    // 3.7.1.5. Find age category closest to 20 yr old woman

    POP.index_age_20 = 0;

    double age_diff = (POP.age_mids[0] - 20.0*365.0)*(POP.age_mids[0] - 20.0*365.0);

    for (int i = 1; i<N_age; i++)
    {
        if ((POP.age_mids[i] - 20.0*365.0)*(POP.age_mids[i] - 20.0*365.0) < age_diff)
        {
            age_diff = (POP.age_mids[i] - 20.0*365.0)*(POP.age_mids[i] - 20.0*365.0);
            POP.index_age_20 = i;
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
            lam_eq[i][j] = theta.EIR_equil*theta.bb*POP.x_age_het[i][j];
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
        HH_bb[0] = POP.w_het[j] * theta.mu_H;


        for (int k1 = 0; k1 < (K_max + 1); k1++)
        {
            for (int k2 = 0; k2 < (K_max + 1); k2++)
            {
                HH_mat[k1][k2] = lam_eq[0][j] * theta.H_MAT[k1][k2] + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - POP.r_age[0] * theta.D_MAT[k1][k2];
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
                HH_bb[k] = -POP.r_age[i - 1] * HH_xx[k];
            }

            for (int k1 = 0; k1 < (K_max + 1); k1++)
            {
                for (int k2 = 0; k2 < (K_max + 1); k2++)
                {
                    HH_mat[k1][k2] = lam_eq[i][j] * theta.H_MAT[k1][k2] + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - POP.r_age[i] * theta.D_MAT[k1][k2];
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
            G_VEC[k] = lam_eq[0][j] * HH_eq[0][j][k - 1] / HH_eq[0][j][k] + theta.ff*((double)k);
        }
        G_VEC[K_max] = G_VEC[K_max] + lam_eq[0][j];

        for (int k = 0; k < (K_max + 1); k++)
        {
            G_VEC[k] = G_VEC[k] / (G_VEC[k] * theta.u_par + 1.0);
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
                ODE_eq_MAT[k1][k2] = -(lam_eq[0][j] * LAM_MAT[k1][k2] + theta.gamma_L*GAM_MAT[k1][k2] -
                    theta.r_par*theta.D_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - POP.r_age[0] * theta.D_MAT[k1][k2]);
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
                G_VEC[k] = lam_eq[i][j] * HH_eq[i][j][k - 1] / HH_eq[i][j][k] + theta.ff*((double)k);
            }
            G_VEC[K_max] = G_VEC[K_max] + lam_eq[i][j];

            for (int k = 0; k < (K_max + 1); k++)
            {
                G_VEC[k] = G_VEC[k] / (G_VEC[k] * theta.u_par + 1.0);
            }


            for (int k = 0; k < (K_max + 1); k++)
            {
                ODE_eq_vec[k] = G_VEC[k] + POP.r_age[i - 1] * A_par_eq[i - 1][j][k] * HH_eq[i - 1][j][k] / HH_eq[i][j][k];
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
                    ODE_eq_MAT[k1][k2] = -(lam_eq[i][j] * LAM_MAT[k1][k2] + theta.gamma_L*GAM_MAT[k1][k2] -
                        theta.r_par*theta.D_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - POP.r_age[i] * theta.D_MAT[k1][k2]);
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
            G_VEC[k] = lam_eq[0][j] * (HH_eq[0][j][k - 1] / HH_eq[0][j][k]) + theta.ff*((double)k);
        }
        G_VEC[K_max] = G_VEC[K_max] + lam_eq[0][j];

        for (int k = 0; k < (K_max + 1); k++)
        {
            G_VEC[k] = G_VEC[k] / (G_VEC[k] * theta.u_clin + 1.0);
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
                ODE_eq_MAT[k1][k2] = -(lam_eq[0][j] * LAM_MAT[k1][k2] + theta.gamma_L*GAM_MAT[k1][k2] -
                    theta.r_clin*theta.D_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - POP.r_age[0] * theta.D_MAT[k1][k2]);
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
                G_VEC[k] = lam_eq[i][j] * HH_eq[i][j][k - 1] / HH_eq[i][j][k] + theta.ff*((double)k);
            }
            G_VEC[K_max] = G_VEC[K_max] + lam_eq[i][j];

            for (int k = 0; k < (K_max + 1); k++)
            {
                G_VEC[k] = G_VEC[k] / (G_VEC[k] * theta.u_clin + 1.0);
            }


            for (int k = 0; k < (K_max + 1); k++)
            {
                ODE_eq_vec[k] = G_VEC[k] + POP.r_age[i - 1] * A_clin_eq[i - 1][j][k] * HH_eq[i - 1][j][k] / HH_eq[i][j][k];
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
                    ODE_eq_MAT[k1][k2] = -(lam_eq[i][j] * LAM_MAT[k1][k2] + theta.gamma_L*GAM_MAT[k1][k2] -
                        theta.r_clin*theta.D_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - POP.r_age[i] * theta.D_MAT[k1][k2]);
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
                A_par_eq[i][j][k] = A_par_eq[i][j][k] + A_par_eq_mean[POP.index_age_20][j] * theta.P_mat*exp(-POP.age_mids[i] / theta.d_mat);
                A_clin_eq[i][j][k] = A_clin_eq[i][j][k] + A_clin_eq_mean[POP.index_age_20][j] * theta.P_mat*exp(-POP.age_mids[i] / theta.d_mat);
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
                phi_LM_eq[i][j][k] = theta.phi_LM_min + (theta.phi_LM_max - theta.phi_LM_min) / (1.0 + pow(A_par_eq[i][j][k] / theta.A_LM_50pc, theta.K_LM));
                phi_D_eq[i][j][k] = theta.phi_D_min + (theta.phi_D_max - theta.phi_D_min) / (1.0 + pow(A_clin_eq[i][j][k] / theta.A_D_50pc, theta.K_D));
                r_PCR_eq[i][j][k] = 1.0 / (theta.d_PCR_min + (theta.d_PCR_max - theta.d_PCR_min) / (1.0 + pow(A_par_eq[i][j][k] / theta.A_PCR_50pc, theta.K_PCR)));
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

        bb[0] = -POP.w_het[j] * theta.mu_H;

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
                bb[c] = -POP.r_age[i - 1] * xx[c];
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
        theta.lam_M[g] = 0.0;

        for (int i = 0; i<N_age; i++)
        {
            for (int j = 0; j<N_het; j++)
            {
                for (int k = 0; k < (K_max + 1); k++)
                {
                    theta.lam_M[g] = theta.lam_M[g] + POP.x_age_het[i][j] * theta.aa[g] * (theta.c_PCR*yH_eq[i][j][k][1] + theta.c_LM*yH_eq[i][j][k][2] +
                        theta.c_D*yH_eq[i][j][k][3] + theta.c_T*yH_eq[i][j][k][4]);
                }
            }
        }
    }


    double I_M_star[N_spec];
    for (int g = 0; g < N_spec; g++)
    {
        I_M_star[g] = theta.lam_M[g] * exp(-theta.mu_M[g] * theta.tau_M[g]) / (theta.lam_M[g] + theta.mu_M[g]);
    }

    double a_I_M_sum;

    if (theta.Prop_mosq[0] > 0.0)
    {
        a_I_M_sum = theta.aa[0] * I_M_star[0];

        for (int g = 1; g < N_spec; g++)
        {
            a_I_M_sum = a_I_M_sum + (theta.Prop_mosq[g] / theta.Prop_mosq[0])*theta.aa[g] * I_M_star[g];
        }

        theta.mm_0[0] = theta.EIR_equil / a_I_M_sum;

        for (int g = 1; g < N_spec; g++)
        {
            theta.mm_0[g] = (theta.Prop_mosq[g] / theta.Prop_mosq[0])*theta.mm_0[0];
        }
    }

    if ((theta.Prop_mosq[0] < 1.0e-10) && (theta.Prop_mosq[1] > 0.0))
    {
        theta.mm_0[0] = 0.0;

        a_I_M_sum = theta.aa[1] * I_M_star[1];

        for (int g = 2; g < N_spec; g++)
        {
            a_I_M_sum = a_I_M_sum + (theta.Prop_mosq[g] / theta.Prop_mosq[1])*theta.aa[g] * I_M_star[g];
        }

        theta.mm_0[1] = theta.EIR_equil / a_I_M_sum;

        for (int g = 2; g < N_spec; g++)
        {
            theta.mm_0[g] = (theta.Prop_mosq[g] / theta.Prop_mosq[1])*theta.mm_0[1];
        }
    }


    if ((theta.Prop_mosq[0] < 1.0e-10) && (theta.Prop_mosq[1] < 1.0e-10))
    {
        theta.mm_0[0] = 0.0;
        theta.mm_0[1] = 0.0;
        theta.mm_0[2] = theta.EIR_equil / (theta.aa[2] * I_M_star[2]);

        //theta.mm_0[2] = theta.EIR_equil / (theta.aa[2]*(theta.lam_M[2] / (theta.lam_M[2] + theta.mu_M[2]))*exp(-theta.mu_M[2]*theta.tau_M[2]));
    }



    for (int g = 0; g < N_spec; g++)
    {
        POP.yM[g][0] = 2.0*theta.omega_larvae[g] * theta.mu_M[g] * theta.d_L_larvae*(1.0 + theta.d_pupae*theta.mu_P)*theta.mm_0[g];
        POP.yM[g][1] = 2.0*theta.mu_M[g] * theta.d_L_larvae*(1.0 + theta.d_pupae*theta.mu_P)*theta.mm_0[g];
        POP.yM[g][2] = 2.0*theta.d_pupae*theta.mu_M[g] * theta.mm_0[g];
        POP.yM[g][3] = theta.mm_0[g] * (theta.mu_M[g] / (theta.lam_M[g] + theta.mu_M[g]));
        POP.yM[g][4] = theta.mm_0[g] * (theta.lam_M[g] / (theta.lam_M[g] + theta.mu_M[g]))*(1.0 - exp(-theta.mu_M[g] * theta.tau_M[g]));
        POP.yM[g][5] = theta.mm_0[g] * (theta.lam_M[g] / (theta.lam_M[g] + theta.mu_M[g]))*exp(-theta.mu_M[g] * theta.tau_M[g]);

        theta.Karry[g] = theta.mm_0[g] * 2.0*theta.d_L_larvae*theta.mu_M[g] * (1.0 + theta.d_pupae*theta.mu_P)*theta.gamma_larvae*(theta.omega_larvae[g] + 1.0) /
            (theta.omega_larvae[g] / (theta.mu_L0*theta.d_E_larvae) - 1.0 / (theta.mu_L0*theta.d_L_larvae) - 1.0);   // Larval carry capacity


        if (theta.Karry[g] < 1.0e-10) { theta.Karry[g] = 1.0e-10; } //
    }



    for (int g = 0; g < N_spec; g++)
    {
        if (g == 0) { cout << "An. farauti:  " << 100.0 * theta.Prop_mosq[0] << "%" << endl; }
        if (g == 1) { cout << "An. punctulatus:  " << 100.0 * theta.Prop_mosq[1] << "%" << endl; }
        if (g == 2) { cout << "An. koliensis:  " << 100.0 * theta.Prop_mosq[2] << "%" << endl; }

        cout << "EL_M  " << POP.yM[g][0] << endl;
        cout << "LL_M  " << POP.yM[g][1] << endl;
        cout << "P_M  " << POP.yM[g][2] << endl;
        cout << "S_M  " << POP.yM[g][3] << endl;
        cout << "E_M  " << POP.yM[g][4] << endl;
        cout << "I_M  " << POP.yM[g][5] << endl;

        cout << "lam_M = " << theta.lam_M[g] << endl;

        cout << "I_M = " << POP.yM[g][5] << endl;

        cout << "mm = " << theta.mm_0[g] << endl;

        cout << endl;
    }
    cout << endl;

    cout << "lam_H = " << theta.bb*theta.EIR_equil << endl;

    double EIR_out = 0.0;
    for (int g = 0; g < N_spec; g++)
    {
        EIR_out = EIR_out + 365.0*theta.aa[g] * POP.yM[g][5];
    }

    cout << "EIR = " << EIR_out << endl;


    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    // 3.7.4. Proportion in each age and heterogeneity stratified category

    //////////////////////////////////////////
    // Fill out vector of lagged lam_M*S_M

    theta.lam_S_M_track.resize(N_spec);

    for (int g = 0; g<N_spec; g++)
    {
        for (int k = 0; k < theta.M_track; k++)
        {
            theta.lam_S_M_track[g].push_back(theta.lam_M[g] * POP.yM[g][3]);
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

    for (int n = 0; n<POP.N_pop; n++)
    {
        //////////////////////////////////////////////////////////////////
        // 3.7.4.2.1. Assign age and heterogeneity 

        age_start = genexp(theta.age_mean);

        while (age_start > theta.age_max)
        {
            age_start = genexp(theta.age_mean);
        }

        zeta_start = exp(gennor(-0.5*theta.sig_het*theta.sig_het, theta.sig_het));

        while (zeta_start > theta.het_max)
        {
            zeta_start = exp(gennor(-0.5*theta.sig_het*theta.sig_het, theta.sig_het));
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
            if ((zeta_start > POP.x_het_bounds[j]) && (zeta_start < POP.x_het_bounds[j + 1]))
            {
                j_index = j;
            }
        }

        //////////////////////////////////////////////////////////////////
        // 3.7.4.2.3. Construct a new individual

        double q_rand;

        Individual HH(age_start, zeta_start);

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
            if (genunf(0.0, 1.0) < theta.G6PD_prev)
            {
                HH.G6PD_def = 1;
            }
            else {
                HH.G6PD_def = 0;
            }
        }
        else {
            q_rand = genunf(0.0, 1.0);

            if (q_rand <= theta.G6PD_prev*theta.G6PD_prev)
            {
                HH.G6PD_def = 2;
            }

            if ((q_rand > theta.G6PD_prev*theta.G6PD_prev) && (q_rand <= theta.G6PD_prev*theta.G6PD_prev + 2 * theta.G6PD_prev*(1.0 - theta.G6PD_prev)))
            {
                HH.G6PD_def = 1;
            }

            if (q_rand >  theta.G6PD_prev*theta.G6PD_prev + 2 * theta.G6PD_prev*(1.0 - theta.G6PD_prev))
            {
                HH.G6PD_def = 0;
            }
        }


        if (genunf(0.0, 1.0) < theta.CYP2D6_prev)
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

        HH.A_par_mat = A_par_eq_mean[POP.index_age_20][j_index] * theta.P_mat*exp(-POP.age_mids[i_index] / theta.d_mat);
        HH.A_clin_mat = A_clin_eq_mean[POP.index_age_20][j_index] * theta.P_mat*exp(-POP.age_mids[i_index] / theta.d_mat);

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

        for (int k = 0; k<theta.H_track; k++)
        {
            HH.lam_bite_track.push_back(lam_eq[i_index][j_index]);
        }

        for (int k = 0; k<theta.H_track; k++)
        {
            HH.lam_rel_track.push_back(HH.Hyp*theta.ff);
        }


        ////////////////////////////////////////////////////////
        // 3.7.4.2.7. Give individuals their life-long intervention access score

        for (int p = 0; p<N_int; p++)
        {
            for (int q = 0; q<N_int; q++)
            {
                theta.V_int_dummy[p][q] = theta.V_int[p][q];
            }
        }

        setgmn(GMN_zero, *theta.V_int_dummy, N_int, GMN_parm);

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

        POP.people.push_back(move(HH));
    }

    ///////////////////////////////////////////////////////////
    // 3.7.5.1. Proportion of bites received by each person

    POP.pi_n.resize(POP.N_pop);
    for (int n = 0; n < POP.N_pop; n++)
    {
        POP.pi_n[n].resize(N_spec);
    }

    for (int n = 0; n<POP.N_pop; n++)
    {
        for (int g = 0; g < N_spec; g++)
        {
            POP.pi_n[n][g] = POP.people[n].zeta_het*(1.0 - theta.rho_age*exp(-POP.people[n].age*theta.age_0_inv));
        }
    }



    double SIGMA_PI[N_spec];

    for (int g = 0; g < N_spec; g++)
    {
        SIGMA_PI[g] = 0.0;
        for (int n = 0; n<POP.N_pop; n++)
        {
            SIGMA_PI[g] = SIGMA_PI[g] + POP.pi_n[n][g];
        }

        for (int n = 0; n<POP.N_pop; n++)
        {
            POP.pi_n[n][g] = POP.pi_n[n][g] / SIGMA_PI[g];
        }
    }

    ///////////////////////////////////////////////////////////
    // 3.7.5. Initialise population level quantitites

    //////////////////////////////////////////////
    // 3.7.5.1. Vector control quantities

    for (int g = 0; g < N_spec; g++)
    {
        POP.SUM_pi_w[g] = 0;
        POP.SUM_pi_z[g] = 0;


        for (int n = 0; n < POP.N_pop; n++)
        {
            POP.SUM_pi_w[g] = POP.SUM_pi_w[g] + POP.pi_n[n][g] * POP.people[n].w_VC[g];
            POP.SUM_pi_z[g] = POP.SUM_pi_z[g] + POP.pi_n[n][g] * POP.people[n].z_VC[g];
        }
    }


    for (int g = 0; g < N_spec; g++)
    {
        POP.W_VC[g] = 1.0 - theta.Q_0[g] + theta.Q_0[g] * POP.SUM_pi_w[g];
        POP.Z_VC[g] = theta.Q_0[g] * POP.SUM_pi_z[g];

        POP.delta_1_VC[g] = theta.delta_1 / (1.0 - POP.Z_VC[g]);
        POP.delta_VC[g] = POP.delta_1_VC[g] + theta.delta_2;

        POP.p_1_VC[g] = exp(-theta.mu_M[g] * POP.delta_1_VC[g]);
        POP.mu_M_VC[g] = -log(POP.p_1_VC[g] * theta.p_2[g]) / POP.delta_VC[g];

        POP.Q_VC[g] = 1.0 - (1.0 - theta.Q_0[g]) / POP.W_VC[g];

        POP.aa_VC[g] = POP.Q_VC[g] / POP.delta_VC[g];

        POP.exp_muM_tauM_VC[g] = exp(-POP.mu_M_VC[g] * theta.tau_M[g]);
        POP.beta_VC[g] = theta.eps_max[g] * POP.mu_M_VC[g] / (exp(POP.delta_VC[g] * POP.mu_M_VC[g]) - 1.0);
    }


    /////////////////////////////////////////////////////////////////////////
    // 3.7.5.2. The rate at which person n is bitten by a single mosquito

    POP.lam_n.resize(POP.N_pop);
    for (int n = 0; n < POP.N_pop; n++)
    {
        POP.lam_n[n].resize(N_spec);
    }


    for (int n = 0; n < POP.N_pop; n++)
    {
        for (int g = 0; g < N_spec; g++)
        {
            POP.lam_n[n][g] = theta.aa[g] * POP.pi_n[n][g];
        }
    }


    /////////////////////////////////////////////////////////////////////////
    // 3.7.5.3. Output an overview of the initial set up

    cout << "Equilibrium set up......." << endl;

    double S_ind = 0.0, I_PCR_ind = 0.0, I_LM_ind = 0.0, I_D_ind = 0.0, T_ind = 0.0, P_ind = 0.0;
    double S_eqq = 0.0, I_PCR_eqq = 0.0, I_LM_eqq = 0.0, I_D_eqq = 0.0, T_eqq = 0.0, P_eqq = 0.0;

    for (int n = 0; n<POP.N_pop; n++)
    {
        S_ind = S_ind + POP.people[n].S;
        I_PCR_ind = I_PCR_ind + POP.people[n].I_PCR;
        I_LM_ind = I_LM_ind + POP.people[n].I_LM;
        I_D_ind = I_D_ind + POP.people[n].I_D;
        T_ind = T_ind + POP.people[n].T;
        P_ind = P_ind + POP.people[n].P;
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


    cout << "S = " << ((double)S_ind) / POP.N_pop << "\t" << S_eqq << endl;
    cout << "I_PCR = " << ((double)I_PCR_ind) / POP.N_pop << "\t" << I_PCR_eqq << endl;
    cout << "I_LM = " << ((double)I_LM_ind) / POP.N_pop << "\t" << I_LM_eqq << endl;
    cout << "I_D = " << ((double)I_D_ind) / POP.N_pop << "\t" << I_D_eqq << endl;
    cout << "T = " << ((double)T_ind) / POP.N_pop << "\t" << T_eqq << endl;
    cout << "P = " << ((double)P_ind) / POP.N_pop << "\t" << P_eqq << endl;

}
