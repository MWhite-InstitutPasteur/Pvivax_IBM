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

#include "Simulation.hpp"

#include <iostream>
#include <fstream>
#include <cmath>


Simulation::Simulation(SimTimes times):
    times(times)
{
    /////////////////////////////////////////////////////////////////////////
    // 1.9.1. Vector of simulation times

    // Number of time steps for simulation:
    N_time = (1 / t_step)*(times.burnin + times.end - times.start) * 365;

    for (int i = 0; i<N_time; i++)
    {
        t_vec.push_back((double)(times.start * 365 - times.burnin * 365 + i*t_step));
    }


    /////////////////////////////////////////////////////////////////////////
    // 1.9.2. Create storage for output

    yH_t.resize(N_time);
    for (int i = 0; i<N_time; i++)
    {
        yH_t[i].resize(N_H_comp);
    }


    yM_t.resize(N_time);
    for (int i = 0; i<N_time; i++)
    {
        yM_t[i].resize(N_spec);
        for (int g = 0; g < N_spec; g++)
        {
            yM_t[i][g].resize(N_M_comp);
        }
    }


    prev_all.resize(N_time);
    for (int i = 0; i<N_time; i++)
    {
        prev_all[i].resize(11);
    }

    prev_U5.resize(N_time);
    for (int i = 0; i<N_time; i++)
    {
        prev_U5[i].resize(11);
    }

    prev_U10.resize(N_time);
    for (int i = 0; i<N_time; i++)
    {
        prev_U10[i].resize(11);
    }

    EIR_t.resize(N_time);

    LLIN_cov_t.resize(N_time);
    IRS_cov_t.resize(N_time);
    ACT_treat_t.resize(N_time);
    PQ_treat_t.resize(N_time);
    pregnant_t.resize(N_time);

    PQ_overtreat_t.resize(N_time);
    PQ_overtreat_9m_t.resize(N_time);

    A_par_mean_t.resize(N_time);
    A_clin_mean_t.resize(N_time);
}


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.6. Simulate the model and store the output in SIM                     //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void Simulation::run(Params& theta, Population& POP, Intervention& INTVEN)
{

    for (int i = 0; i<N_time; i++)
    {
        if (t_vec[i] / 365.0 - floor(t_vec[i] / 365.0) < 0.5*t_step / 365.0)
        {
            cout << "time = " << t_vec[i] / 365.0 << "\t" << 100.0*(t_vec[i] - t_vec[0]) / (double(t_step*N_time)) << "% complete" << endl;
        }

        human_step(theta, POP);

        mosquito_step(t_vec[i], theta, POP);

        INTVEN.distribute(t_vec[i], theta, POP);

        POP.summary();

        //////////////////////////////////////
        // Fill out Simulation object

        for (int k = 0; k<N_H_comp; k++)
        {
            yH_t[i][k] = POP.yH[k];
        }

        for (int k = 0; k<N_M_comp; k++)
        {
            for (int g = 0; g < N_spec; g++)
            {
                yM_t[i][g][k] = POP.yM[g][k];
            }
        }

        for (int k = 0; k<11; k++)
        {
            prev_all[i][k] = POP.prev_all[k];
            prev_U5[i][k] = POP.prev_U5[k];
            prev_U10[i][k] = POP.prev_U10[k];
        }


        LLIN_cov_t[i] = POP.LLIN_cov_t;
        IRS_cov_t[i] = POP.IRS_cov_t;
        ACT_treat_t[i] = POP.ACT_treat_t;
        PQ_treat_t[i] = POP.PQ_treat_t;
        pregnant_t[i] = POP.pregnant_t;

        PQ_overtreat_t[i] = POP.PQ_overtreat_t;
        PQ_overtreat_9m_t[i] = POP.PQ_overtreat_9m_t;


        EIR_t[i] = 0.0;
        for (int g = 0; g < N_spec; g++)
        {
            EIR_t[i] = EIR_t[i] + POP.aa_VC[g] * POP.yM[g][5];
        }

        A_par_mean_t[i] = POP.A_par_mean_t;
        A_clin_mean_t[i] = POP.A_clin_mean_t;

    }
}


void Simulation::write_output(const char *output_File)
{
    cout << "Start writing output to file......" << endl;
    cout << endl;

    ofstream output_Stream(output_File);

    for (int i = (int) (1/t_step)*(times.burnin)*365; i<N_time; i++)
    {
        output_Stream << t_vec[i] << "\t";

        for (int k = 0; k<N_H_comp; k++)
        {
            output_Stream << yH_t[i][k] << "\t";
        }

        for (int g = 0; g < N_spec; g++)
        {
            // Write only compartments S, E and I in mosquitoes
            for (int k = 3; k < N_M_comp; k++)
            // Write all compartments in mosquitoes
            // for (int k = 0; k < N_M_comp; k++)
            {
                output_Stream << yM_t[i][g][k] << "\t";
            }
        }

        for (int k = 0; k<10; k++)
        {
            output_Stream << prev_all[i][k] << "\t";
        }

        // Write output for age categories U5 and U10
        /*for (int k = 0; k<10; k++)
        {
            output_Stream << prev_U5[i][k] << "\t";
        }

        for (int k = 0; k<10; k++)
        {
            output_Stream << prev_U10[i][k] << "\t";
        }*/

        output_Stream << EIR_t[i] << "\t";
        output_Stream << LLIN_cov_t[i] << "\t";
        output_Stream << IRS_cov_t[i] << "\t";
        output_Stream << ACT_treat_t[i] << "\t";
        output_Stream << PQ_treat_t[i] << "\t";
        // Write number of pregnant women
        // output_Stream << pregnant_t[i] << "\t";

        output_Stream << PQ_overtreat_t[i] << "\t";
        output_Stream << PQ_overtreat_9m_t[i] << "\t";

        // Write A_par_mean_t and A_clin_mean_t
        /*output_Stream << A_par_mean_t[i] << "\t";
        output_Stream << A_clin_mean_t[i] << "\t";*/

        output_Stream << endl;
    }

    output_Stream.close();


    cout << "Output successfully written to file......" << endl;
    cout << endl;
}
