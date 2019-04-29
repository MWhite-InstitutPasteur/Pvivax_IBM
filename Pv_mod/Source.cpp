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

#include "Population.hpp"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <time.h>
#include "randlib.h"
#include <omp.h>
#include <algorithm>
#include <regex>


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


/////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
// 0.4. Define a structure for details of interventions                                //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

struct Intervention
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

struct Simulation
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

void mosq_derivs(const double t, double(&yM)[N_spec][N_M_comp], double(&dyMdt)[N_spec][N_M_comp], Params& theta, Population& POP);
void mosq_rk4(const double t, const double t_step_mosq, double(&yM)[N_spec][N_M_comp], Params& theta, Population& POP);
void mosquito_step(double t, Params& theta, Population& POP);
void human_step(Params& theta, Population& POP);
void intervention_dist(double t, Params& theta, Population& POP, Intervention& INTVEN);
void POP_summary(Population& POP);
void model_simulator(Params& theta, Population& POP, Intervention& INTVEN, Simulation& SIM);
double phi_inv(double pp, double mu, double sigma);


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

    Pv_mod_par.read(parameter_File, mosquito_File);
    PNG_pop.N_pop = Pv_mod_par.N_pop;

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
    //        There's very likely a much more effective way to do this.

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

    string discard;

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
    // 1.7.2. Fill out Intervention structure

    Intervention PNG_intven;

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

    equi_pop_setup(PNG_pop, Pv_mod_par);

    cout << "Population of size " << PNG_pop.N_pop << " initialised!" << endl;
    cout << endl;


    /////////////////////////////////////////////////////////////////////////
    //                                                                     //
    // 1.9. Create Simulation object                                       //
    //                                                                     //
    /////////////////////////////////////////////////////////////////////////

    Simulation PNG_sim;


    /////////////////////////////////////////////////////////////////////////
    // 1.9.1. Vector of simulation times

    // Number of time steps for simulation:
    int N_time = (1 / t_step)*(Pv_mod_par.burnin_time + Pv_mod_par.time_end - Pv_mod_par.time_start) * 365;

    PNG_sim.N_time = N_time;

    for (int i = 0; i<N_time; i++)
    {
        PNG_sim.t_vec.push_back((double)(Pv_mod_par.time_start * 365 - Pv_mod_par.burnin_time * 365 + i*t_step));
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

    model_simulator(Pv_mod_par, PNG_pop, PNG_intven, PNG_sim);

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

    for (int i = (int) (1/t_step)*(Pv_mod_par.burnin_time)*365; i<N_time; i++)
    {
        output_Stream << PNG_sim.t_vec[i] << "\t";

        for (int k = 0; k<N_H_comp; k++)
        {
            output_Stream << PNG_sim.yH_t[i][k] << "\t";
        }

        for (int g = 0; g < N_spec; g++)
        {
            // Write only compartments S, E and I in mosquitoes
            for (int k = 3; k < N_M_comp; k++)
            // Write all compartments in mosquitoes
            // for (int k = 0; k < N_M_comp; k++)
            {
                output_Stream << PNG_sim.yM_t[i][g][k] << "\t";
            }
        }

        for (int k = 0; k<10; k++)
        {
            output_Stream << PNG_sim.prev_all[i][k] << "\t";
        }

        // Write output for age categories U5 and U10
        /*for (int k = 0; k<10; k++)
        {
            output_Stream << PNG_sim.prev_U5[i][k] << "\t";
        }

        for (int k = 0; k<10; k++)
        {
            output_Stream << PNG_sim.prev_U10[i][k] << "\t";
        }*/

        output_Stream << PNG_sim.EIR_t[i] << "\t";
        output_Stream << PNG_sim.LLIN_cov_t[i] << "\t";
        output_Stream << PNG_sim.IRS_cov_t[i] << "\t";
        output_Stream << PNG_sim.ACT_treat_t[i] << "\t";
        output_Stream << PNG_sim.PQ_treat_t[i] << "\t";
        // Write number of pregnant women
        // output_Stream << PNG_sim.pregnant_t[i] << "\t";

        output_Stream << PNG_sim.PQ_overtreat_t[i] << "\t";
        output_Stream << PNG_sim.PQ_overtreat_9m_t[i] << "\t";

        // Write A_par_mean_t and A_clin_mean_t
        /*output_Stream << PNG_sim.A_par_mean_t[i] << "\t";
        output_Stream << PNG_sim.A_clin_mean_t[i] << "\t";*/

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


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.4. Update the vector of human classes                                 //
//                                                                          // 
//       THINK CAREFULLY ABOUT THE ORDERING OF EVENTS                       //
//////////////////////////////////////////////////////////////////////////////

void human_step(Params& theta, Population& POP)
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

    for (int n = 0; n<POP.N_pop; n++)
    {
        POP.people[n].ager(theta);
    }


    ///////////////////////////////////////////////
    // 2.4.3. Deaths
    //
    // Look again at how things are erased from vectors.

    int N_dead = 0;

    for (size_t n = 0; n<POP.people.size(); n++)
    {
        /////////////////////////////////////////////
        // Everyone has an equal probability of dying

        if (theta.P_dead > genunf(0, 1))
        {
            POP.people.erase(POP.people.begin() + n);

            POP.pi_n.erase(POP.pi_n.begin() + n);
            POP.lam_n.erase(POP.lam_n.begin() + n);

            N_dead = N_dead + 1;
            n = n - 1;      // If we erase something, the next one moves into it's place so we don't want to step forward.
        }
        else {

            ///////////////////////////////////////////
            // People die once they reach the maximum age

            if (POP.people[n].age > theta.age_max)
            {
                POP.people.erase(POP.people.begin() + n);

                POP.pi_n.erase(POP.pi_n.begin() + n);
                POP.lam_n.erase(POP.lam_n.begin() + n);

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

        for (size_t j = 0; j<POP.people.size(); j++)
        {
            if (POP.people[j].preg_age == 1)
            {
                if (abs(HH.zeta_het - POP.people[j].zeta_het) < het_dif_track)
                {
                    HH.A_par_mat = theta.P_mat*POP.people[j].A_par_mat;
                    HH.A_clin_mat = theta.P_mat*POP.people[j].A_clin_mat;

                    het_dif_track = (HH.zeta_het - POP.people[j].zeta_het)*(HH.zeta_het - POP.people[j].zeta_het);
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

        POP.people.push_back(HH);

        POP.pi_n.push_back(zero_push);
        POP.lam_n.push_back(zero_push);
    }



    ///////////////////////////////////////////////////
    // 2.4.6. Update individual-level vector control

    for (int n = 0; n<POP.N_pop; n++)
    {
        POP.people[n].intervention_updater(theta);
    }


    ///////////////////////////////////////////////////
    // 2.4.7. Update proportion of bites
    //
    //        Note the ordering of n and g loops. Need to 
    //        check if this makes a difference for speed.
    //
    //        Should be able to make this quicker


    for (int n = 0; n<POP.N_pop; n++)
    {
        for (int g = 0; g < N_spec; g++)
        {
            POP.pi_n[n][g] = POP.people[n].zeta_het*(1.0 - theta.rho_age*exp(-POP.people[n].age*theta.age_0_inv));

            //POP.pi_n[n][g] = POP.people[n].zeta_het - (POP.people[n].zeta_het - POP.people[n].zeta_het)*POP.P_age_bite;   // Slightly quicker - no calling of exponentials
        }
    }

    double SIGMA_PI[N_spec];
    for (int g = 0; g < N_spec; g++)
    {
        SIGMA_PI[g] = 0.0;
    }

    for (int n = 0; n < POP.N_pop; n++)
    {
        for (int g = 0; g < N_spec; g++)
        {
            SIGMA_PI[g] = SIGMA_PI[g] + POP.pi_n[n][g];
        }
    }

    for (int g = 0; g < N_spec; g++)
    {
        SIGMA_PI[g] = 1.0 / SIGMA_PI[g];
    }

    for (int n = 0; n < POP.N_pop; n++)
    {
        for (int g = 0; g < N_spec; g++)
        {
            POP.pi_n[n][g] = POP.pi_n[n][g] * SIGMA_PI[g];
        }
    }


    ///////////////////////////////////////////////////
    // 2.4.8 Update population-level vector control quantities

    for (int g = 0; g < N_spec; g++)
    {
        POP.SUM_pi_w[g] = 0;
    }

    for (int n = 0; n < POP.N_pop; n++)
    {
        for (int g = 0; g < N_spec; g++)
        {
            POP.SUM_pi_w[g] = POP.SUM_pi_w[g] + POP.pi_n[n][g] * POP.people[n].w_VC[g];
        }
    }


    for (int g = 0; g < N_spec; g++)
    {
        POP.W_VC[g] = 1.0 - theta.Q_0[g] + theta.Q_0[g] * POP.SUM_pi_w[g];
        POP.Z_VC[g] = theta.Q_0[g] * POP.SUM_pi_z[g];

        POP.delta_1_VC[g] = theta.delta_1 / (1.0 - POP.Z_VC[g]);
        POP.delta_VC[g] = POP.delta_1_VC[g] + theta.delta_2;

        POP.p_1_VC[g] = theta.p_1[g] * POP.W_VC[g] / (1.0 - POP.Z_VC[g] * theta.p_1[g]);

        POP.mu_M_VC[g] = -log(POP.p_1_VC[g] * theta.p_2[g]) / POP.delta_VC[g];

        POP.Q_VC[g] = 1.0 - (1.0 - theta.Q_0[g]) / POP.W_VC[g];

        POP.aa_VC[g] = POP.Q_VC[g] / POP.delta_VC[g];

        POP.exp_muM_tauM_VC[g] = exp(-POP.mu_M_VC[g] * theta.tau_M[g]);
        POP.beta_VC[g] = theta.eps_max[g] * POP.mu_M_VC[g] / (exp(POP.delta_VC[g] * POP.mu_M_VC[g]) - 1.0);
    }


    ///////////////////////////////////////////////////
    // 2.4.9. Update individual-level force of infection on humans

    for (int n = 0; n < POP.N_pop; n++)
    {
        for (int g = 0; g < N_spec; g++)
        {
            POP.lam_n[n][g] = POP.aa_VC[g] * POP.pi_n[n][g] * POP.people[n].w_VC[g] / POP.SUM_pi_w[g];
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
        lam_bite_base[g] = (double(POP.N_pop))*theta.bb*POP.yM[g][5];
    }

    for (int n = 0; n<POP.N_pop; n++)
    {
        lam_bite_n = 0.0;

        for (int g = 0; g < N_spec; g++)
        {
            lam_bite_n = lam_bite_n + POP.lam_n[n][g] * lam_bite_base[g];
        }

        POP.people[n].state_mover(theta, lam_bite_n);
    }

}


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.5. Summarise the output from the population                           //
//                                                                          // 
//////////////////////////////////////////////////////////////////////////////

void POP_summary(Population& POP)
{
    for (int k = 0; k<N_H_comp; k++)
    {
        POP.yH[k] = 0.0;
    }

    for (int k = 0; k<10; k++)
    {
        POP.prev_all[k] = 0.0;
        POP.prev_U5[k] = 0.0;
        POP.prev_U10[k] = 0.0;
    }


    for (int n = 0; n<POP.N_pop; n++)
    {
        ////////////////////////////////////////
        // Numbers in each compartment

        POP.yH[0] = POP.yH[0] + POP.people[n].S;
        POP.yH[1] = POP.yH[1] + POP.people[n].I_PCR;
        POP.yH[2] = POP.yH[2] + POP.people[n].I_LM;
        POP.yH[3] = POP.yH[3] + POP.people[n].I_D;
        POP.yH[4] = POP.yH[4] + POP.people[n].T;
        POP.yH[5] = POP.yH[5] + POP.people[n].P;


        //////////////////////////////////////////////
        //////////////////////////////////////////////
        // Summary - full population

        ////////////////////////////////////////
        // Prevalence

        POP.prev_all[0] = POP.prev_all[0] + 1;                                                                        // Numbers - denominator
        POP.prev_all[1] = POP.prev_all[1] + POP.people[n].I_PCR + POP.people[n].I_LM +
                                            + POP.people[n].I_D + POP.people[n].T;                                      // PCR detectable infections
        POP.prev_all[2] = POP.prev_all[2] + POP.people[n].I_LM + POP.people[n].I_D + POP.people[n].T;                // LM detectable infections
        POP.prev_all[3] = POP.prev_all[3] + POP.people[n].I_D + POP.people[n].T;                                      // Clinical episodes

        if (POP.people[n].Hyp > 0)
        {
            POP.prev_all[4] = POP.prev_all[4] + 1;                     // Hypnozoite positive

            POP.prev_all[5] = POP.prev_all[5] + POP.people[n].Hyp;    // Number of batches of hypnozoites
        }


        ////////////////////////////////////////
        // Incidence

        POP.prev_all[6]  = POP.prev_all[6]  + POP.people[n].I_PCR_new;
        POP.prev_all[7]  = POP.prev_all[7]  + POP.people[n].I_LM_new;
        POP.prev_all[8]  = POP.prev_all[8]  + POP.people[n].I_D_new;
        POP.prev_all[9]  = POP.prev_all[9]  + POP.people[n].ACT_treat;
        POP.prev_all[10] = POP.prev_all[10] + POP.people[n].PQ_treat;


        //////////////////////////////////////////////
        //////////////////////////////////////////////
        // Summary - under 5's

        if (POP.people[n].age < 1825.0)
        {
            ////////////////////////////////////////
            // Prevalence

            POP.prev_U5[0] = POP.prev_U5[0] + 1;                                                                // Numbers - denominator
            POP.prev_U5[1] = POP.prev_U5[1] + POP.people[n].I_PCR + POP.people[n].I_LM
                                              + POP.people[n].I_D + POP.people[n].T;                              // PCR detectable infections
            POP.prev_U5[2] = POP.prev_U5[2] + POP.people[n].I_LM + POP.people[n].I_D + POP.people[n].T;        // LM detectable infections
            POP.prev_U5[3] = POP.prev_U5[3] + POP.people[n].I_D + POP.people[n].T;                              // Clinical episodes

            if (POP.people[n].Hyp > 0)
            {
                POP.prev_U5[4] = POP.prev_U5[4] + 1;                     // Hypnozoite positive

                POP.prev_U5[5] = POP.prev_U5[5] + POP.people[n].Hyp;    // Number of batches of hypnozoites
            }


            ////////////////////////////////////////
            // Incidence

            POP.prev_U5[6] = POP.prev_U5[6] + POP.people[n].I_PCR_new;
            POP.prev_U5[7] = POP.prev_U5[7] + POP.people[n].I_LM_new;
            POP.prev_U5[8] = POP.prev_U5[8] + POP.people[n].I_D_new;
            POP.prev_U5[9] = POP.prev_U5[9] + POP.people[n].ACT_treat;
            POP.prev_U5[10] = POP.prev_U5[10] + POP.people[n].PQ_treat;
        }

        //////////////////////////////////////////////
        //////////////////////////////////////////////
        // Summary - under 10's

        if (POP.people[n].age < 3650.0)
        {
            ////////////////////////////////////////
            // Prevalence

            POP.prev_U10[0] = POP.prev_U10[0] + 1;                                                            // Numbers - denominator
            POP.prev_U10[1] = POP.prev_U10[1] + POP.people[n].I_PCR + POP.people[n].I_LM
                                                + POP.people[n].I_D + POP.people[n].T;                          // PCR detectable infections
            POP.prev_U10[2] = POP.prev_U10[2] + POP.people[n].I_LM + POP.people[n].I_D + POP.people[n].T;    // LM detectable infections
            POP.prev_U10[3] = POP.prev_U10[3] + POP.people[n].I_D + POP.people[n].T;                          // Clinical episodes

            if (POP.people[n].Hyp > 0)
            {
                POP.prev_U10[4] = POP.prev_U10[4] + 1;                     // Hypnozoite positive

                POP.prev_U10[5] = POP.prev_U10[5] + POP.people[n].Hyp;    // Number of batches of hypnozoites
            }


            ////////////////////////////////////////
            // Incidence

            POP.prev_U10[6]  = POP.prev_U10[6]  + POP.people[n].I_PCR_new;
            POP.prev_U10[7]  = POP.prev_U10[7]  + POP.people[n].I_LM_new;
            POP.prev_U10[8]  = POP.prev_U10[8]  + POP.people[n].I_D_new;
            POP.prev_U10[9]  = POP.prev_U10[9]  + POP.people[n].ACT_treat;
            POP.prev_U10[10] = POP.prev_U10[10] + POP.people[n].PQ_treat;
        }
    }


    //////////////////////////////
    // Intervention coverage

    POP.LLIN_cov_t = 0;
    POP.IRS_cov_t = 0;
    POP.ACT_treat_t = 0;
    POP.PQ_treat_t = 0;
    POP.pregnant_t = 0;

    POP.PQ_overtreat_t = 0;
    POP.PQ_overtreat_9m_t = 0;


    for (int n = 0; n<POP.N_pop; n++)
    {
        POP.LLIN_cov_t  = POP.LLIN_cov_t  + POP.people[n].LLIN;
        POP.IRS_cov_t   = POP.IRS_cov_t   + POP.people[n].IRS;
        POP.ACT_treat_t = POP.ACT_treat_t + POP.people[n].ACT_treat;
        POP.PQ_treat_t  = POP.PQ_treat_t  + POP.people[n].PQ_treat;
        POP.pregnant_t  = POP.pregnant_t  + POP.people[n].pregnant;

        POP.PQ_overtreat_t    = POP.PQ_overtreat_t    + POP.people[n].PQ_overtreat;
        POP.PQ_overtreat_9m_t = POP.PQ_overtreat_9m_t + POP.people[n].PQ_overtreat_9m;
    }


    //////////////////////////////
    // Immunity

    double A_par_mean = 0.0, A_clin_mean = 0.0;

    for (int n = 0; n<POP.N_pop; n++)
    {
        A_par_mean = A_par_mean + POP.people[n].A_par;
        A_clin_mean = A_clin_mean + POP.people[n].A_clin;
    }

    POP.A_par_mean_t = A_par_mean / ((double)POP.N_pop);
    POP.A_clin_mean_t = A_clin_mean / ((double)POP.N_pop);
}



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.6. Simulate the model and store the output in SIM                     //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void model_simulator(Params& theta, Population& POP, Intervention& INTVEN, Simulation& SIM)
{

    for (int i = 0; i<SIM.N_time; i++)
    {
        if (SIM.t_vec[i] / 365.0 - floor(SIM.t_vec[i] / 365.0) < 0.5*t_step / 365.0)
        {
            cout << "time = " << SIM.t_vec[i] / 365.0 << "\t" << 100.0*(SIM.t_vec[i] - SIM.t_vec[0]) / (double(t_step*SIM.N_time)) << "% complete" << endl;
        }

        human_step(theta, POP);

        mosquito_step(SIM.t_vec[i], theta, POP);

        intervention_dist(SIM.t_vec[i], theta, POP, INTVEN);

        POP_summary(POP);

        //////////////////////////////////////
        // Fill out Simulation object

        for (int k = 0; k<N_H_comp; k++)
        {
            SIM.yH_t[i][k] = POP.yH[k];
        }

        for (int k = 0; k<N_M_comp; k++)
        {
            for (int g = 0; g < N_spec; g++)
            {
                SIM.yM_t[i][g][k] = POP.yM[g][k];
            }
        }

        for (int k = 0; k<11; k++)
        {
            SIM.prev_all[i][k] = POP.prev_all[k];
            SIM.prev_U5[i][k] = POP.prev_U5[k];
            SIM.prev_U10[i][k] = POP.prev_U10[k];
        }


        SIM.LLIN_cov_t[i] = POP.LLIN_cov_t;
        SIM.IRS_cov_t[i] = POP.IRS_cov_t;
        SIM.ACT_treat_t[i] = POP.ACT_treat_t;
        SIM.PQ_treat_t[i] = POP.PQ_treat_t;
        SIM.pregnant_t[i] = POP.pregnant_t;

        SIM.PQ_overtreat_t[i] = POP.PQ_overtreat_t;
        SIM.PQ_overtreat_9m_t[i] = POP.PQ_overtreat_9m_t;


        SIM.EIR_t[i] = 0.0;
        for (int g = 0; g < N_spec; g++)
        {
            SIM.EIR_t[i] = SIM.EIR_t[i] + POP.aa_VC[g] * POP.yM[g][5];
        }

        SIM.A_par_mean_t[i] = POP.A_par_mean_t;
        SIM.A_clin_mean_t[i] = POP.A_clin_mean_t;

    }

}



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2.7. Vector control distribution                                        //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////


void intervention_dist(double t, Params& theta, Population& POP, Intervention& INTVEN)
{
    double QQ;

    bool BS_effective;
    bool PQ_treat;
    bool PQ_effective;

    bool MSAT_pos;
    bool SSAT_pos;

    //////////////////////////////////////////////////////////
    // Intervention 1: LLINS

    for (size_t m = 0; m<INTVEN.LLIN_year.size(); m++)
    {
        if ((t > INTVEN.LLIN_year[m] - 0.5*t_step) &&
            (t < INTVEN.LLIN_year[m] + 0.51*t_step))
        {
            cout << "LLIN distribution" << endl;

            try 
            {
                QQ = phi_inv(INTVEN.LLIN_cover[m], 0.0, sqrt(1.0 + theta.sig_round_LLIN*theta.sig_round_LLIN));
            }
            catch (const char* e)
            {
                std::cerr << e << std::endl;
                exit (1);
            }

            for (int n = 0; n<POP.N_pop; n++)
            {
                if (gennor(POP.people[n].zz_int[0], theta.sig_round_LLIN) < QQ)
                {
                    POP.people[n].LLIN = 1;
                    POP.people[n].LLIN_age = 0.0;

                    for (int g = 0; g < N_spec; g++)
                    {
                        POP.people[n].d_LLIN[g] = theta.d_LLIN_0[g];
                        POP.people[n].r_LLIN[g] = theta.r_LLIN_0[g];
                        POP.people[n].s_LLIN[g] = 1.0 - POP.people[n].d_LLIN[g] - POP.people[n].r_LLIN[g];
                    }
                }
            }
        }
    }


    //////////////////////////////////////////////////////////
    // Intervention 2: IRS

    for (size_t m = 0; m<INTVEN.IRS_year.size(); m++)
    {
        if ((t > INTVEN.IRS_year[m] - 0.5*t_step) &&
            (t < INTVEN.IRS_year[m] + 0.51*t_step))
        {
            cout << "IRS distribution" << endl;

            try 
            {
                QQ = phi_inv(INTVEN.IRS_cover[m], 0.0, sqrt(1.0 + theta.sig_round_IRS*theta.sig_round_IRS));
            }
            catch (const char* e)
            {
                std::cerr << e << std::endl;
                exit (1);
            }

            for (int n = 0; n<POP.N_pop; n++)
            {
                if (gennor(POP.people[n].zz_int[1], theta.sig_round_IRS) < QQ)
                {
                    POP.people[n].IRS = 1;
                    POP.people[n].IRS_age = 0.0;

                    for (int g = 0; g < N_spec; g++)
                    {
                        POP.people[n].d_IRS[g] = theta.d_IRS_0[g];
                        POP.people[n].r_IRS[g] = theta.r_IRS_0[g];
                        POP.people[n].s_IRS[g] = 1.0 - POP.people[n].d_IRS[g] - POP.people[n].r_IRS[g];
                    }
                }
            }
        }
    }


    //////////////////////////////////////////////////////////
    // Intervention 3: first-line treatment (blood-stage)

    for (size_t m = 0; m<INTVEN.BS_treat_year_on.size(); m++)
    {
        if ((t > INTVEN.BS_treat_year_on[m] - 0.5*t_step) &&
            (t < INTVEN.BS_treat_year_on[m] + 0.51*t_step))
        {
            cout << "New front-line BS treatment" << endl;

            theta.BS_treat_BScover = INTVEN.BS_treat_BScover[m];
            theta.BS_treat_BSeff   = INTVEN.BS_treat_BSeff[m];
            theta.BS_treat_BSproph = INTVEN.BS_treat_BSproph[m];

            theta.treat_BScover = theta.BS_treat_BScover;
            theta.treat_BSeff   = theta.BS_treat_BSeff;
            theta.treat_PQavail = 0.0;
            theta.r_P           = 1.0 / theta.BS_treat_BSproph;
        }
    }


    //////////////////////////////////////////////////////////
    // Switching back to baseline.

    for (size_t m = 0; m<INTVEN.BS_treat_year_on.size(); m++)
    {
        if ((t > INTVEN.BS_treat_year_off[m] - 0.5*t_step) &&
            (t < INTVEN.BS_treat_year_off[m] + 0.51*t_step))
        {
            cout << "End of changing front-line BS treatment" << endl;

            theta.treat_BScover = theta.BS_treat_BScover_base;
            theta.treat_BSeff   = theta.BS_treat_BSeff_base;
            theta.treat_PQavail = 0.0;
            theta.r_P           = 1.0 / theta.BS_treat_BSproph_base;
        }
    }


    //////////////////////////////////////////////////////////
    // Intervention 4: first-line treatment (primaquine)

    for (size_t m = 0; m<INTVEN.PQ_treat_year_on.size(); m++)
    {
        if ((t > INTVEN.PQ_treat_year_on[m] - 0.5*t_step) &&
            (t < INTVEN.PQ_treat_year_on[m] + 0.51*t_step))
        {
            cout << "New front-line PQ treatment" << endl;

            theta.PQ_treat_BScover     = INTVEN.PQ_treat_BScover[m];
            theta.PQ_treat_BSeff       = INTVEN.PQ_treat_BSeff[m];
            theta.PQ_treat_BSproph     = INTVEN.PQ_treat_BSproph[m];
            theta.PQ_treat_PQavail     = INTVEN.PQ_treat_PQavail[m];
            theta.PQ_treat_PQeff       = INTVEN.PQ_treat_PQeff[m];
            theta.PQ_treat_PQproph     = INTVEN.PQ_treat_PQproph[m];
            theta.PQ_treat_G6PD_risk   = INTVEN.PQ_treat_G6PD_risk[m];
            theta.PQ_treat_CYP2D6_risk = INTVEN.PQ_treat_CYP2D6_risk[m];
            theta.PQ_treat_preg_risk   = INTVEN.PQ_treat_preg_risk[m];
            theta.PQ_treat_low_age     = INTVEN.PQ_treat_low_age[m];

            theta.treat_BScover = theta.PQ_treat_BScover;
            theta.treat_BSeff   = theta.PQ_treat_BSeff;
            theta.treat_PQavail = theta.PQ_treat_PQavail;
            theta.r_P           = 1.0 / theta.PQ_treat_BSproph;
        }
    }


    //////////////////////////////////////////////////////////
    // Switching back to baseline.

    for (size_t m = 0; m<INTVEN.PQ_treat_year_on.size(); m++)
    {
        if ((t > INTVEN.PQ_treat_year_off[m] - 0.5*t_step) &&
            (t < INTVEN.PQ_treat_year_off[m] + 0.51*t_step))
        {
            cout << "End of changing front-line PQ treatment" << endl;

            theta.PQ_treat_BScover     = 0.0;
            theta.PQ_treat_BSeff       = 0.0;
            theta.PQ_treat_BSproph     = 10.0;
            theta.PQ_treat_PQavail     = 0.0;
            theta.PQ_treat_PQeff       = 0.0;
            theta.PQ_treat_PQproph     = 10.0;
            theta.PQ_treat_G6PD_risk   = 1;
            theta.PQ_treat_CYP2D6_risk = 1;
            theta.PQ_treat_preg_risk   = 1;
            theta.PQ_treat_low_age     = 180.0;

            theta.treat_BScover = theta.BS_treat_BScover_base;
            theta.treat_BSeff   = theta.BS_treat_BSeff_base;
            theta.treat_PQavail = 0.0;
            theta.r_P           = 1.0 / theta.BS_treat_BSproph_base;
        }
    }


    //////////////////////////////////////////////////////////
    // Intervention 5: MDA (blood-stage)

    for (size_t m = 0; m<INTVEN.MDA_BS_year.size(); m++)
    {
        if ((t > INTVEN.MDA_BS_year[m] - 0.5*t_step) &&
            (t < INTVEN.MDA_BS_year[m] + 0.51*t_step))
        {
            cout << "MDA (BS) distribution" << endl;

            theta.MDA_BS_BScover = INTVEN.MDA_BS_BScover[m];
            theta.MDA_BS_BSeff   = INTVEN.MDA_BS_BSeff[m];
            theta.MDA_BS_BSproph = INTVEN.MDA_BS_BSproph[m];

            try 
            {
                QQ = phi_inv(theta.MDA_BS_BScover, 0.0, sqrt(1.0 + theta.sig_round_MDA*theta.sig_round_MDA));
            }
            catch (const char* e)
            {
                std::cerr << e << std::endl;
                exit (1);
            }

            for (int n = 0; n<POP.N_pop; n++)
            {
                if (gennor(POP.people[n].zz_int[2], theta.sig_round_MDA) < QQ)
                {
                    POP.people[n].ACT_treat = 1;

                    if (genunf(0.0, 1.0) < theta.MDA_BS_BSeff)
                    {
                        if (POP.people[n].S == 1    ) { POP.people[n].S = 0;     POP.people[n].P = 1; }
                        if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P = 1; }
                        if (POP.people[n].I_LM == 1 ) { POP.people[n].I_LM = 0;  POP.people[n].P = 1; }
                        if (POP.people[n].I_D == 1  ) { POP.people[n].I_D = 0;   POP.people[n].T = 1; }
                    }
                }
            }
        }
    }


    //////////////////////////////////////////////////////////
    // Intervention 6: MDA (blood-stage and liver-stage)

    for (size_t m = 0; m<INTVEN.MDA_PQ_year.size(); m++)
    {
        if ((t > INTVEN.MDA_PQ_year[m] - 0.5*t_step) &&
            (t < INTVEN.MDA_PQ_year[m] + 0.51*t_step))
        {
            cout << "MDA (BS+PQ) distribution" << endl;

            theta.MDA_PQ_BScover     = INTVEN.MDA_PQ_BScover[m];
            theta.MDA_PQ_BSeff       = INTVEN.MDA_PQ_BSeff[m];
            theta.MDA_PQ_BSproph     = INTVEN.MDA_PQ_BSproph[m];
            theta.MDA_PQ_PQavail     = INTVEN.MDA_PQ_PQavail[m];
            theta.MDA_PQ_PQeff       = INTVEN.MDA_PQ_PQeff[m];
            theta.MDA_PQ_PQproph     = INTVEN.MDA_PQ_PQproph[m];
            theta.MDA_PQ_G6PD_risk   = INTVEN.MDA_PQ_G6PD_risk[m];
            theta.MDA_PQ_CYP2D6_risk = INTVEN.MDA_PQ_CYP2D6_risk[m];
            theta.MDA_PQ_preg_risk   = INTVEN.MDA_PQ_preg_risk[m];
            theta.MDA_PQ_low_age     = INTVEN.MDA_PQ_low_age[m];

            try 
            {
                QQ = phi_inv(theta.MDA_PQ_BScover, 0.0, sqrt(1.0 + theta.sig_round_MDA*theta.sig_round_MDA));
            }
            catch (const char* e)
            {
                std::cerr << e << std::endl;
                exit (1);
            }

            for (int n = 0; n<POP.N_pop; n++)
            {
                if (gennor(POP.people[n].zz_int[3], theta.sig_round_MDA) < QQ)
                {
                    /////////////////////////////////////////////////////
                    // Blood-stage treatment is always administered
                    // Is blood-stage treatment effective

                    BS_effective = 0;

                    if (genunf(0.0, 1.0) < theta.MDA_PQ_BSeff)
                    {
                        BS_effective = 1;
                    }


                    /////////////////////////////////////////////////////////////////////
                    // Is PQ administered?

                    PQ_treat = 0;

                    if( genunf(0.0, 1.0) < theta.MDA_PQ_PQavail )
                    {
                        PQ_treat = 1;
                    }
                    

                    /////////////////////////////////////////////////////////////////////
                    // Exclude PQ because of G6PD deficiency

                    if( (theta.MDA_PQ_G6PD_risk == 1) && (POP.people[n].G6PD_def == 1) )
                    {
                        PQ_treat = 0;
                    }


                    /////////////////////////////////////////////////////////////////////
                    // Exclude PQ because of pregancy

                    if( (theta.MDA_PQ_preg_risk == 1) && (POP.people[n].pregnant == 1) )
                    {
                        PQ_treat = 0;
                    }


                    /////////////////////////////////////////////////////////////////////
                    // Exclude PQ because of young age

                    if (POP.people[n].age < theta.MDA_PQ_low_age)
                    {
                        PQ_treat = 0;
                    }

                    if (PQ_treat == 1)
                    {
                        POP.people[n].PQ_treat = 1;
                    }

                    /////////////////////////////////////////////////////////////////////
                    // Is PQ effective?

                    PQ_effective = 1;

                    if( genunf(0.0, 1.0) > theta.MDA_PQ_PQeff )
                    {
                        PQ_effective = 0;
                    }

                    if( (theta.MDA_PQ_CYP2D6_risk == 1) && (POP.people[n].CYP2D6 == 1) )
                    {
                        PQ_effective = 0;
                    }


                    /////////////////////////////////////////////////////////////////////
                    // Was there PQ overtreatment?

                    if( (PQ_treat == 1) && (POP.people[n].Hyp == 0) )
                    {
                        POP.people[n].PQ_overtreat = 1;
                    }

                    if( (PQ_treat == 1) && (POP.people[n].T_last_BS > 270.0) )
                    {
                        POP.people[n].PQ_overtreat_9m = 1;
                    }


                    /////////////////////////////////////////////////////////////////////
                    // ACTION: administer blood-stage drug

                    POP.people[n].ACT_treat = 1;

                    if (BS_effective == 1)
                    {
                        if (POP.people[n].S == 1) {     POP.people[n].S = 0;     POP.people[n].P = 1; }
                        if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P = 1; }
                        if (POP.people[n].I_LM == 1) {  POP.people[n].I_LM = 0;  POP.people[n].P = 1; }
                        if (POP.people[n].I_D == 1) {   POP.people[n].I_D = 0;   POP.people[n].T = 1; }
                    }


                    /////////////////////////////////////////////////////////////////////
                    // ACTION: administer primaquine

                    if ((PQ_treat == 1) && (PQ_effective == 1))
                    {
                        POP.people[n].Hyp = 0;

                        POP.people[n].PQ_proph = 1;
                        POP.people[n].PQ_proph_timer = theta.MDA_PQ_PQproph;
                    }
                }
            }
        }
    }


    //////////////////////////////////////////////////////////
    // Intervention 7: MSAT (blood-stage and liver-stage)

    for (size_t m = 0; m < INTVEN.MSAT_PQ_year.size(); m++)
    {
        if( (t > INTVEN.MSAT_PQ_year[m] - 0.5*t_step) &&
            (t < INTVEN.MSAT_PQ_year[m] + 0.51*t_step) )
        {
            cout << "MSAT (BS+PQ) distribution" << endl;

            theta.MSAT_PQ_BScover     = INTVEN.MSAT_PQ_BScover[m];
            theta.MSAT_PQ_RDT_PCR     = INTVEN.MSAT_PQ_RDT_PCR[m];
            theta.MSAT_PQ_sens        = INTVEN.MSAT_PQ_sens[m];
            theta.MSAT_PQ_BSeff       = INTVEN.MSAT_PQ_BSeff[m];
            theta.MSAT_PQ_BSproph     = INTVEN.MSAT_PQ_BSproph[m];
            theta.MSAT_PQ_PQavail     = INTVEN.MSAT_PQ_PQavail[m];
            theta.MSAT_PQ_PQeff       = INTVEN.MSAT_PQ_PQeff[m];
            theta.MSAT_PQ_PQproph     = INTVEN.MSAT_PQ_PQproph[m];
            theta.MSAT_PQ_G6PD_risk   = INTVEN.MSAT_PQ_G6PD_risk[m];
            theta.MSAT_PQ_CYP2D6_risk = INTVEN.MSAT_PQ_CYP2D6_risk[m];
            theta.MSAT_PQ_preg_risk   = INTVEN.MSAT_PQ_preg_risk[m];
            theta.MSAT_PQ_low_age     = INTVEN.MSAT_PQ_low_age[m];

            try 
            {
                QQ = phi_inv(theta.MSAT_PQ_BScover, 0.0, sqrt(1.0 + theta.sig_round_MDA*theta.sig_round_MDA));
            }
            catch (const char* e)
            {
                std::cerr << e << std::endl;
                exit (1);
            }

            for (int n = 0; n < POP.N_pop; n++)
            {
                if (gennor(POP.people[n].zz_int[4], theta.sig_round_MDA) < QQ)
                {
                    /////////////////////////////////////////////////////
                    // Blood-stage treatment is always administered

                    MSAT_pos = 0;

                    ////////////////////////////////////////////////
                    // Diagnosis by RDT, assumed the same as LM 

                    if (theta.MSAT_PQ_RDT_PCR == 1)
                    {
                        if ((POP.people[n].I_LM == 1) || (POP.people[n].I_D == 1) || (POP.people[n].T == 1))
                        {
                            if (genunf(0.0, 1.0) < theta.MSAT_PQ_sens)
                            {
                                MSAT_pos = 1;
                            }
                        }
                    }


                    ////////////////////////////////////////////////
                    // Diagnosis by PCR 

                    if (theta.MSAT_PQ_RDT_PCR == 2)
                    {
                        if ((POP.people[n].I_PCR == 1) || (POP.people[n].I_LM == 1) || (POP.people[n].I_D == 1) || (POP.people[n].T == 1))
                        {
                            if (genunf(0.0, 1.0) < theta.MSAT_PQ_sens)
                            {
                                MSAT_pos = 1;
                            }
                        }
                    }


                    /////////////////////////////////////////////////////
                    // Is blood-stage treatment effective

                    BS_effective = 0;

                    if (genunf(0.0, 1.0) < theta.MSAT_PQ_BSeff)
                    {
                        BS_effective = 1;
                    }


                    /////////////////////////////////////////////////////////////////////
                    // Is PQ administered?

                    PQ_treat = 0;

                    if (MSAT_pos == 1)
                    {
                        if (genunf(0.0, 1.0) < theta.MSAT_PQ_PQavail)
                        {
                            PQ_treat = 1;
                        }
                    }


                    /////////////////////////////////////////////////////////////////////
                    // Exclude PQ because of G6PD deficiency

                    if ( (theta.MSAT_PQ_G6PD_risk == 1) && (POP.people[n].G6PD_def == 1) )
                    {
                        PQ_treat = 0;
                    }

                    /////////////////////////////////////////////////////////////////////
                    // Exclude PQ because of pregancy

                    if ( (theta.MSAT_PQ_preg_risk == 1) && (POP.people[n].pregnant == 1) )
                    {
                        PQ_treat = 0;
                    }

                    /////////////////////////////////////////////////////////////////////
                    // Exclude PQ because of young age

                    if (POP.people[n].age < theta.MSAT_PQ_low_age)
                    {
                        PQ_treat = 0;
                    }

                    if (PQ_treat == 1)
                    {
                        POP.people[n].PQ_treat = 1;
                    }

                    /////////////////////////////////////////////////////////////////////
                    // Is PQ effective?

                    PQ_effective = 1;

                    if (genunf(0.0, 1.0) > theta.MSAT_PQ_PQeff)
                    {
                        PQ_effective = 0;
                    }

                    if ((theta.MSAT_PQ_CYP2D6_risk == 1) && (POP.people[n].CYP2D6 == 1))
                    {
                        PQ_effective = 0;
                    }


                    /////////////////////////////////////////////////////////////////////
                    // Was there PQ overtreatment?

                    if ((PQ_treat == 1) && (POP.people[n].Hyp == 0))
                    {
                        POP.people[n].PQ_overtreat = 1;
                    }

                    if ((PQ_treat == 1) && (POP.people[n].T_last_BS > 270.0))
                    {
                        POP.people[n].PQ_overtreat_9m = 1;
                    }


                    /////////////////////////////////////////////////////////////////////
                    // ACTION: administer blood-stage drug

                    if (MSAT_pos == 1)
                    {
                        POP.people[n].ACT_treat = 1;

                        if (BS_effective == 1)
                        {
                            if (POP.people[n].S == 1) {     POP.people[n].S = 0;     POP.people[n].P = 1; }
                            if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P = 1; }
                            if (POP.people[n].I_LM == 1) {  POP.people[n].I_LM = 0;  POP.people[n].P = 1; }
                            if (POP.people[n].I_D == 1) {   POP.people[n].I_D = 0;   POP.people[n].T = 1; }
                        }
                    }

                    
                    /////////////////////////////////////////////////////////////////////
                    // ACTION: administer primaquine

                    if ((PQ_treat == 1) && (PQ_effective == 1))
                    {
                        POP.people[n].Hyp = 0;

                        POP.people[n].PQ_proph = 1;
                        POP.people[n].PQ_proph_timer = theta.MSAT_PQ_PQproph;
                    }
                }
            }
        }
    }


    //////////////////////////////////////////////////////////
    // Intervention 8: SSAT (blood-stage and liver-stage)

    for (size_t m = 0; m<INTVEN.SSAT_PQ_year.size(); m++)
    {
        if ((t > INTVEN.SSAT_PQ_year[m] - 0.5*t_step) &&
            (t < INTVEN.SSAT_PQ_year[m] + 0.51*t_step))
        {
            cout << "SSAT (BS+PQ) distribution" << endl;

            theta.SSAT_PQ_BScover     = INTVEN.SSAT_PQ_BScover[m];
            theta.SSAT_PQ_sens        = INTVEN.SSAT_PQ_sens[m];
            theta.SSAT_PQ_spec        = INTVEN.SSAT_PQ_spec[m];
            theta.SSAT_PQ_BSeff       = INTVEN.SSAT_PQ_BSeff[m];
            theta.SSAT_PQ_BSproph     = INTVEN.SSAT_PQ_BSproph[m];
            theta.SSAT_PQ_PQavail     = INTVEN.SSAT_PQ_PQavail[m];
            theta.SSAT_PQ_PQeff       = INTVEN.SSAT_PQ_PQeff[m];
            theta.SSAT_PQ_PQproph     = INTVEN.SSAT_PQ_PQproph[m];
            theta.SSAT_PQ_G6PD_risk   = INTVEN.SSAT_PQ_G6PD_risk[m];
            theta.SSAT_PQ_CYP2D6_risk = INTVEN.SSAT_PQ_CYP2D6_risk[m];
            theta.SSAT_PQ_preg_risk   = INTVEN.SSAT_PQ_preg_risk[m];
            theta.SSAT_PQ_low_age     = INTVEN.SSAT_PQ_low_age[m];

            try 
            {
                QQ = phi_inv(theta.SSAT_PQ_BScover, 0.0, sqrt(1.0 + theta.sig_round_MDA*theta.sig_round_MDA));
            }
            catch (const char* e)
            {
                std::cerr << e << std::endl;
                exit (1);
            }

            for (int n = 0; n < POP.N_pop; n++)
            {
                if (gennor(POP.people[n].zz_int[5], theta.sig_round_MDA) < QQ)
                {
                    /////////////////////////////////////////////////////
                    // Blood-stage treatment is always administered

                    POP.people[n].ACT_treat = 1;


                    /////////////////////////////////////////////////////
                    // Is blood-stage treatment effective

                    BS_effective = 0;

                    if( genunf(0.0, 1.0) < theta.SSAT_PQ_BSeff )
                    {
                        BS_effective = 1;
                    }

                    /////////////////////////////////////////////////////////////////////
                    // SSAT screening for blood-stage infection in the last 9 months
                    //
                    // There are two options here. Option 1 define over-treatment on the 
                    // basis of blood-stage infection with the last 9 month. Option 2
                    // defines over-treatment on the basis of presence of hypnozoites.

                    SSAT_pos = 0;

                    // OPTION 1
/*
                    if( (POP.people[n].T_last_BS <= 270.0) && (genunf(0.0, 1.0) < theta.SSAT_PQ_sens) )
                    {
                        SSAT_pos = 1;
                    }

                    if( (POP.people[n].T_last_BS > 270.0) && (genunf(0.0, 1.0) > theta.SSAT_PQ_spec) )
                    {
                        SSAT_pos = 1;
                    }
*/

                    // OPTION 2

                    if ((POP.people[n].Hyp > 0) && (genunf(0.0, 1.0) < theta.SSAT_PQ_sens))
                    {
                        SSAT_pos = 1;
                    }

                    if ((POP.people[n].Hyp == 0) && (genunf(0.0, 1.0) > theta.SSAT_PQ_spec))
                    {
                        SSAT_pos = 1;
                    }


                    /////////////////////////////////////////////////////////////////////
                    // Is PQ administered?

                    PQ_treat = 0;

                    if( SSAT_pos == 1 )
                    {
                        if( genunf(0.0, 1.0) < theta.SSAT_PQ_PQavail )
                        {
                            PQ_treat = 1;
                        }
                    }


                    /////////////////////////////////////////////////////////////////////
                    // Exclude PQ because of G6PD deficiency

                    if( (theta.SSAT_PQ_G6PD_risk == 1) && (POP.people[n].G6PD_def == 1) )
                    {
                        PQ_treat = 0;
                    }

                    /////////////////////////////////////////////////////////////////////
                    // Exclude PQ because of pregancy

                    if( (theta.SSAT_PQ_preg_risk == 1) && (POP.people[n].pregnant == 1) )
                    {
                        PQ_treat = 0;
                    }

                    /////////////////////////////////////////////////////////////////////
                    // Exclude PQ because of young age

                    if( POP.people[n].age < theta.SSAT_PQ_low_age )
                    {
                        PQ_treat = 0;
                    }

                    if( PQ_treat == 1 )
                    {
                        POP.people[n].PQ_treat = 1;
                    }

                    /////////////////////////////////////////////////////////////////////
                    // Is PQ effective?

                    PQ_effective = 0;

                    if( genunf(0.0, 1.0) < theta.SSAT_PQ_PQeff )
                    {
                        PQ_effective = 1;
                    }

                    if( (theta.SSAT_PQ_CYP2D6_risk == 1) && (POP.people[n].CYP2D6 == 1) )
                    {
                        PQ_effective = 0;
                    }


                    /////////////////////////////////////////////////////////////////////
                    // Was there PQ overtreatment?

                    if ((PQ_treat == 1) && (POP.people[n].Hyp == 0))
                    {
                        POP.people[n].PQ_overtreat = 1;
                    }

                    if ((PQ_treat == 1) && (POP.people[n].T_last_BS > 270.0))
                    {
                        POP.people[n].PQ_overtreat_9m = 1;
                    }


                    /////////////////////////////////////////////////////////////////////
                    // ACTION: administer blood-stage drug

                    if( BS_effective == 1 )
                    {
                        if (POP.people[n].S     == 1) { POP.people[n].S = 0;     POP.people[n].P = 1; }
                        if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P = 1; }
                        if (POP.people[n].I_LM  == 1) { POP.people[n].I_LM = 0;  POP.people[n].P = 1; }
                        if (POP.people[n].I_D   == 1) { POP.people[n].I_D = 0;   POP.people[n].T = 1; }
                    }

                    /////////////////////////////////////////////////////////////////////
                    // ACTION: administer primaquine

                    if( (PQ_treat == 1) && (PQ_effective == 1) )
                    {
                        POP.people[n].Hyp = 0;

                        POP.people[n].PQ_proph = 1;
                        POP.people[n].PQ_proph_timer = theta.SSAT_PQ_PQproph;
                    }
                }
            }


        }
    }

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
    if (pp < 0.0 || pp > 1.0)
    {
        throw("bad vlaue of pp (coverage) in phi_inv");
    }
    else if (pp == 0.0) 
    {
        return -std::numeric_limits<double>::infinity();
    }
    else if (pp == 1.0)
    {
        return std::numeric_limits<double>::infinity();
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
