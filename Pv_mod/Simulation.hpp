/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///                                                                       ///
///  Individual-based Plasmodium vivax transmission model.                ///
///                                                                       ///
///  Dr Michael White                                                     ///
///  Institut Pasteur                                                     ///
///  michael.white@pasteur.fr                                             ///
///                                                                       ///
///  This file and accompanying .cpp contain:                             ///
///                                                                       ///
///  1. MAIN - SIMULATION                                                 ///
///     Here we read in parameters from files (model parameters,          ///
///     mosquito parameters, intervention parameters).                    ///
///     A population of individuals is created at equilibrium and         ///
///     then simulated.                                                   ///
///                                                                       ///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Include guard
#ifndef PVIVAX_MODEL_SIMULATION
#define PVIVAX_MODEL_SIMULATION

#include "Intervention.hpp"


/////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
// 0.5. Define a structure for storing the output of a simulation                      //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

class Simulation
{
public:
    //////////////////////////////////////////////////////////////////////////
    //  Functions
    //////////////////////////////////////////////////////////////////////////
    
    
    /////////////////////////////////////
    // Set up simulation times and storage for outputs
    Simulation(SimTimes times);
    
    /////////////////////////////////////
    //  2.6. Simulate the model and store the output
    void run(Params& theta, Population& POP, Intervention& INTVEN);

    /////////////////////////////////////
    // Write output file
    void write_output(const char *output_File);
    

    //////////////////////////////////////////////////////////////////////////
    //  Data
    //////////////////////////////////////////////////////////////////////////

    
    //////////////////////////////////////////
    // 0.5.1. Vector of simulation times

    int N_time;
    SimTimes times;

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

#endif
