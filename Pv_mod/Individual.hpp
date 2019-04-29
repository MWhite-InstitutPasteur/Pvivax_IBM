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
///  4. INDIVIDUAL-BASED MODEL                                            ///
///     Details of the stochastic individual-based model for each         ///
///     person. Transitions occur with a fixed time step according to     ///
///     compting hazards                                                  ///
///                                                                       ///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Include guard
#ifndef PVIVAX_MODEL_INDIVIDUAL
#define PVIVAX_MODEL_INDIVIDUAL

#include "Params.hpp"


/////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
// 0.2. Define a  for humans                                                      //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

class Individual
{
public:
    //////////////////////////////////////////////////////////////////////////
    //  Class constructors and destructors
    //////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////
    // 0.2.1. Class constructor
    Individual(double a, double zeta)
    {
        age = a;
        zeta_het = zeta;
    }
    
    
    ////////////////////////////////////////////////////
    // Copy and move constructors

    // Delete unwanted copy constructors
    Individual(Individual&) =delete;
    void operator= (Individual&) =delete;
    // Allow default move constructors
    Individual(Individual&&) = default;
    Individual& operator= (Individual&&) = default;


    //////////////////////////////////////////////////////////////////////////
    //  Class member functions
    //////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////
    // 0.2.2. Function declarations within the human class

    void state_mover(Params& theta, double lam_bite);
    void ager(Params& theta);
    void intervention_updater(Params& theta);


    //////////////////////////////////////////////////////////////////////////
    //  Data
    //////////////////////////////////////////////////////////////////////////

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

    double T_last_BS;         // Tracking of time since last PCR-detectable blood-stage infection

    ///////////////////////////////////////////////////
    // Individual-level effect of vector control

    double z_VC[N_spec];      // probability of mosquito being repelled from this individual during a single feeding attempt
    double y_VC[N_spec];      // probability of mosquito feeding on this individual during a single attempt
    double w_VC[N_spec];      // probability of mosquito feeding and surviving on this individual during a single feeding attempt
};

#endif
