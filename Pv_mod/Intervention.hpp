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
///     A structure called Population stores all individuals.             ///
///                                                                       ///
///  3. EQUILIBRIUM SETUP                                                 ///
///     This set of functions calculates the equilibrium set up of the    ///
///     population. It is only called once while the population is        ///
///     being initialised.                                                ///
///                                                                       ///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Include guard
#ifndef PVIVAX_MODEL_INTERVENTION
#define PVIVAX_MODEL_INTERVENTION

#include "Population.hpp"


/////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
// 0.4. Define a structure for details of interventions                                //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

class Intervention
{
public:
    //////////////////////////////////////////////////////////////////////////
    //  Class constructors and destructors
    //////////////////////////////////////////////////////////////////////////
    
    /////////////////////////////////////
    // Constructor: read intervention data from input files
    Intervention(const char *coverage_File);
    
    
    ////////////////////////////////////////////////////
    // Copy and move constructors

    // Delete unwanted copy constructors
    Intervention(Intervention&) =delete;
    void operator= (Intervention&) =delete;
    // Allow default move constructors
    Intervention(Intervention&&) = default;
    Intervention& operator= (Intervention&&) = default;


    //////////////////////////////////////////////////////////////////////////
    //  Class member functions
    //////////////////////////////////////////////////////////////////////////
    
    /////////////////////////////////////
    // Distribute interventions
    void distribute(double t, Params& theta, Population& POP);
    

private:
    //////////////////////////////////////////////////////////////////////////
    //  Data
    //////////////////////////////////////////////////////////////////////////


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

#endif
