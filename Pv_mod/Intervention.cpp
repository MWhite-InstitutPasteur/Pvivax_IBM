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

#include "Intervention.hpp"

#include <iostream>
#include <fstream>


/////////////////////////////////////
// Read intervention data from input files
void Intervention::read(const char *coverage_File)
{
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

    for (int j = 0; j<N_cov_rounds; j++)
    {
        //////////////////////////////////////////////////////////////
        // LLINs

        if ((coverage[0][j] > -0.5) && (coverage[1][j] > -0.5))
        {
            LLIN_year.push_back(  coverage[0][j]*365.0 );
            LLIN_cover.push_back( coverage[1][j] );
        }


        //////////////////////////////////////////////////////////////
        // IRS

        if ((coverage[0][j] > -0.5) && (coverage[2][j] > -0.5))
        {
            IRS_year.push_back(  coverage[0][j]*365.0 );
            IRS_cover.push_back( coverage[2][j] );
        }


        //////////////////////////////////////////////////////////////
        // Front-line treatment - blood-stage drugs

        if ((coverage[0][j] > -0.5) && (coverage[4][j] > -0.5))
        {
            BS_treat_year_on.push_back(  coverage[0][j]*365.0 );
            BS_treat_year_off.push_back( coverage[3][j]*365.0 );
            BS_treat_BScover.push_back(  coverage[4][j] );
            BS_treat_BSeff.push_back(    coverage[5][j] );
            BS_treat_BSproph.push_back(  coverage[6][j] );
        }


        //////////////////////////////////////////////////////////////
        // Front-line treatment - primaquine

        if ((coverage[0][j] > -0.5) && (coverage[8][j] > -0.5))
        {
            PQ_treat_year_on.push_back(     coverage[0][j]*365.0 );
            PQ_treat_year_off.push_back(    coverage[7][j]*365.0 );
            PQ_treat_BScover.push_back(     coverage[8][j] );
            PQ_treat_BSeff.push_back(       coverage[9][j] );
            PQ_treat_BSproph.push_back(     coverage[10][j] );
            PQ_treat_PQavail.push_back(     coverage[11][j] );
            PQ_treat_PQeff.push_back(       coverage[12][j] );
            PQ_treat_PQproph.push_back(     coverage[13][j] );
            PQ_treat_G6PD_risk.push_back(   (int)(coverage[14][j]) );
            PQ_treat_CYP2D6_risk.push_back( (int)(coverage[15][j]) );
            PQ_treat_preg_risk.push_back(   (int)(coverage[16][j]) );
            PQ_treat_low_age.push_back(     coverage[17][j] );
        }


        //////////////////////////////////////////////////////////////
        // MDA - blood-stage drugs

        if ((coverage[0][j] > -0.5) && (coverage[18][j] > -0.5))
        {
            MDA_BS_year.push_back(    coverage[0][j]*365.0 );
            MDA_BS_BScover.push_back( coverage[18][j] );
            MDA_BS_BSeff.push_back(   coverage[19][j] );
            MDA_BS_BSproph.push_back( coverage[20][j] );
        }


        //////////////////////////////////////////////////////////////
        // MDA - blood-stage drugs plus primaquine

        if ((coverage[0][j] > -0.5) && (coverage[21][j] > -0.5))
        {
            MDA_PQ_year.push_back(        coverage[0][j]*365.0 );
            MDA_PQ_BScover.push_back(     coverage[21][j] );
            MDA_PQ_BSeff.push_back(       coverage[22][j] );
            MDA_PQ_BSproph.push_back(     coverage[23][j] );
            MDA_PQ_PQavail.push_back(     coverage[24][j] );
            MDA_PQ_PQeff.push_back(       coverage[25][j] );
            MDA_PQ_PQproph.push_back(     coverage[26][j] );
            MDA_PQ_G6PD_risk.push_back(   (int)(coverage[27][j]) );
            MDA_PQ_CYP2D6_risk.push_back( (int)(coverage[28][j]) );
            MDA_PQ_preg_risk.push_back(   (int)(coverage[29][j]) );
            MDA_PQ_low_age.push_back(     coverage[30][j] );
        }


        //////////////////////////////////////////////////////////////
        // MSAT - blood-stage drugs plus primaquine

        if ((coverage[0][j] > -0.5) && (coverage[31][j] > -0.5))
        {
            MSAT_PQ_year.push_back(        coverage[0][j]*365.0 );
            MSAT_PQ_BScover.push_back(     coverage[31][j]);
            MSAT_PQ_RDT_PCR.push_back(     coverage[32][j] );
            MSAT_PQ_sens.push_back(        coverage[33][j] );
            MSAT_PQ_BSeff.push_back(       coverage[34][j] );
            MSAT_PQ_BSproph.push_back(     coverage[35][j] );
            MSAT_PQ_PQavail.push_back(     coverage[36][j] );
            MSAT_PQ_PQeff.push_back(       coverage[37][j] );
            MSAT_PQ_PQproph.push_back(     coverage[38][j] );
            MSAT_PQ_G6PD_risk.push_back(   (int)(coverage[39][j]) );
            MSAT_PQ_CYP2D6_risk.push_back( (int)(coverage[40][j]) );
            MSAT_PQ_preg_risk.push_back(   (int)(coverage[41][j]) );
            MSAT_PQ_low_age.push_back(     coverage[42][j] );
        }

        //////////////////////////////////////////////////////////////
        // SSAT - blood-stage drugs plus primaquine

        if ((coverage[0][j] > -0.5) && (coverage[43][j] > -0.5))
        {
            SSAT_PQ_year.push_back(        coverage[0][j]*365.0 );
            SSAT_PQ_BScover.push_back(     coverage[43][j] );
            SSAT_PQ_sens.push_back(        coverage[44][j] );
            SSAT_PQ_spec.push_back(        coverage[45][j] );
            SSAT_PQ_BSeff.push_back(       coverage[46][j] );
            SSAT_PQ_BSproph.push_back(     coverage[47][j] );
            SSAT_PQ_PQavail.push_back(     coverage[48][j] );
            SSAT_PQ_PQeff.push_back(       coverage[49][j] );
            SSAT_PQ_PQproph.push_back(     coverage[50][j] );
            SSAT_PQ_G6PD_risk.push_back(   (int)(coverage[51][j]) );
            SSAT_PQ_CYP2D6_risk.push_back( (int)(coverage[52][j]) );
            SSAT_PQ_preg_risk.push_back(   (int)(coverage[53][j]) );
            SSAT_PQ_low_age.push_back(     coverage[54][j] );
        }

    }
}
