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
#include <cmath>
#include <limits>
#include "randlib.h"


////////////////////////////////////////////////////////////
//                                                        //
//  Function declarations                                 //
//                                                        //
////////////////////////////////////////////////////////////

double phi_inv(double pp, double mu, double sigma);


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
