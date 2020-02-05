
generate_intervention_file <- function() {
  #######################
  ## LLINs
  ## At each timepoint, x% of the population receive a new net.
  ## If they have an existing net it is replaced.

  LLIN_years <- c() # c(2002, 2005, 2008, 2011, 2015, 2018, 2021, 2024, 2027, 2030, 2033, 2036, 2039, 2042)
  LLIN_cover <- c() # c(0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)


  #######################
  ## IRS
  ## At each timepoint, x% of the population is protected by spraying
  ## Note the model does not have household structure and protection
  ## status is on the individual level

  IRS_years <- c() ##c(2002, 2010)
  IRS_cover <- c() ##c(0.5, 0.5)


  ###########################
  ## Mass drug administration with a blood-stage drug.
  ## Includes information for coverage, efficacy and
  ## duration of prophylaxis

  MDA_BS_years   <- c() ##c(2002, 2011)
  MDA_BS_cover   <- c() ##c(0.8, 0.8)
  MDA_BS_BSeff   <- c() ##c(1.0, 1.0)
  MDA_BS_BSproph <- c() ##c(10, 10)


  ###########################
  ## Mass drug administration with a blood-stage drug.
  ## and a hypnozoiticidal drug (primaquine or tafenoquine).
  ## Includes information for coverage, efficacy and
  ## duration of prophylaxis

  MDA_PQ_years   <- c( 2000, 2020 )
  MDA_PQ_cover   <- c( 0.8, 0.9 )
  MDA_PQ_BSeff   <- c( 1, 1 )
  MDA_PQ_PQeff   <- c( 0.95, 0.95 )
  MDA_PQ_BSproph <- c( 14, 60 )
  MDA_PQ_PQproph <- c( 14, 60 )
  MDA_PQ_CYP2D6  <- c( 1, 0 )


  ###############################
  ## First-line treatment of clinical episodes with 
  ## a blood-stage drug.
  ## Includes information for coverage, efficacy and
  ## duration of prophylaxis

  BS_treat_year_on  <-  c() ## c(2000 ,2010)
  BS_treat_year_off <-  c() ##c(2010, 2018.5)
  BS_treat_cover    <-  c() ##c(0.5, 0.8)
  BS_treat_BSeff    <-  c() ##c(0.95, 0.95)
  BS_treat_BSproph  <-  c() ##c(10, 28)


  ###############################
  ## First-line treatment of clinical episodes with 
  ## a hypnozoiticidal drug (tafenoquine or primaquine).
  ## Includes information for coverage, efficacy and
  ## duration of prophylaxis

  PQ_treat_year_on  <- c( 2010, 2021 )
  PQ_treat_year_off <- c( 2020, 2025 )
  PQ_treat_cover    <- c( 0.9, 0.8 )
  PQ_treat_PQcover  <- c( 1, 0.9 )
  PQ_treat_BSeff    <- c( 0.95, 0.7 )
  PQ_treat_PQeff    <- c( 0.95, 0.8 )
  PQ_treat_BSproph  <- c( 14, 60 )
  PQ_treat_PQproph  <- c( 14, 60 ) 
  PQ_treat_CYP2D6   <- c( 1, 0 )

  years <- sort(unique( c(LLIN_years, IRS_years, MDA_BS_years, MDA_PQ_years,
                          BS_treat_year_on, PQ_treat_year_on) ))


  INT_cov <- matrix(-1, nrow=length(years), ncol=24)
  colnames(INT_cov) <- c("year_on", "LLIN_cov", "IRS_cov",
                         "MDA_BS_cover", "MDA_BS_BSeff",                 "MDA_BS_BSproph", 
                         "MDA_PQ_cover", "MDA_PQ_BSeff", "MDA_PQ_PQeff", "MDA_PQ_BSproph",  "MDA_PQ_PQproph", "MDA_PQ_CYP2D6",
                         "BS_treat_cover",                     "BS_treat_BSeff",                   "BS_treat_BSproph",                     "BS_year_off",
                         "PQ_treat_cover", "PQ_treat_PQcover", "PQ_treat_BSeff", "PQ_treat_PQeff", "PQ_treat_BSproph", "PQ_treat_PQproph", "PQ_treat_CYP2D6", "PQ_year_off") 


  INT_cov[,1] = years

  for(i in 1:nrow(INT_cov))
  {
    if( INT_cov[i,1] %in% LLIN_years )
    {
      INT_cov[i,2] = LLIN_cover[which(LLIN_years==INT_cov[i,1])]
    }

    if( INT_cov[i,1] %in% IRS_years )
    {
      INT_cov[i,3] = IRS_cover[which(IRS_years==INT_cov[i,1])]
    }

    if( INT_cov[i,1] %in% MDA_BS_years )
    {
      INT_cov[i,4] = MDA_BS_cover[which(MDA_BS_years==INT_cov[i,1])]
      INT_cov[i,5] = MDA_BS_BSeff[which(MDA_BS_years==INT_cov[i,1])]
      INT_cov[i,6] = MDA_BS_BSproph[which(MDA_BS_years==INT_cov[i,1])]
    }

    if( INT_cov[i,1] %in% MDA_PQ_years )
    {
      INT_cov[i,7]  = MDA_PQ_cover[which(MDA_PQ_years==INT_cov[i,1])]
      INT_cov[i,8]  = MDA_PQ_BSeff[which(MDA_PQ_years==INT_cov[i,1])]
      INT_cov[i,9]  = MDA_PQ_PQeff[which(MDA_PQ_years==INT_cov[i,1])]
      INT_cov[i,10] = MDA_PQ_BSproph[which(MDA_PQ_years==INT_cov[i,1])]
      INT_cov[i,11] = MDA_PQ_PQproph[which(MDA_PQ_years==INT_cov[i,1])]
      INT_cov[i,12] = MDA_PQ_CYP2D6[which(MDA_PQ_years==INT_cov[i,1])]
    }

    if( INT_cov[i,1] %in% BS_treat_year_on )
    {
      INT_cov[i,13] = BS_treat_cover[which(BS_treat_year_on==INT_cov[i,1])]
      INT_cov[i,14] = BS_treat_BSeff[which(BS_treat_year_on==INT_cov[i,1])]
      INT_cov[i,15] = BS_treat_BSproph[which(BS_treat_year_on==INT_cov[i,1])]
      INT_cov[i,16] = BS_treat_year_off[which(BS_treat_year_on==INT_cov[i,1])]
    }

    if( INT_cov[i,1] %in% PQ_treat_year_on )
    {
      INT_cov[i,17] = PQ_treat_cover[which(PQ_treat_year_on==INT_cov[i,1])]
      INT_cov[i,18] = PQ_treat_PQcover[which(PQ_treat_year_on==INT_cov[i,1])]
      INT_cov[i,19] = PQ_treat_BSeff[which(PQ_treat_year_on==INT_cov[i,1])]
      INT_cov[i,20] = PQ_treat_PQeff[which(PQ_treat_year_on==INT_cov[i,1])]
      INT_cov[i,21] = PQ_treat_BSproph[which(PQ_treat_year_on==INT_cov[i,1])]
      INT_cov[i,22] = PQ_treat_PQproph[which(PQ_treat_year_on==INT_cov[i,1])]
      INT_cov[i,23] = PQ_treat_CYP2D6[which(PQ_treat_year_on==INT_cov[i,1])]
      INT_cov[i,24] = PQ_treat_year_off[which(PQ_treat_year_on==INT_cov[i,1])]
    }
  }



  INT_cov <- rbind( INT_cov, rep(-1,ncol(INT_cov)) )


  write.table(INT_cov, file="intervention_coverage.txt", row.names=FALSE, col.names=FALSE)
}
