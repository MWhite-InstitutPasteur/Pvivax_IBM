default_intervention_params <- list(
  #######################
  ## LLINs
  ## At each timepoint, x% of the population receive a new net.
  ## If they have an existing net it is replaced.
  LLIN_years = c(), #eg, c(2002, 2005, 2008, 2011, 2015, 2018, 2021, 2024, 2027, 2030, 2033, 2036, 2039, 2042)
  LLIN_cover = c(), #eg, c(0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)

  #######################
  ## IRS
  ## At each timepoint, x% of the population is protected by spraying
  ## Note the model does not have household structure and protection
  ## status is on the individual level
  IRS_years = c(), #eg, c(2002, 2010)
  IRS_cover = c(), #eg, c(0.5, 0.5)

  ###############################
  ## First-line treatment of clinical episodes with 
  ## a blood-stage drug.
  ## Includes information for coverage, efficacy and
  ## duration of prophylaxis
  BS_treat_year_on  =  c(), #eg  c(2000 ,2010)
  BS_treat_year_off =  c(), #eg c(2010, 2018.5)
  BS_treat_cover    =  c(), #eg c(0.5, 0.8)
  BS_treat_BSeff    =  c(), #eg c(0.95, 0.95)
  BS_treat_BSproph  =  c(), #eg c(10, 28)

  ###############################
  ## First-line treatment of clinical episodes with 
  ## a hypnozoiticidal drug (tafenoquine or primaquine).
  ## Includes information for coverage, efficacy and
  ## duration of prophylaxis
  PQ_treat_year_on  = c(), #eg c( 2010, 2021 )
  PQ_treat_year_off = c(), #eg c( 2020, 2025 )
  PQ_treat_cover    = c(), #eg c( 0.9, 0.8 )
  PQ_treat_PQcover  = c(), #eg c( 1, 0.9 )
  PQ_treat_BSeff    = c(), #eg c( 0.95, 0.7 )
  PQ_treat_PQeff    = c(), #eg c( 0.95, 0.8 )
  PQ_treat_BSproph  = c(), #eg c( 14, 60 )
  PQ_treat_PQproph  = c(), #eg c( 14, 60 ) 
  PQ_treat_CYP2D6   = c()  #eg c( 1, 0 )
)

#'@title Generate the intervention file
#'@description to write a file to parameterise the interventions in the model
#'
#'@param overrides, the intervention model params that take precedence over the
#'default parameters
#'@param outpath, the path to write the intervention parameters to
generate_intervention_file <- function(overrides, outpath) {
  params <- default_intervention_params
  for (name in names(overrides)) {
    params[[name]] <- overrides[[name]]
  }

  if(length(params$LLIN_years) != length(params$LLIN_cover)) {
    stop('LLIN vectors are not the same size')
  }

  if(length(params$IRS_years) != length(params$IRS_cover)) {
    stop('IRS vectors are not the same size')
  }

  years <- sort(unique(c(
    params$LLIN_years,
    params$IRS_years,
    params$BS_treat_year_on,
    params$PQ_treat_year_on
  )))

  if (length(years) == 0) {
    INT_cov <- matrix(-1, nrow=1, ncol=24)
    write.table(
      INT_cov,
      file=outpath,
      row.names=FALSE,
      col.names=FALSE
    )
    return()
  }

  INT_cov <- matrix(-1, nrow=length(years), ncol=24)

  colnames(INT_cov) <- c(
    "year_on",
    "LLIN_cov",
    "IRS_cov",
    "MDA_BS_cover", #Unsupported
    "MDA_BS_BSeff", #Unsupported
    "MDA_BS_BSproph", #Unsupported
    "MDA_PQ_cover", #Unsupported
    "MDA_PQ_BSeff", #Unsupported
    "MDA_PQ_PQeff", #Unsupported
    "MDA_PQ_BSproph", #Unsupported
    "MDA_PQ_PQproph", #Unsupported
    "MDA_PQ_CYP2D6", #Unsupported
    "BS_treat_cover",
    "BS_treat_BSeff",
    "BS_treat_BSproph",
    "BS_year_off",
    "PQ_treat_cover",
    "PQ_treat_PQcover",
    "PQ_treat_BSeff",
    "PQ_treat_PQeff",
    "PQ_treat_BSproph",
    "PQ_treat_PQproph",
    "PQ_treat_CYP2D6",
    "PQ_year_off"
  )

  INT_cov[,1] = years

  for(i in 1:nrow(INT_cov))
  {
    if( INT_cov[i,1] %in% params$LLIN_years )
    {
      INT_cov[i,2] = params$LLIN_cover[which(params$LLIN_years==INT_cov[i,1])]
    }

    if( INT_cov[i,1] %in% params$IRS_years )
    {
      INT_cov[i,3] = params$IRS_cover[which(params$IRS_years==INT_cov[i,1])]
    }

    if( INT_cov[i,1] %in% params$BS_treat_year_on )
    {
      INT_cov[i,13] = params$BS_treat_cover[which(params$BS_treat_year_on==INT_cov[i,1])]
      INT_cov[i,14] = params$BS_treat_BSeff[which(params$BS_treat_year_on==INT_cov[i,1])]
      INT_cov[i,15] = params$BS_treat_BSproph[which(params$BS_treat_year_on==INT_cov[i,1])]
      INT_cov[i,16] = params$BS_treat_year_off[which(params$BS_treat_year_on==INT_cov[i,1])]
    }

    if( INT_cov[i,1] %in% params$PQ_treat_year_on )
    {
      INT_cov[i,17] = params$PQ_treat_cover[which(params$PQ_treat_year_on==INT_cov[i,1])]
      INT_cov[i,18] = params$PQ_treat_PQcover[which(params$PQ_treat_year_on==INT_cov[i,1])]
      INT_cov[i,19] = params$PQ_treat_BSeff[which(params$PQ_treat_year_on==INT_cov[i,1])]
      INT_cov[i,20] = params$PQ_treat_PQeff[which(params$PQ_treat_year_on==INT_cov[i,1])]
      INT_cov[i,21] = params$PQ_treat_BSproph[which(params$PQ_treat_year_on==INT_cov[i,1])]
      INT_cov[i,22] = params$PQ_treat_PQproph[which(params$PQ_treat_year_on==INT_cov[i,1])]
      INT_cov[i,23] = params$PQ_treat_CYP2D6[which(params$PQ_treat_year_on==INT_cov[i,1])]
      INT_cov[i,24] = params$PQ_treat_year_off[which(params$PQ_treat_year_on==INT_cov[i,1])]
    }
  }

  INT_cov <- rbind( INT_cov, rep(-1, ncol(INT_cov)) )

  write.table(
    INT_cov,
    file=outpath,
    row.names=FALSE,
    col.names=FALSE
  )
}
