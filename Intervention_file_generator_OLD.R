rm(list=ls())

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

MDA_PQ_years   <- c( 2022.1 )
MDA_PQ_cover   <- c( 0.8 )
MDA_PQ_BSeff   <- c( 1)
MDA_PQ_PQeff   <- c( 0.95 )
MDA_PQ_BSproph <- c( 14 )
MDA_PQ_PQproph <- c( 14 )
MDA_PQ_G6PD    <- c( 1 )
MDA_PQ_CYP2D6  <- c( 1 )
MDA_PQ_preg    <- c( 1 )
MDA_PQ_infant  <- c( 1 )


###########################
## Mass drug administration with a blood-stage drug.
## and a hypnozoiticidal drug (primaquine or tafenoquine).
## Includes information for coverage, efficacy and
## duration of prophylaxis

MSAT_PQ_years   <-  c( 2020.1 )
MSAT_PQ_RDT_PCR <-  c( 2 )   ## 1 for RDT; 2 for PCR
MSAT_PQ_cover   <-  c( 0.8 )
MSAT_PQ_sens    <-  c( 1 )
MSAT_PQ_spec    <-  c( 1 )
MSAT_PQ_BSeff   <-  c( 1 )
MSAT_PQ_PQeff   <-  c( 0.95 )
MSAT_PQ_BSproph <-  c( 14 )
MSAT_PQ_PQproph <-  c( 14 )
MSAT_PQ_G6PD    <-  c( 1 )
MSAT_PQ_CYP2D6  <-  c( 1 )
MSAT_PQ_preg    <-  c( 1 )
MSAT_PQ_infant  <-  c( 1 )


###########################
## Mass drug administration with a blood-stage drug.
## and a hypnozoiticidal drug (primaquine or tafenoquine).
## Includes information for coverage, efficacy and
## duration of prophylaxis

SSAT_PQ_years   <- c( 2018.1 )    ## .1 year ensures mass treatment round after the rainy season peak
SSAT_PQ_cover   <- c( 0.8 )       ## mass treatment coverage
SSAT_PQ_sens    <- c( 0.98 )      ## SSAT sensitivity 
SSAT_PQ_spec    <- c( 0.50 )      ## SSAT specificity
SSAT_PQ_BSeff   <- c( 1 )         ## efficacy of blood-stage druges
SSAT_PQ_PQeff   <- c( 0.95 )      ## efficacy of hypnozoite drugs: CASE 1: PQ = 70%; CASE 2: TQ = 100%; CASE 3: perfect = 100%
SSAT_PQ_BSproph <- c( 14 )        ## duration of blood-stage prophylaxis: CASE 1: PQ = 14 days; CASE 2: TQ = 60 days; CASE 3: perfect = 60 days
SSAT_PQ_PQproph <- c( 14 )        ## duration of liver-stage prophylaxis: CASE 1: PQ = 14 days; CASE 2: TQ = 60 days; CASE 3: perfect = 60 days 
SSAT_PQ_G6PD    <- c( 1 )         ## are there problems in G6PD deficient people: CASE 1: PQ = 1; CASE 2: TQ = 1; CASE 3: perfect = 0
SSAT_PQ_CYP2D6  <- c( 1 )         ## are there drug problems in CYP2D6 low metabolizers: CASE 1: PQ = 1; CASE 2: TQ = 0; CASE 3: perfect = 0
SSAT_PQ_preg    <- c( 1 )         ## are there drug problems in pregnant women: CASE 1: PQ = 1; CASE 2: TQ = 1; CASE 3: perfect = 0
SSAT_PQ_infant  <- c( 1 )         ## are there drug problems in infants <6 months: CASE 1: PQ = 1; CASE 2: TQ = 1; CASE 3: perfect = 0



###############################
## First-line treatment of clinical episodes with 
## a blood-stage drug.
## Includes information for coverage, efficacy and
## duration of prophylaxis

BS_treat_year_on  <-  c() ## c(2000 ,2010)
BS_treat_year_off <-  c() ## c(2010, 2018.5)
BS_treat_cover    <-  c() ## c(0.5, 0.8)
BS_treat_BSeff    <-  c() ## c(0.95, 0.95)
BS_treat_BSproph  <-  c() ## c(10, 28)


###############################
## First-line treatment of clinical episodes with 
## a hypnozoiticidal drug (tafenoquine or primaquine).
## Includes information for coverage, efficacy and
## duration of prophylaxis

PQ_treat_year_on  <- c( 2025, 2028 )
PQ_treat_year_off <- c( 2027, 2040 )
PQ_treat_cover    <- c( 0.9, 0.8 )
PQ_treat_PQcover  <- c( 1, 0.9 )
PQ_treat_BSeff    <- c( 0.95, 0.7 )
PQ_treat_PQeff    <- c( 0.95, 0.8 )
PQ_treat_BSproph  <- c( 14, 60 )
PQ_treat_PQproph  <- c( 14, 60 ) 
PQ_treat_G6PD     <- c( 1, 0 )
PQ_treat_CYP2D6   <- c( 1, 0 )
PQ_treat_preg     <- c( 1, 0 )
PQ_treat_infant   <- c( 1, 0 )


years <- sort(unique( c(LLIN_years, IRS_years, MDA_BS_years, MDA_PQ_years, MSAT_PQ_years, SSAT_PQ_years,
                        BS_treat_year_on, PQ_treat_year_on) ))


INT_cov <- matrix(-1, nrow=length(years), ncol=53)
colnames(INT_cov) <- c("year_on", "LLIN_cov", "IRS_cov",
                       "MDA_BS_cover", "MDA_BS_BSeff",                 "MDA_BS_BSproph", 
                       "MDA_PQ_cover", "MDA_PQ_BSeff", "MDA_PQ_PQeff", "MDA_PQ_BSproph",  "MDA_PQ_PQproph", "MDA_PQ_G6PD", "MDA_PQ_CYP2D6", "MDA_PQ_preg",  "MDA_PQ_infant",
                       "MSAT_PQ_cover", "MSAT_PQ_RDT_PCR", "MSAT_PQ_sens", "MSAT_PQ_spec", "MSAT_PQ_BSeff", "MSAT_PQ_PQeff", "MSAT_PQ_BSproph",  "MSAT_PQ_PQproph", "MSAT_PQ_G6PD", "MSAT_PQ_CYP2D6", "MSAT_PQ_preg", "MSAT_PQ_infant",
                       "SSAT_PQ_cover", "SSAT_PQ_sens", "SSAT_PQ_spec", "SSAT_PQ_BSeff", "SSAT_PQ_PQeff", "SSAT_PQ_BSproph",  "SSAT_PQ_PQproph", "SSAT_PQ_G6PD", "SSAT_PQ_CYP2D6", "SSAT_PQ_preg", "SSAT_PQ_infant",
                       "BS_treat_cover",                     "BS_treat_BSeff",                   "BS_treat_BSproph",                     "BS_year_off",
                       "PQ_treat_cover", "PQ_treat_PQcover", "PQ_treat_BSeff", "PQ_treat_PQeff", "PQ_treat_BSproph", "PQ_treat_PQproph", "PQ_treat_G6PD", "PQ_treat_CYP2D6", "PQ_treat_preg", "PQ_treat_infant", "PQ_year_off") 


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
		INT_cov[i,12] = MDA_PQ_G6PD[which(MDA_PQ_years==INT_cov[i,1])]
		INT_cov[i,13] = MDA_PQ_CYP2D6[which(MDA_PQ_years==INT_cov[i,1])]
		INT_cov[i,14] = MDA_PQ_preg[which(MDA_PQ_years==INT_cov[i,1])]
		INT_cov[i,15] = MDA_PQ_infant[which(MDA_PQ_years==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% MSAT_PQ_years )
	{
		INT_cov[i,16] = MSAT_PQ_cover[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,17] = MSAT_PQ_RDT_PCR[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,18] = MSAT_PQ_sens[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,19] = MSAT_PQ_spec[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,20] = MSAT_PQ_BSeff[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,21] = MSAT_PQ_PQeff[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,22] = MSAT_PQ_BSproph[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,23] = MSAT_PQ_PQproph[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,24] = MSAT_PQ_G6PD[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,25] = MSAT_PQ_CYP2D6[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,26] = MSAT_PQ_preg[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,27] = MSAT_PQ_infant[which(MSAT_PQ_years==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% SSAT_PQ_years )
	{
		INT_cov[i,28] = SSAT_PQ_cover[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,29] = SSAT_PQ_sens[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,30] = SSAT_PQ_spec[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,31] = SSAT_PQ_BSeff[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,32] = SSAT_PQ_PQeff[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,33] = SSAT_PQ_BSproph[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,34] = SSAT_PQ_PQproph[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,35] = SSAT_PQ_G6PD[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,36] = SSAT_PQ_CYP2D6[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,37] = SSAT_PQ_preg[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,38] = SSAT_PQ_infant[which(SSAT_PQ_years==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% BS_treat_year_on )
	{
		INT_cov[i,39] = BS_treat_cover[which(BS_treat_year_on==INT_cov[i,1])]
		INT_cov[i,40] = BS_treat_BSeff[which(BS_treat_year_on==INT_cov[i,1])]
		INT_cov[i,41] = BS_treat_BSproph[which(BS_treat_year_on==INT_cov[i,1])]
		INT_cov[i,42] = BS_treat_year_off[which(BS_treat_year_on==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% PQ_treat_year_on )
	{
		INT_cov[i,43] = PQ_treat_cover[which(PQ_treat_year_on==INT_cov[i,1])]
		INT_cov[i,44] = PQ_treat_PQcover[which(PQ_treat_year_on==INT_cov[i,1])]
		INT_cov[i,45] = PQ_treat_BSeff[which(PQ_treat_year_on==INT_cov[i,1])]
		INT_cov[i,46] = PQ_treat_PQeff[which(PQ_treat_year_on==INT_cov[i,1])]
		INT_cov[i,47] = PQ_treat_BSproph[which(PQ_treat_year_on==INT_cov[i,1])]
		INT_cov[i,48] = PQ_treat_PQproph[which(PQ_treat_year_on==INT_cov[i,1])]
		INT_cov[i,49] = PQ_treat_G6PD[which(PQ_treat_year_on==INT_cov[i,1])]
		INT_cov[i,50] = PQ_treat_CYP2D6[which(PQ_treat_year_on==INT_cov[i,1])]
		INT_cov[i,51] = PQ_treat_preg[which(PQ_treat_year_on==INT_cov[i,1])]
		INT_cov[i,52] = PQ_treat_infant[which(PQ_treat_year_on==INT_cov[i,1])]
		INT_cov[i,53] = PQ_treat_year_off[which(PQ_treat_year_on==INT_cov[i,1])]
	}
}



INT_cov <- rbind( INT_cov, rep(-1,ncol(INT_cov)) )


write.table(INT_cov, file="intervention_SSAT_98_50.txt", row.names=FALSE, col.names=FALSE)



	