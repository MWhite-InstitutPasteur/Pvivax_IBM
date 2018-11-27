rm(list=ls())

#######################
## LLINs
## At each timepoint, x% of the population receive a new net.
## If they have an existing net it is replaced.

LLIN_years <- c(2008, 2011, 2015, 2018, 2021, 2024, 2027, 2030, 2033, 2036, 2039, 2042)
LLIN_cover <- c( 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)


#######################
## IRS
## At each timepoint, x% of the population is protected by spraying
## Note the model does not have household structure and protection
## status is on the individual level

IRS_years <- c(1995, 2010)
IRS_cover <- c(0.5, 0.5)



###############################
## First-line treatment of clinical episodes with 
## a blood-stage drug.
## Includes information for coverage, efficacy and
## duration of prophylaxis
 
BS_treat_year_on  <-  c() ## c(2000, 2010)      ## Start of first-line treatment regimen
BS_treat_year_off <-  c() ## c(2010, 2018.5)    ## End of first-line treatment regimen
BS_treat_BScover  <-  c() ## c(0.5, 0.8)        ## Coverage of blood-stage drugs (accounts for proportion seeking care)
BS_treat_BSeff    <-  c() ## c(0.95, 0.95)      ## Efficacy of blood-stage drugs 
BS_treat_BSproph  <-  c() ## c(10, 28)          ## Duration of blood-stage prophylaxis


###############################
## First-line treatment of clinical episodes with 
## a hypnozoiticidal drug (tafenoquine or primaquine).
## Includes information for coverage, efficacy and
## duration of prophylaxis

PQ_treat_year_on     <- c( 2015, 2018 )   ## Start of first-line treatment regimen
PQ_treat_year_off    <- c( 2017, 2020 )   ## End of first-line treatment regimen
PQ_treat_BScover     <- c( 0.8, 0.8 )     ## Coverage of blood-stage drugs (accounts for proportion seeking care)
PQ_treat_BSeff       <- c( 1, 1 )         ## Efficacy of blood-stage drugs
PQ_treat_BSproph     <- c( 14, 60 )       ## Duration of blood-stage prophylaxis
PQ_treat_PQavail     <- c( 1, 0.9 )       ## Proportion of people receiving BS drugs where PQ is available
PQ_treat_PQeff       <- c( 1, 0.8 )       ## Efficacy of primaquine 
PQ_treat_PQproph     <- c( 14, 60 )       ## Duration of primaquine prophylaxis (prevents new hypnozoites)
PQ_treat_G6PD_risk   <- c( 1, 1 )         ## Risk in G6PD-deficient individuals
PQ_treat_CYP2D6_risk <- c( 1, 0 )         ## Risk of not working in low CYP2D6 metabolizers
PQ_treat_preg_risk   <- c( 1, 1 )         ## Risk in pregnant women
PQ_treat_low_age     <- c( 0.5, 18 )      ## Lower age limit for treatment (in years) 


###########################
## Mass drug administration with a blood-stage drug.
## Includes information for coverage, efficacy and
## duration of prophylaxis

MDA_BS_years   <- c(2002, 2011)     ## Time of MDA programme
MDA_BS_BScover <- c(0.8, 0.8)       ## Coverage of blood-stage drugs        
MDA_BS_BSeff   <- c(1.0, 1.0)       ## Efficacy of blood-stage drugs
MDA_BS_BSproph <- c(10, 10)         ## Duration of blood-stage prophylaxis 


###########################
## Mass drug administration with a blood-stage drug.
## and a hypnozoiticidal drug (primaquine or tafenoquine).
## Includes information for coverage, efficacy and
## duration of prophylaxis

MDA_PQ_years       <- c( 2022.1 )         ## Time of MDA programme
MDA_PQ_BScover     <- c( 0.8 )            ## Coverage of blood-stage drugs  
MDA_PQ_BSeff       <- c( 1 )              ## Efficacy of blood-stage drugs 
MDA_PQ_BSproph     <- c( 14 )             ## Duration of blood-stage prophylaxis
MDA_PQ_PQavail     <- c( 1 )              ## Availability of primaquine  
MDA_PQ_PQeff       <- c( 0.95 )           ## Efficacy of primaquine
MDA_PQ_PQproph     <- c( 14 )             ## Duration of primaquine prophylaxis (prevents new hypnozoites)
MDA_PQ_G6PD_risk   <- c( 1 )              ## Risk in G6PD-deficient individuals
MDA_PQ_CYP2D6_risk <- c( 1 )              ## Risk of not working in low CYP2D6 metabolizers
MDA_PQ_preg_risk   <- c( 1 )              ## Risk in pregnant women 
MDA_PQ_low_age     <- c( 0.5 )            ## Lower age limit for treatment (in years) 


###########################
## Mass drug administration with a blood-stage drug.
## and a hypnozoiticidal drug (primaquine or tafenoquine).
## Includes information for coverage, efficacy and
## duration of prophylaxis

MSAT_PQ_years       <-  c( 2024.1 )      ## Time of MDA programme
MSAT_PQ_BScover     <-  c( 0.8 )         ## Coverage of blood-stage drugs  
MSAT_PQ_RDT_PCR     <-  c( 2 )           ## What diagnostic tool: 1 for RDT (=LM); 2 for PCR
MSAT_PQ_sens        <-  c( 1 )           ## Sensitivity of diagnostic tool
MSAT_PQ_BSeff       <-  c( 1 )           ## Efficacy of blood-stage drugs 
MSAT_PQ_BSproph     <-  c( 14 )          ## Duration of blood-stage prophylaxis
MSAT_PQ_PQavail     <-  c( 1 )           ## Availability of primaquine  
MSAT_PQ_PQeff       <-  c( 1 )           ## Efficacy of primaquine
MSAT_PQ_PQproph     <-  c( 14 )          ## Duration of primaquine prophylaxis (prevents new hypnozoites) 
MSAT_PQ_G6PD_risk   <-  c( 1 )           ## Risk in G6PD-deficient individuals
MSAT_PQ_CYP2D6_risk <-  c( 1 )           ## Risk of not working in low CYP2D6 metabolizers
MSAT_PQ_preg_risk   <-  c( 1 )           ## Risk in pregnant women 
MSAT_PQ_low_age     <-  c( 1 )           ## Lower age limit for treatment (in years) 



###########################
## Mass drug administration with a blood-stage drug.
## and a hypnozoiticidal drug (primaquine or tafenoquine).
## Includes information for coverage, efficacy and
## duration of prophylaxis

SSAT_PQ_years       <-  c( 2025.1 )      ## Time of MDA programme
SSAT_PQ_BScover     <-  c( 0.8 )         ## Coverage of blood-stage drugs  
SSAT_PQ_sens        <-  c( 1 )           ## Sensitivity of diagnostic tool
SSAT_PQ_spec        <-  c( 1 )           ## Sensitivity of diagnostic tool
SSAT_PQ_BSeff       <-  c( 1 )           ## Efficacy of blood-stage drugs 
SSAT_PQ_BSproph     <-  c( 14 )          ## Duration of blood-stage prophylaxis
SSAT_PQ_PQavail     <-  c( 1 )           ## Availability of primaquine  
SSAT_PQ_PQeff       <-  c( 1 )           ## Efficacy of primaquine
SSAT_PQ_PQproph     <-  c( 14 )          ## Duration of primaquine prophylaxis (prevents new hypnozoites) 
SSAT_PQ_G6PD_risk   <-  c( 1 )           ## Risk in G6PD-deficient individuals
SSAT_PQ_CYP2D6_risk <-  c( 1 )           ## Risk of not working in low CYP2D6 metabolizers
SSAT_PQ_preg_risk   <-  c( 1 )           ## Risk in pregnant women 
SSAT_PQ_low_age     <-  c( 1 )           ## Lower age limit for treatment (in years) 




years <- sort(unique( c(LLIN_years, IRS_years, 
                        BS_treat_year_on, PQ_treat_year_on,
			      MDA_BS_years, MDA_PQ_years, 
                        MSAT_PQ_years, SSAT_PQ_years) ))


INT_cov <- matrix(-1, nrow=length(years), ncol=55)
colnames(INT_cov) <- c("year_on", "LLIN_cov", "IRS_cov",
                       "BS_treat_year_off", "BS_treat_BScover", "BS_treat_BSeff", "BS_treat_BSproph",  
			     "PQ_treat_year_off", "PQ_treat_BScover", "PQ_treat_BSeff", "PQ_treat_BSproph",    
                       "PQ_treat_PQavail", "PQ_treat_PQeff", "PQ_treat_PQproph", "PQ_treat_G6PD_risk", "PQ_treat_CYP2D6_risk", "PQ_treat_preg_risk", "PQ_treat_low_age",    
			     "MDA_BS_BScover", "MDA_BS_BSeff", "MDA_BS_BSproph",
                       "MDA_PQ_BScover", "MDA_PQ_BSeff", "MDA_PQ_BSproph", 
                       "MDA_PQ_PQavail", "MDA_PQ_PQeff", "MDA_PQ_PQproph", "MDA_PQ_G6PD_risk", "MDA_PQ_CYP2D6_risk", "MDA_PQ_preg_risk", "MDA_PQ_low_age",     
			     "MSAT_PQ_BScover", "MSAT_PQ_RDT_PCR", "MSAT_PQ_sens", "MSAT_PQ_BSeff", "MSAT_PQ_BSproph", 
                       "MSAT_PQ_PQavail", "MSAT_PQ_PQeff", "MSAT_PQ_PQproph", "MSAT_PQ_G6PD_risk", "MSAT_PQ_CYP2D6_risk", "MSAT_PQ_preg_risk", "MSAT_PQ_low_age", 
		  	     "SSAT_PQ_BScover", "SSAT_PQ_sens", "SSAT_PQ_spec", "SSAT_PQ_BSeff", "SSAT_PQ_BSproph",
                       "SSAT_PQ_PQavail", "SSAT_PQ_PQeff", "SSAT_PQ_PQproph", "SSAT_PQ_G6PD_risk", "SSAT_PQ_CYP2D6_risk", "SSAT_PQ_preg_risk", "SSAT_PQ_low_age")    

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

	if( INT_cov[i,1] %in% BS_treat_year_on )
	{
		INT_cov[i,4] = BS_treat_year_off[which(BS_treat_year_on==INT_cov[i,1])]
		INT_cov[i,5] = BS_treat_BScover[which(BS_treat_year_on==INT_cov[i,1])]
		INT_cov[i,6] = BS_treat_BSeff[which(BS_treat_year_on==INT_cov[i,1])]
		INT_cov[i,7] = BS_treat_BSproph[which(BS_treat_year_on==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% PQ_treat_year_on )
	{
		INT_cov[i,8]  = PQ_treat_year_off[which(PQ_treat_year_on==INT_cov[i,1])]
		INT_cov[i,9]  = PQ_treat_BScover[which(PQ_treat_year_on==INT_cov[i,1])]
		INT_cov[i,10] = PQ_treat_BSeff[which(PQ_treat_year_on==INT_cov[i,1])]
		INT_cov[i,11] = PQ_treat_BSproph[which(PQ_treat_year_on==INT_cov[i,1])]
		INT_cov[i,12] = PQ_treat_PQavail[which(PQ_treat_year_on==INT_cov[i,1])]
		INT_cov[i,13] = PQ_treat_PQeff[which(PQ_treat_year_on==INT_cov[i,1])]
		INT_cov[i,14] = PQ_treat_PQproph[which(PQ_treat_year_on==INT_cov[i,1])]
		INT_cov[i,15] = PQ_treat_G6PD_risk[which(PQ_treat_year_on==INT_cov[i,1])]
		INT_cov[i,16] = PQ_treat_CYP2D6_risk[which(PQ_treat_year_on==INT_cov[i,1])]
		INT_cov[i,17] = PQ_treat_preg_risk[which(PQ_treat_year_on==INT_cov[i,1])]
		INT_cov[i,18] = PQ_treat_low_age[which(PQ_treat_year_on==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% MDA_BS_years )
	{
		INT_cov[i,19] = MDA_BS_BScover[which(MDA_BS_years==INT_cov[i,1])]
		INT_cov[i,20] = MDA_BS_BSeff[which(MDA_BS_years==INT_cov[i,1])]
		INT_cov[i,21] = MDA_BS_BSproph[which(MDA_BS_years==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% MDA_PQ_years )
	{
		INT_cov[i,22] = MDA_PQ_BScover[which(MDA_PQ_years==INT_cov[i,1])]
		INT_cov[i,23] = MDA_PQ_BSeff[which(MDA_PQ_years==INT_cov[i,1])]
		INT_cov[i,24] = MDA_PQ_BSproph[which(MDA_PQ_years==INT_cov[i,1])]
		INT_cov[i,25] = MDA_PQ_PQavail[which(MDA_PQ_years==INT_cov[i,1])]
		INT_cov[i,26] = MDA_PQ_PQeff[which(MDA_PQ_years==INT_cov[i,1])]
		INT_cov[i,27] = MDA_PQ_PQproph[which(MDA_PQ_years==INT_cov[i,1])]
		INT_cov[i,28] = MDA_PQ_G6PD_risk[which(MDA_PQ_years==INT_cov[i,1])]
		INT_cov[i,29] = MDA_PQ_CYP2D6_risk[which(MDA_PQ_years==INT_cov[i,1])]
		INT_cov[i,30] = MDA_PQ_preg_risk[which(MDA_PQ_years==INT_cov[i,1])]
		INT_cov[i,31] = MDA_PQ_low_age[which(MDA_PQ_years==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% MSAT_PQ_years )
	{
		INT_cov[i,32] = MSAT_PQ_BScover[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,33] = MSAT_PQ_RDT_PCR[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,34] = MSAT_PQ_sens[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,35] = MSAT_PQ_BSeff[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,36] = MSAT_PQ_BSproph[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,37] = MSAT_PQ_PQavail[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,38] = MSAT_PQ_PQeff[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,39] = MSAT_PQ_PQproph[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,40] = MSAT_PQ_G6PD_risk[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,41] = MSAT_PQ_CYP2D6_risk[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,42] = MSAT_PQ_preg_risk[which(MSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,43] = MSAT_PQ_low_age[which(MSAT_PQ_years==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% SSAT_PQ_years )
	{
		INT_cov[i,44] = SSAT_PQ_BScover[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,45] = SSAT_PQ_sens[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,46] = SSAT_PQ_spec[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,47] = SSAT_PQ_BSeff[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,48] = SSAT_PQ_BSproph[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,49] = SSAT_PQ_PQavail[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,50] = SSAT_PQ_PQeff[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,51] = SSAT_PQ_PQproph[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,52] = SSAT_PQ_G6PD_risk[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,53] = SSAT_PQ_CYP2D6_risk[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,54] = SSAT_PQ_preg_risk[which(SSAT_PQ_years==INT_cov[i,1])]
		INT_cov[i,55] = SSAT_PQ_low_age[which(SSAT_PQ_years==INT_cov[i,1])]
	}

}
	


INT_cov <- rbind( INT_cov, rep(-1,ncol(INT_cov)) )
INT_cov <- t(INT_cov)

write.table(INT_cov, file="intervention_coverage.txt", col.names=FALSE)



