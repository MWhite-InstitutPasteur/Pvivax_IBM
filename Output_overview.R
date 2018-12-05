
OUTPUT <- read.table("Output\\model_output.txt")

N_spec = 3

mosq_species = c("faruati", "punctulatus", "koliensis")

mosq_comp = c( "S_M_far", "E_M_far", "I_M_far",
		   "S_M_pun", "E_M_pun", "I_M_pun",
               "S_M_kol", "E_M_kol", "I_M_kol")

mosq_comp = mosq_comp[1:(N_spec*3)]



colnames(OUTPUT) <- c("time",
                      "S", "I_PCR", "I_LM", "D", "T", "P",
                      mosq_comp,
                      "N_pop", "PvPR_PCR", "PvPR_LM", "Pv_clin", "PvHR",
                      "PvHR_batch", "new_PCR", "new_LM", "new_D", "new_T",
                      "EIR", "LLIN_cov", "IRS_cov", "ACT_treat", "PQ_treat", "PQ_overtreat", "PQ_overtreat_9m")

par(ask=TRUE)

#############################################
#############################################
##          ##                             ##
##  PLOT 1  ##  Human dynamics             ##
##          ##                             ##
#############################################
#############################################

par(mfrow=c(2,3))

plot(x=OUTPUT$time/365, y=OUTPUT$S, type='l',
xlab="time (years)", ylab="number", main="S")

plot(x=OUTPUT$time/365, y=OUTPUT$I_PCR, type='l',
xlab="time (years)", ylab="number", main="I_PCR")

plot(x=OUTPUT$time/365, y=OUTPUT$I_LM, type='l',
xlab="time (years)", ylab="number", main="I_LM")

plot(x=OUTPUT$time/365, y=OUTPUT$D, type='l',
xlab="time (years)", ylab="number", main="D")

plot(x=OUTPUT$time/365, y=OUTPUT$T, type='l',
xlab="time (years)", ylab="number", main="T")

plot(x=OUTPUT$time/365, y=OUTPUT$P, type='l',
xlab="time (years)", ylab="number", main="P")


#############################################
#############################################
##          ##                             ##
##  PLOT 2  ##  Mosquito dynamics          ##
##          ##                             ##
#############################################
#############################################

par(mfrow=c(3,3))

for(g in 1:N_spec)
{
	plot(x=OUTPUT$time/365, y=OUTPUT[,which( colnames(OUTPUT)==mosq_comp[1+(g-1)*3] )], type='l',
	xlab="time (years)", ylab="number", main=mosq_comp[1+(g-1)*3])

	plot(x=OUTPUT$time/365, y=OUTPUT[,which( colnames(OUTPUT)==mosq_comp[2+(g-1)*3] )], type='l',
	xlab="time (years)", ylab="number", main=mosq_comp[2+(g-1)*3])

	plot(x=OUTPUT$time/365, y=OUTPUT[,which( colnames(OUTPUT)==mosq_comp[3+(g-1)*3] )], type='l',
	xlab="time (years)", ylab="number", main=mosq_comp[3+(g-1)*3])
}




#############################################
#############################################
##          ##                             ##
##  PLOT 3  ##  Prevalence                 ##
##          ##                             ##
#############################################
#############################################


par(mfrow=c(2,4))


#############################
## PANEL 1: PCR prevalence

plot(x=OUTPUT$time/365, y=OUTPUT$PvPR_PCR/OUTPUT$N_pop, type='l', ylim=c(0,0.6),
xlab="time (years)", ylab="prevalence", main="Prevalence by PCR")


#############################
## PANEL 2: LM prevalence

plot(x=OUTPUT$time/365, y=OUTPUT$PvPR_LM/OUTPUT$N_pop, type='l', ylim=c(0,0.2),
xlab="time (years)", ylab="prevalence", main="Prevalence by LM")


#############################
## PANEL 3: Clinical incdience

plot(x=OUTPUT$time/365, y=OUTPUT$Pv_clin/OUTPUT$N_pop, type='l', ylim=c(0,0.05),
xlab="time (years)", ylab="prevalence", main="Prevalence of clinical disease")


#############################
## PANEL 4: EIR

plot(x=OUTPUT$time/365, y=365*OUTPUT$EIR, type='l', 
xlab="time (years)", ylab="EIR (ibppy)", main="EIR")


#############################
## PANEL 5: Hypnozoite prevalence

plot(x=OUTPUT$time/365, y=OUTPUT$PvHR/OUTPUT$N_pop, type='l', ylim=c(0,1),
xlab="time (years)", ylab="prevalence", main="Prevalence of hypnozoites")


#############################
## PANEL 6: Hypnozoite broods

plot(x=OUTPUT$time/365, y=OUTPUT$PvHR_batch/OUTPUT$N_pop, type='l', ylim=c(0,max(OUTPUT$PvHR_batch/OUTPUT$N_pop)),
xlab="time (years)", ylab="", main="Mean number of hypnozoite batches per person")


#############################
## PANEL 7: Intervention coverage

INT_COV <- read.table("intervention_coverage.txt")

plot(x=OUTPUT$time/365, y=OUTPUT$LLIN_cov/OUTPUT$N_pop, type='l', ylim=c(0,1),
xlab="time (years)", ylab="coverage", main="Intervention coverage")

points(x=OUTPUT$time/365, y=OUTPUT$IRS_cov/OUTPUT$N_pop, type='l', col="green")





MDA_PQ_t   = OUTPUT$time[which(OUTPUT$PQ_treat > 0.05*OUTPUT$N_pop)]
MDA_PQ_cov = (OUTPUT$PQ_treat/OUTPUT$N_pop)[which(OUTPUT$PQ_treat > 0.05*OUTPUT$N_pop)]

OUTPUT$PQ_treat[which(OUTPUT$PQ_treat > 0.05*OUTPUT$N_pop)] = 0


MDA_ACT_t   = OUTPUT$time[which(OUTPUT$ACT_treat > 0.05*OUTPUT$N_pop)]
MDA_ACT_cov = (OUTPUT$ACT_treat/OUTPUT$N_pop)[which(OUTPUT$ACT_treat > 0.05*OUTPUT$N_pop)]

OUTPUT$ACT_treat[which(OUTPUT$ACT_treat > 0.05*OUTPUT$N_pop)] = 0





N_window = 50

N_tt = length(OUTPUT$ACT_treat)


ACT_smooth = rep(NA, N_tt )
PQ_smooth = rep(NA, N_tt )

for(i in 1:N_tt)
{
	ACT_smooth[i] = 365*mean( OUTPUT$ACT_treat[max(1,i-N_window):min(i+N_window,N_tt)]/OUTPUT$N_pop[max(1,i-N_window):min(i+N_window,N_tt)])
	PQ_smooth[i]  = 365*mean( OUTPUT$PQ_treat[max(1,i-N_window):min(i+N_window,N_tt)]/OUTPUT$N_pop[max(1,i-N_window):min(i+N_window,N_tt)])
}

points(x=OUTPUT$time/365, y=ACT_smooth, type='l', col="blue")

points(x=OUTPUT$time/365, y=PQ_smooth, type='l', col="magenta")


for(j in 1:length(MDA_ACT_t))
{
	points(x=rep(MDA_ACT_t[j]/365,2), y=c(0,MDA_ACT_cov[j]), type='l', col="blue", lwd=2)
}

for(j in 1:length(MDA_PQ_t))
{
	points(x=rep(MDA_PQ_t[j]/365,2), y=c(0,MDA_PQ_cov[j]), type='l', col="magenta", lwd=2)
}



legend(x="topleft", legend=c("LLIN", "IRS", "MDA (BS)", "MDA (BS + PQ)"),
fill=c("black", "green", "blue", "red"),
border=c("black", "green", "blue", "red") )




#############################
## PANEL 8: Primaquine overtreatment

plot(x=OUTPUT$time/365, y=OUTPUT$PQ_overtreat/OUTPUT$N_pop, type='l', ylim=c(0,1),
xlab="time (years)", ylab="", main="Primaquine over-treatment")

points(x=OUTPUT$time/365, y=OUTPUT$PQ_overtreat_9m/OUTPUT$N_pop, type='l', col="red")




