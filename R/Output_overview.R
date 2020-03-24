N_spec = 3

mosq_comp = c("EL_M_far", "LL_M_far", "P_M_far", "S_M_far", "E_M_far", "I_M_far",
        "EL_M_pun", "LL_M_pun", "P_M_pun", "S_M_pun", "E_M_pun", "I_M_pun",
                "EL_M_kol", "LL_M_kol", "P_M_kol", "S_M_kol", "E_M_kol", "I_M_kol")

mosq_comp = mosq_comp[1:(N_spec*6)]


#'@title present output
#' Adds column names to the model output
#'@param output, the model output
present_output <- function(output) {
  colnames(output) <- c("time",
                        "S", "I_PCR", "I_LM", "D", "T", "P",
                        mosq_comp,
                        "N_pop", "PvPR_PCR", "PvPR_LM", "Pv_clin", "PvHR",
                        "PvHR_batch", "new_PCR", "new_LM", "new_D", "new_T",
                        "N_pop_U5", "PvPR_PCR_U5", "PvPR_LM_U5", "Pv_clin_U5", "PvHR_U5",
                        "PvHR_batch_U5", "new_PCR_U5", "new_LM_U5", "new_D_U5", "new_T_U5",
                        "N_pop_U10", "PvPR_PCR_U10", "PvPR_LM_U10", "Pv_clin_U10", "PvHR_U10",
                        "PvHR_batch_U10", "new_PCR_U10", "new_LM_U10", "new_D_U10", "new_T_U10",
            "EIR", "LLIN_cov", "IRS_cov", "ACT_treat", "PQ_treat", "pregnant", "A_par", "A_clin")
  output
}

#'@title Make some standard plots for the simulation
#'@description takes the simulation output and plots standard graphs. For use in
#'interactive mode only!
plot_outputs <- function(OUTPUT) {
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

  par(mfrow=c(2,3))

  for(g in 1:N_spec)
  {
    plot(x=OUTPUT$time/365, y=OUTPUT[,which( colnames(OUTPUT)==mosq_comp[1+(g-1)*6] )], type='l',
    xlab="time (years)", ylab="number", main=mosq_comp[1+(g-1)*6])

    plot(x=OUTPUT$time/365, y=OUTPUT[,which( colnames(OUTPUT)==mosq_comp[2+(g-1)*6] )], type='l',
    xlab="time (years)", ylab="number", main=mosq_comp[2+(g-1)*6])

    plot(x=OUTPUT$time/365, y=OUTPUT[,which( colnames(OUTPUT)==mosq_comp[3+(g-1)*6] )], type='l',
    xlab="time (years)", ylab="number", main=mosq_comp[3+(g-1)*6])

    plot(x=OUTPUT$time/365, y=OUTPUT[,which( colnames(OUTPUT)==mosq_comp[4+(g-1)*6] )], type='l',
    xlab="time (years)", ylab="number", main=mosq_comp[4+(g-1)*6])

    plot(x=OUTPUT$time/365, y=OUTPUT[,which( colnames(OUTPUT)==mosq_comp[5+(g-1)*6] )], type='l',
    xlab="time (years)", ylab="number", main=mosq_comp[5+(g-1)*6])

    plot(x=OUTPUT$time/365, y=OUTPUT[,which( colnames(OUTPUT)==mosq_comp[6+(g-1)*6] )], type='l',
    xlab="time (years)", ylab="number", main=mosq_comp[6+(g-1)*6])
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

  points(x=OUTPUT$time/365, y=OUTPUT$PvPR_PCR_U10/OUTPUT$N_pop_U10, type='l', col="orange")

  points(x=OUTPUT$time/365, y=OUTPUT$PvPR_PCR_U5/OUTPUT$N_pop_U5, type='l', col="red")


  #############################
  ## PANEL 2: LM prevalence

  plot(x=OUTPUT$time/365, y=OUTPUT$PvPR_LM/OUTPUT$N_pop, type='l', ylim=c(0,0.2),
  xlab="time (years)", ylab="prevalence", main="Prevalence by LM")

  points(x=OUTPUT$time/365, y=OUTPUT$PvPR_LM_U10/OUTPUT$N_pop_U10, type='l', col="orange")

  points(x=OUTPUT$time/365, y=OUTPUT$PvPR_LM_U5/OUTPUT$N_pop_U5, type='l', col="red")


  #############################
  ## PANEL 3: Clinical incdience

  plot(x=OUTPUT$time/365, y=OUTPUT$Pv_clin/OUTPUT$N_pop, type='l', ylim=c(0,0.05),
  xlab="time (years)", ylab="prevalence", main="Prevalence of clinical disease")

  points(x=OUTPUT$time/365, y=OUTPUT$Pv_clin_U10/OUTPUT$N_pop_U10, type='l', col="orange")

  points(x=OUTPUT$time/365, y=OUTPUT$Pv_clin_U5/OUTPUT$N_pop_U5, type='l', col="red")


  #############################
  ## PANEL 4: EIR


  plot(x=OUTPUT$time/365, y=365*OUTPUT$EIR, type='l', 
  xlab="time (years)", ylab="EIR (ibppy)", main="EIR")


  #############################
  ## PANEL 5: Hypnozoite prevalence


  plot(x=OUTPUT$time/365, y=OUTPUT$PvHR/OUTPUT$N_pop, type='l', ylim=c(0,1),
  xlab="time (years)", ylab="prevalence", main="Prevalence of hypnozoites")

  points(x=OUTPUT$time/365, y=OUTPUT$PvHR_U10/OUTPUT$N_pop_U10, type='l', col="orange")

  points(x=OUTPUT$time/365, y=OUTPUT$PvHR_U5/OUTPUT$N_pop_U5, type='l', col="red")


  #############################
  ## PANEL 6: Hypnozoite broods

  plot(x=OUTPUT$time/365, y=OUTPUT$PvHR_batch/OUTPUT$N_pop, type='l', ylim=c(0,max(OUTPUT$PvHR_batch/OUTPUT$N_pop)),
  xlab="time (years)", ylab="", main="Mean number of hypnozoite batches per person")

  points(x=OUTPUT$time/365, y=OUTPUT$PvHR_batch_U10/OUTPUT$N_pop_U10, type='l', col="orange")

  points(x=OUTPUT$time/365, y=OUTPUT$PvHR_batch_U5/OUTPUT$N_pop_U5, type='l', col="red")


  #############################
  ## PANEL 7: Intervention coverage

  # Unsupported graph
  #INT_COV <- read.table("intervention_coverage.txt")

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



  #if( any( INT_COV[,4] > 0 ) )
  #{
  #	MDA_BS_index <- which(INT_COV[,4]>0)
  #	
  #	for(j in 1:length(MDA_BS_index))
  #	{
  #
  #		points( x=rep(INT_COV[MDA_BS_index[j],1],2), y=c(0,INT_COV[MDA_BS_index[j],4]),
  #                    type='l', lwd=2, col="blue")
  #	}
  #}


  #if( any( INT_COV[,5] > 0 ) )
  #{
  #	MDA_PQ_index <- which(INT_COV[,5]>0)
  #	
  #	for(j in 1:length(MDA_PQ_index))
  #	{
  #
  #		points( x=rep(INT_COV[MDA_PQ_index[j],1],2), y=c(0,INT_COV[MDA_PQ_index[j],5]),
  #                    type='l', lwd=2, col="red")
  #	}
  #}

  legend(x="topleft", legend=c("LLIN", "IRS", "MDA (BS)", "MDA (BS + PQ)"),
  fill=c("black", "green", "blue", "red"),
  border=c("black", "green", "blue", "red") )



  #############################
  ## PANEL 9: Immunity


  plot(x=OUTPUT$time/365, y=OUTPUT$A_par, type='l', col="blue",
  ylim=c(0,max(OUTPUT$A_par, OUTPUT$A_clin) ),
  xlab="time (years)", ylab="immunity", main="Anti-parasite & clinical immunity")

  points(x=OUTPUT$time/365, y=OUTPUT$A_clin, type='l', col="orange")
}
