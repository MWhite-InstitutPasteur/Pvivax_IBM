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
                        "N_pop_2_10", "PvPR_PCR_2_10", "PvPR_LM_2_10", "Pv_clin_2_10", "PvHR_2_10",
                        "PvHR_batch_2_10", "new_PCR_2_10", "new_LM_2_10", "new_D_2_10", "new_T_2_10",
            "EIR", "LLIN_cov", "IRS_cov", "ACT_treat", "PQ_treat", "pregnant", "A_par", "A_clin")
  output
}
