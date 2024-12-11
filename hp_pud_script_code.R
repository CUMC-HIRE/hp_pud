---
  title: "HP Bleeding Ulcers Notebook"
output: html_document
date: "2024-12-09"
authors: "Michael Artin, Josephine Soddano"
---
  
  ```{r Import Libraries}
knitr::opts_chunk$set(echo = TRUE)
library(heemod)
library(data.table)
library(ggplot2)
library(tidyverse)
library(httpgd)
library(diagram)
library(dplyr)
library(rgho)
hgd()
```
```{r Number of Patients in Model}
n_patients <- 10000
```

```{r Transition Probabilities}
param <- define_parameters(
  # Kanotra et al. 2016
  age_init = 65,
  age = age_init + model_time,
  
  # Kanotra et al. 2016
  p_ulcer_male = 0.523,
  
  # Since there is not way to specify population parameters, mortality rate
  # will be weighted for the presumed sex proportion of patients admitted to
  # the hospital for PUD.
  # Change country to USA for US-specific results
  mr = (1 - p_ulcer_male) *get_who_mr(age, sex = "FMLE", country = NULL, local = TRUE) +
    p_ulcer_male *get_who_mr(age, sex = "MLE", country = NULL, local = TRUE),
  
  # For annual incidence of new HP infections we've got Logan and Walker 2001,
  # (0.3%-0.7%), and Luzza et al. 2014 (0.25%, SD 0.10-0.63). The first
  # doesn't cite where they got their statistics from, is for developed
  # countries in general, and seems like it might include both new infection
  # and reinfection. The second is from a prospective cohort study in Italy
  # (developed country), and seems to be novel infection. Given this I think
  # it would be better to use the second source. The first source is higher
  # but that would make sense if it is lumping together reinfection, which
  # presumably is more likely than primary infection.
  p_infect = 0.0025,
  
  # Niv et al. 2008
  # We made the assumption the recrudescence is simply treatment failure,
  # which is captured in non-eradication, so no need to consider a different
  # rate for recrudescence.
  p_reinfect = 0.0145,
  
  # Lau et al. 2011 (95% CI 5.8 - 11.4)
  # This is 30-day mortality, which means that we make the assumption
  # that if they die from the ulcer they are dying within the 1st 30 days,
  # and other causes are from comorbid conditions (e.g. Malmi 2016)
  p_death_bleed = 0.086,
  
  # All sensitivity and specificity for UGIB from Gisbert et al. 2006.,
  # which specifically studied the sensitivity and specificity of these tests
  # after a GI bleed. Note that in the US we are using UBT-14,
  # but Gisbert et al. studied UBT-13. So we will assume the sensitivity and
  # specificity values of the UBT-13 but use the price of the UBT-14 to
  # reflect that we are considering the US.
  ##########################################
  sens_sero = 0.88,
  spec_sero = 0.69,
  
  p_tp_sero = sens_sero,
  p_tn_sero = spec_sero,
  p_fp_sero = 1 - spec_sero,
  p_fn_sero = 1 - sens_sero,
  ##########################################
  sens_stool = 0.87,
  spec_stool = 0.70,
  
  p_tp_stool = sens_stool,
  p_tn_stool = spec_stool,
  p_fp_stool = 1 - spec_stool,
  p_fn_stool = 1 - sens_stool,
  ##########################################
  sens_ubt = 0.93,
  spec_ubt = 0.92,
  
  p_tp_ubt = sens_ubt,
  p_tn_ubt = spec_ubt,
  p_fp_ubt = 1 - spec_ubt,
  p_fn_ubt = 1 - sens_ubt,
  ##########################################
  sens_rut = 0.67,
  spec_rut = 0.93,
  
  p_tp_rut = sens_rut,
  p_tn_rut = spec_rut,
  p_fp_rut = 1 - spec_rut,
  p_fn_rut = 1 - sens_rut,
  ##########################################
  sens_hist = 0.70,
  spec_hist = 0.90,
  p_tp_hist = sens_hist,
  p_tn_hist = spec_hist,
  p_fp_hist = 1 - spec_hist,
  p_fn_hist = 1 - sens_hist,
  ##########################################
  # Gisbert and Calvet 2011, Kowada 2018
  p_erad_first = 0.7,
  # Gisbert 2009
  p_erad_second = 0.81,
  
  # Clinical judgment.
  p_fu_stool = 1,
  
  # Here we need to create transition probabilities from the bleed states to
  # non-death states. In order to accomplish this we need to first consider
  # that if the patient dies in the hospital they won't get any eradication
  # treatment. So once we confirm they don't die, we collapse the rest of the
  # probabilities including whether they get follow-up for stool antigen or
  # UBT, and whether they succeed first-line or second-line eradication
  # therapy. If they fail second-line we assume they remain infected.
  #########################################################
  p_non_erad_bleed_pos_rut = (1 - p_death_bleed) *
    (p_fn_rut + (p_tp_rut * (1 - p_erad_first) * (1 - p_erad_second))),
  p_erad_bleed_pos_rut = 1 - p_non_erad_bleed_pos_rut,
  #########################################################
  p_non_erad_bleed_pos_hist = (1 - p_death_bleed) *
    (p_fn_hist + (p_tp_hist * (1 - p_erad_first) * (1 - p_erad_second))),
  p_erad_bleed_pos_hist = 1 - p_non_erad_bleed_pos_hist,
  #########################################################
  p_non_erad_bleed_pos_stool = (1 - p_death_bleed) *
    ((1 - p_fu_stool) + (p_fu_stool * p_fn_stool) +
       (p_fu_stool * p_tp_stool * (1 - p_erad_first) * (1 - p_erad_second))),
  p_erad_bleed_pos_stool = 1 - p_non_erad_bleed_pos_stool,
  #########################################################
  p_non_erad_bleed_pos_ubt = (1 - p_death_bleed) *
    (p_fn_ubt + (p_tp_ubt * (1 - p_erad_first) * (1 - p_erad_second))),
  p_erad_bleed_pos_ubt = 1 - p_non_erad_bleed_pos_ubt,
  #########################################################
  p_non_erad_bleed_pos_sero = (1 - p_death_bleed) *
    (p_fn_sero + (p_tp_sero * (1 - p_erad_first) * (1 - p_erad_second))),
  p_erad_bleed_pos_sero = 1 - p_non_erad_bleed_pos_sero,
  #########################################################
  
  p_bleed_erad_first_year = 0.029,
  p_bleed_erad_sub_year = 0.0015,
  
  p_bleed_non_erad_first_year = 0.2,
  p_bleed_non_erad_sub_year =  0.058,
  
  # There seems to be some evidence that HP negative ulcers have a higher rate of
  # rebleeding than HP positive with eradication, but they are limited to, for
  # example, non-HP non-NSAID ideopathic ulcers (Hung et al. 2005). The first
  # set of numbers below is from that data, but that fails to capture what's
  # likely the majority of HP negative ulcers secondary to NSAID use. For now
  # will use the rates from HP eradicated. This is unlikely to play a major
  # role in the cost-effectiveness of the diagnostic strategies though.
  
  # p_bleed_still_neg_first_year = 0.134,
  # p_bleed_still_neg_sub_year = 0.076,
  p_bleed_still_neg_first_year = 0.029,
  p_bleed_still_neg_sub_year = 0.0015,
  
  # Use 1-2% annual risk
  # of complications of PUD (Milosavljevic et al. 2011) times the 15% lifetime
  # risk of PUD among HP positive rescaled from a lifetime of 80 years
  p_bleed_pud = 0.015,
  p_pud_pos_lifetime = 0.15,
  p_bleed_pos = p_bleed_pud * rescale_prob(
    p_pud_pos_lifetime, to = 1, from = 80
  ),
  
  # Costs
  # CMS (Histology and RUT assumed to be cost of biopsy)
  c_egd = 312,
  c_egd_biopsy = 401,
  c_biopsy = 89, # c_egd_biopsy - c_egd
  c_stool = 14.38,
  c_sero = 16.85,
  
  # Cost of UBT-14, as UBT-13 is not used in the US.
  c_ubt_13 = 105.47,
  
  # Cost of first and second-line HP eradication therapies
  # (Boklage et al. 2016)
  c_first_line = 119.18,
  c_second_line = 136.62,
  
  # Kanotra et al. 2016
  # Note that for the very first cycle we will not consider the hospital
  # costs, as every patient is assumed admitted to the hospital with a
  # bleeding peptic ulcer, so we don't want to factor that initial sunk cost
  # into the cost-effectiveness of our decisions.
  
  c_hospital = ifelse(model_time == 1, 0, 17193.03),
  
  u_bleed = 0.958,
  u_erad = 1,
  u_pos = 0.9,
  u_neg = 1,
  
  # Transition probabilites #
  
  # We need to calculate the following transition probabilities from the
  # assumption that the patient doesn't die from the bleed. This makes sense
  # because if the patient dies within the year they cannot
  # transition to another state other than death. It is also necessary
  # because without this assumption we could have transition rows that sum
  # to greater than 1.
  
  
  ############### NATURAL HISTORY TRANSITION PROBS #######################
  p_nat_bleed_pos_to_bleed_non_erad = p_bleed_non_erad_first_year,
  p_nat_bleed_neg_to_bleed_still_neg = p_bleed_still_neg_first_year,
  p_nat_bleed_erad_to_bleed_erad = p_bleed_erad_first_year,
  p_nat_bleed_non_erad_to_bleed_non_erad = p_bleed_non_erad_first_year,
  p_nat_bleed_still_neg_to_bleed_still_neg = p_bleed_still_neg_first_year,
  p_nat_erad_to_bleed_erad = p_bleed_erad_sub_year,
  p_nat_erad_to_non_erad = p_reinfect,
  p_nat_non_erad_to_bleed_non_erad = p_bleed_non_erad_sub_year,
  p_nat_pos_to_bleed_pos = p_bleed_pos,
  p_nat_neg_to_bleed_still_neg = p_bleed_still_neg_sub_year,
  p_nat_neg_to_pos = p_infect,
  
  
  ############### RUT TRANSITION PROBS #######################
  p_rut_bleed_pos_to_bleed_erad = p_erad_bleed_pos_rut *
    p_bleed_erad_first_year,
  p_rut_bleed_pos_to_bleed_non_erad = p_non_erad_bleed_pos_rut *
    p_bleed_non_erad_first_year,
  p_rut_bleed_pos_to_non_erad = p_non_erad_bleed_pos_rut,
  p_rut_bleed_neg_to_bleed_still_neg = p_bleed_still_neg_first_year,
  p_rut_bleed_erad_to_bleed_erad = p_bleed_erad_first_year,
  p_rut_bleed_non_erad_to_bleed_erad = p_erad_bleed_pos_rut *
    p_bleed_erad_first_year,
  p_rut_bleed_non_erad_to_bleed_non_erad = p_non_erad_bleed_pos_rut *
    p_bleed_non_erad_first_year,
  p_rut_bleed_non_erad_to_non_erad = p_non_erad_bleed_pos_rut,
  p_rut_bleed_still_neg_to_bleed_still_neg = p_bleed_still_neg_first_year,
  p_rut_erad_to_bleed_erad = p_bleed_erad_sub_year,
  p_rut_erad_to_non_erad = p_reinfect,
  p_rut_non_erad_to_bleed_non_erad = p_bleed_non_erad_sub_year,
  p_rut_pos_to_bleed_pos = p_bleed_pos,
  p_rut_neg_to_bleed_still_neg = p_bleed_still_neg_sub_year,
  p_rut_neg_to_pos = p_infect,
  
  ############### HISTOLOGY TRANSITION PROBS #######################
  p_hist_bleed_pos_to_bleed_erad = p_erad_bleed_pos_hist *
    p_bleed_erad_first_year,
  p_hist_bleed_pos_to_bleed_non_erad = p_non_erad_bleed_pos_hist *
    p_bleed_non_erad_first_year,
  p_hist_bleed_pos_to_non_erad = p_non_erad_bleed_pos_hist,
  p_hist_bleed_neg_to_bleed_still_neg = p_bleed_still_neg_first_year,
  p_hist_bleed_erad_to_bleed_erad = p_bleed_erad_first_year,
  p_hist_bleed_non_erad_to_bleed_erad = p_erad_bleed_pos_hist *
    p_bleed_erad_first_year,
  p_hist_bleed_non_erad_to_bleed_non_erad = p_non_erad_bleed_pos_hist *
    p_bleed_non_erad_first_year,
  p_hist_bleed_non_erad_to_non_erad = p_non_erad_bleed_pos_hist,
  p_hist_bleed_still_neg_to_bleed_still_neg = p_bleed_still_neg_first_year,
  p_hist_erad_to_bleed_erad = p_bleed_erad_sub_year,
  p_hist_erad_to_non_erad = p_reinfect,
  p_hist_non_erad_to_bleed_non_erad = p_bleed_non_erad_sub_year,
  p_hist_pos_to_bleed_pos = p_bleed_pos,
  p_hist_neg_to_bleed_still_neg = p_bleed_still_neg_sub_year,
  p_hist_neg_to_pos = p_infect,
  ############### SAT TRANSITION PROBS #######################
  p_stool_bleed_pos_to_bleed_erad = p_erad_bleed_pos_stool *
    p_bleed_erad_first_year,
  p_stool_bleed_pos_to_bleed_non_erad = p_non_erad_bleed_pos_stool *
    p_bleed_non_erad_first_year,
  p_stool_bleed_pos_to_non_erad = p_non_erad_bleed_pos_stool,
  p_stool_bleed_neg_to_bleed_still_neg = p_bleed_still_neg_first_year,
  p_stool_bleed_erad_to_bleed_erad = p_bleed_erad_first_year,
  p_stool_bleed_non_erad_to_bleed_erad = p_erad_bleed_pos_stool *
    p_bleed_erad_first_year,
  p_stool_bleed_non_erad_to_bleed_non_erad = p_non_erad_bleed_pos_stool *
    p_bleed_non_erad_first_year,
  p_stool_bleed_non_erad_to_non_erad = p_non_erad_bleed_pos_stool,
  p_stool_bleed_still_neg_to_bleed_still_neg = p_bleed_still_neg_first_year,
  p_stool_erad_to_bleed_erad = p_bleed_erad_sub_year,
  p_stool_erad_to_non_erad = p_reinfect,
  p_stool_non_erad_to_bleed_non_erad = p_bleed_non_erad_sub_year,
  p_stool_pos_to_bleed_pos = p_bleed_pos,
  p_stool_neg_to_bleed_still_neg = p_bleed_still_neg_sub_year,
  p_stool_neg_to_pos = p_infect,
  
  ############### UBT TRANSITION PROBS #######################
  p_ubt_bleed_pos_to_bleed_erad = p_erad_bleed_pos_ubt *
    p_bleed_erad_first_year,
  p_ubt_bleed_pos_to_bleed_non_erad = p_non_erad_bleed_pos_ubt *
    p_bleed_non_erad_first_year,
  p_ubt_bleed_pos_to_non_erad = p_non_erad_bleed_pos_ubt,
  p_ubt_bleed_neg_to_bleed_still_neg = p_bleed_still_neg_first_year,
  p_ubt_bleed_erad_to_bleed_erad = p_bleed_erad_first_year,
  p_ubt_bleed_non_erad_to_bleed_erad = p_erad_bleed_pos_ubt *
    p_bleed_erad_first_year,
  p_ubt_bleed_non_erad_to_bleed_non_erad = p_non_erad_bleed_pos_ubt *
    p_bleed_non_erad_first_year,
  p_ubt_bleed_non_erad_to_non_erad = p_non_erad_bleed_pos_ubt,
  p_ubt_bleed_still_neg_to_bleed_still_neg = p_bleed_still_neg_first_year,
  p_ubt_erad_to_bleed_erad = p_bleed_erad_sub_year,
  p_ubt_erad_to_non_erad = p_reinfect,
  p_ubt_non_erad_to_bleed_non_erad = p_bleed_non_erad_sub_year,
  p_ubt_pos_to_bleed_pos = p_bleed_pos,
  p_ubt_neg_to_bleed_still_neg = p_bleed_still_neg_sub_year,
  p_ubt_neg_to_pos = p_infect,
  ############### SEROLOGY TRANSITION PROBS #######################
  p_sero_bleed_pos_to_bleed_erad = p_erad_bleed_pos_sero *
    p_bleed_erad_first_year,
  p_sero_bleed_pos_to_bleed_non_erad = p_non_erad_bleed_pos_sero *
    p_bleed_non_erad_first_year,
  p_sero_bleed_pos_to_non_erad = p_non_erad_bleed_pos_sero,
  p_sero_bleed_neg_to_bleed_still_neg = p_bleed_still_neg_first_year,
  p_sero_bleed_erad_to_bleed_erad = p_bleed_erad_first_year,
  p_sero_bleed_non_erad_to_bleed_erad = p_erad_bleed_pos_sero *
    p_bleed_erad_first_year,
  p_sero_bleed_non_erad_to_bleed_non_erad = p_non_erad_bleed_pos_sero *
    p_bleed_non_erad_first_year,
  p_sero_bleed_non_erad_to_non_erad = p_non_erad_bleed_pos_sero,
  p_sero_bleed_still_neg_to_bleed_still_neg = p_bleed_still_neg_first_year,
  p_sero_erad_to_bleed_erad = p_bleed_erad_sub_year,
  p_sero_erad_to_non_erad = p_reinfect,
  p_sero_non_erad_to_bleed_non_erad = p_bleed_non_erad_sub_year,
  p_sero_pos_to_bleed_pos = p_bleed_pos,
  p_sero_neg_to_bleed_still_neg = p_bleed_still_neg_sub_year,
  p_sero_neg_to_pos = p_infect
  #############################################################
)
```

```{r Natural History Transition Probabilities}
mat_nat <- define_transition(
  state_names = c(
    "state_bleed_pos", "state_bleed_neg", "state_bleed_erad",
    "state_bleed_non_erad", "state_bleed_still_neg", "state_erad",
    "state_non_erad", "state_pos", "state_neg", "state_death"
  ),
  # state_bleed_pos
  0, 0, 0, p_nat_bleed_pos_to_bleed_non_erad, 0, 0, C, 0, 0, p_death_bleed,
  # state_bleed_neg
  0, 0, 0, 0, p_nat_bleed_neg_to_bleed_still_neg, 0, 0, 0, C, p_death_bleed,
  # state_bleed_erad
  0, 0, p_nat_bleed_erad_to_bleed_erad, 0, 0, C, 0, 0, 0, p_death_bleed,
  # state_bleed_non_erad
  0, 0, 0, p_nat_bleed_non_erad_to_bleed_non_erad, 0, 0, C, 0, 0,
  p_death_bleed,
  # state_bleed_still_neg
  0, 0, 0, 0, p_nat_bleed_still_neg_to_bleed_still_neg, 0, 0, 0, C,
  p_death_bleed,
  # state_erad
  0, 0, p_nat_erad_to_bleed_erad, 0, 0, C, p_nat_erad_to_non_erad, 0, 0, mr,
  # state_non_erad
  0, 0, 0, p_nat_non_erad_to_bleed_non_erad, 0, 0, C, 0, 0, mr,
  # state_pos
  p_nat_pos_to_bleed_pos, 0, 0, 0, 0, 0, 0, C, 0, mr,
  # state_neg
  0, 0, 0, 0, p_nat_neg_to_bleed_still_neg, 0, 0, p_nat_neg_to_pos, C, mr,
  # state_death
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1
)
```

```{r RUT Transition Probabilities}
mat_rut <- define_transition(
  state_names = c(
    "state_bleed_pos", "state_bleed_neg", "state_bleed_erad",
    "state_bleed_non_erad", "state_bleed_still_neg", "state_erad",
    "state_non_erad", "state_pos", "state_neg", "state_death"
  ),
  # state_bleed_pos
  0, 0, p_rut_bleed_pos_to_bleed_erad, p_rut_bleed_pos_to_bleed_non_erad, 0,
  C, p_rut_bleed_pos_to_non_erad, 0, 0, p_death_bleed,
  # state_bleed_neg
  0, 0, 0, 0, p_rut_bleed_neg_to_bleed_still_neg, 0, 0, 0, C, p_death_bleed,
  # state_bleed_erad
  0, 0, p_rut_bleed_erad_to_bleed_erad, 0, 0, C, 0, 0, 0, p_death_bleed,
  # state_bleed_non_erad
  0, 0, p_rut_bleed_non_erad_to_bleed_erad,
  p_rut_bleed_non_erad_to_bleed_non_erad, 0, C,
  p_rut_bleed_non_erad_to_non_erad, 0, 0, p_death_bleed,
  # state_bleed_still_neg
  0, 0, 0, 0, p_rut_bleed_still_neg_to_bleed_still_neg, 0, 0, 0, C,
  p_death_bleed,
  # state_erad
  0, 0, p_rut_erad_to_bleed_erad, 0, 0, C, p_rut_erad_to_non_erad, 0, 0, mr,
  # state_non_erad
  0, 0, 0, p_rut_non_erad_to_bleed_non_erad, 0, 0, C, 0, 0, mr,
  # state_pos
  p_rut_pos_to_bleed_pos, 0, 0, 0, 0, 0, 0, C, 0, mr,
  # state_neg
  0, 0, 0, 0, p_rut_neg_to_bleed_still_neg, 0, 0, p_rut_neg_to_pos, C, mr,
  # state_death
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1
)
```

```{r HIST Transition Probabilities }
mat_hist <- define_transition(
  state_names = c(
    "state_bleed_pos", "state_bleed_neg", "state_bleed_erad",
    "state_bleed_non_erad", "state_bleed_still_neg", "state_erad",
    "state_non_erad", "state_pos", "state_neg", "state_death"
  ),
  # state_bleed_pos
  0, 0, p_hist_bleed_pos_to_bleed_erad, p_hist_bleed_pos_to_bleed_non_erad, 0,
  C, p_hist_bleed_pos_to_non_erad, 0, 0, p_death_bleed,
  # state_bleed_neg
  0, 0, 0, 0, p_hist_bleed_neg_to_bleed_still_neg, 0, 0, 0, C, p_death_bleed,
  # state_bleed_erad
  0, 0, p_hist_bleed_erad_to_bleed_erad, 0, 0, C, 0, 0, 0, p_death_bleed,
  # state_bleed_non_erad
  0, 0, p_hist_bleed_non_erad_to_bleed_erad,
  p_hist_bleed_non_erad_to_bleed_non_erad, 0, C,
  p_hist_bleed_non_erad_to_non_erad, 0, 0, p_death_bleed,
  # state_bleed_still_neg
  0, 0, 0, 0, p_hist_bleed_still_neg_to_bleed_still_neg, 0, 0, 0, C,
  p_death_bleed,
  # state_erad
  0, 0, p_hist_erad_to_bleed_erad, 0, 0, C, p_hist_erad_to_non_erad, 0, 0, mr,
  # state_non_erad
  0, 0, 0, p_hist_non_erad_to_bleed_non_erad, 0, 0, C, 0, 0, mr,
  # state_pos
  p_hist_pos_to_bleed_pos, 0, 0, 0, 0, 0, 0, C, 0, mr,
  # state_neg
  0, 0, 0, 0, p_hist_neg_to_bleed_still_neg, 0, 0, p_hist_neg_to_pos, C, mr,
  # state_death
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1
)
```

```{r Stool Transition Probabilities}
mat_stool <- define_transition(
  state_names = c(
    "state_bleed_pos", "state_bleed_neg", "state_bleed_erad",
    "state_bleed_non_erad", "state_bleed_still_neg", "state_erad",
    "state_non_erad", "state_pos", "state_neg", "state_death"
  ),
  # state_bleed_pos
  0, 0, p_stool_bleed_pos_to_bleed_erad, p_stool_bleed_pos_to_bleed_non_erad,
  0, C, p_stool_bleed_pos_to_non_erad, 0, 0, p_death_bleed,
  # state_bleed_neg
  0, 0, 0, 0, p_stool_bleed_neg_to_bleed_still_neg, 0, 0, 0, C, p_death_bleed,
  # state_bleed_erad
  0, 0, p_stool_bleed_erad_to_bleed_erad, 0, 0, C, 0, 0, 0, p_death_bleed,
  # state_bleed_non_erad
  0, 0, p_stool_bleed_non_erad_to_bleed_erad,
  p_stool_bleed_non_erad_to_bleed_non_erad, 0, C,
  p_stool_bleed_non_erad_to_non_erad, 0, 0, p_death_bleed,
  # state_bleed_still_neg
  0, 0, 0, 0, p_stool_bleed_still_neg_to_bleed_still_neg, 0, 0, 0, C,
  p_death_bleed,
  # state_erad
  0, 0, p_stool_erad_to_bleed_erad, 0, 0, C, p_stool_erad_to_non_erad, 0, 0,
  mr,
  # state_non_erad
  0, 0, 0, p_stool_non_erad_to_bleed_non_erad, 0, 0, C, 0, 0, mr,
  # state_pos
  p_stool_pos_to_bleed_pos, 0, 0, 0, 0, 0, 0, C, 0, mr,
  # state_neg
  0, 0, 0, 0, p_stool_neg_to_bleed_still_neg, 0, 0, p_stool_neg_to_pos, C, mr,
  # state_death
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1
)
```

```{r UBT Transition Probabilities}
mat_ubt <- define_transition(
  state_names = c(
    "state_bleed_pos", "state_bleed_neg", "state_bleed_erad",
    "state_bleed_non_erad", "state_bleed_still_neg", "state_erad",
    "state_non_erad", "state_pos", "state_neg", "state_death"
  ),
  # state_bleed_pos
  0, 0, p_ubt_bleed_pos_to_bleed_erad, p_ubt_bleed_pos_to_bleed_non_erad, 0,
  C, p_ubt_bleed_pos_to_non_erad, 0, 0, p_death_bleed,
  # state_bleed_neg
  0, 0, 0, 0, p_ubt_bleed_neg_to_bleed_still_neg, 0, 0, 0, C, p_death_bleed,
  # state_bleed_erad
  0, 0, p_ubt_bleed_erad_to_bleed_erad, 0, 0, C, 0, 0, 0, p_death_bleed,
  # state_bleed_non_erad
  0, 0, p_ubt_bleed_non_erad_to_bleed_erad,
  p_ubt_bleed_non_erad_to_bleed_non_erad, 0, C,
  p_ubt_bleed_non_erad_to_non_erad, 0, 0, p_death_bleed,
  # state_bleed_still_neg
  0, 0, 0, 0, p_ubt_bleed_still_neg_to_bleed_still_neg, 0, 0, 0, C,
  p_death_bleed,
  # state_erad
  0, 0, p_ubt_erad_to_bleed_erad, 0, 0, C, p_ubt_erad_to_non_erad, 0, 0, mr,
  # state_non_erad
  0, 0, 0, p_ubt_non_erad_to_bleed_non_erad, 0, 0, C, 0, 0, mr,
  # state_pos
  p_ubt_pos_to_bleed_pos, 0, 0, 0, 0, 0, 0, C, 0, mr,
  # state_neg
  0, 0, 0, 0, p_ubt_neg_to_bleed_still_neg, 0, 0, p_ubt_neg_to_pos, C, mr,
  # state_death
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1
)
```

```{r Sero Transition Probabilities}
mat_sero <- define_transition(
  state_names = c(
    "state_bleed_pos", "state_bleed_neg", "state_bleed_erad",
    "state_bleed_non_erad", "state_bleed_still_neg", "state_erad",
    "state_non_erad", "state_pos", "state_neg", "state_death"
  ),
  # state_bleed_pos
  0, 0, p_sero_bleed_pos_to_bleed_erad, p_sero_bleed_pos_to_bleed_non_erad, 0,
  C, p_sero_bleed_pos_to_non_erad, 0, 0, p_death_bleed,
  # state_bleed_neg
  0, 0, 0, 0, p_sero_bleed_neg_to_bleed_still_neg, 0, 0, 0, C, p_death_bleed,
  # state_bleed_erad
  0, 0, p_sero_bleed_erad_to_bleed_erad, 0, 0, C, 0, 0, 0, p_death_bleed,
  # state_bleed_non_erad
  0, 0, p_sero_bleed_non_erad_to_bleed_erad,
  p_sero_bleed_non_erad_to_bleed_non_erad, 0, C,
  p_sero_bleed_non_erad_to_non_erad, 0, 0, p_death_bleed,
  # state_bleed_still_neg
  0, 0, 0, 0, p_sero_bleed_still_neg_to_bleed_still_neg, 0, 0, 0, C,
  p_death_bleed,
  # state_erad
  0, 0, p_sero_erad_to_bleed_erad, 0, 0, C, p_sero_erad_to_non_erad, 0, 0, mr,
  # state_non_erad
  0, 0, 0, p_sero_non_erad_to_bleed_non_erad, 0, 0, C, 0, 0, mr,
  # state_pos
  p_sero_pos_to_bleed_pos, 0, 0, 0, 0, 0, 0, C, 0, mr,
  # state_neg
  0, 0, 0, 0, p_sero_neg_to_bleed_still_neg, 0, 0, p_sero_neg_to_pos, C, mr,
  # state_death
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1
)
```

```{r Bleed State Positive}
state_bleed_pos <- define_state(
  c_test_drug = dispatch_strategy(
    strat_nh = (1 - p_death_bleed) * c_biopsy,
    strat_rut = (1 - p_death_bleed) * (c_biopsy +
                                         (p_tp_rut * (c_first_line + ((1 - p_erad_first) * c_second_line)))),
    strat_histology = (1 - p_death_bleed) * (c_biopsy +
                                               (p_tp_hist * (c_first_line + ((1 - p_erad_first) *
                                                                               c_second_line)))),
    strat_stool = (1 - p_death_bleed) * (c_stool +
                                           (p_fu_stool * p_tp_stool *
                                              (c_first_line + ((1 - p_erad_first) * c_second_line)))),
    strat_ubt = (1 - p_death_bleed) * (c_ubt_13 +
                                         (p_tp_ubt * (c_first_line + ((1 - p_erad_first) * c_second_line)))),
    strat_serology = (1 - p_death_bleed) * (c_sero +
                                              (p_tp_sero * (c_first_line + ((1 - p_erad_first) * c_second_line)))),
  ),
  c_test_drug_hospital = discount(c_test_drug + c_hospital, 0.03),
  utility = u_bleed
)
```

```{r Bleed State Negative}
state_bleed_neg <- define_state(
  c_test_drug = dispatch_strategy(
    strat_nh = (1 - p_death_bleed) * c_biopsy,
    strat_rut = (1 - p_death_bleed) * (c_biopsy +
                                         (p_fp_rut * (c_first_line + (p_fp_rut * c_second_line)))),
    strat_histology = (1 - p_death_bleed) * (c_biopsy +
                                               (p_fp_hist * (c_first_line + (p_fp_hist * c_second_line)))),
    strat_stool = (1 - p_death_bleed) * (c_stool +
                                           (p_fu_stool * p_fp_stool *
                                              (c_first_line + (p_fp_stool * c_second_line)))),
    strat_ubt = (1 - p_death_bleed) * (c_ubt_13 +
                                         (p_fp_ubt * (c_first_line + (p_fp_ubt * c_second_line)))),
    strat_serology = (1 - p_death_bleed) * (c_sero +
                                              (p_fp_sero * (c_first_line + (p_fp_sero * c_second_line))))
  ),
  c_test_drug_hospital = discount(c_test_drug + c_hospital, 0.03),
  utility = u_bleed
)
```

```{r Bleed State Erad}
state_bleed_erad <- define_state(
  c_test_drug = dispatch_strategy(
    strat_nh = (1 - p_death_bleed) * c_biopsy,
    strat_rut = (1 - p_death_bleed) * (c_biopsy +
                                         (p_fp_rut * (c_first_line + (p_fp_rut * c_second_line)))),
    strat_histology = (1 - p_death_bleed) * (c_biopsy +
                                               (p_fp_hist * (c_first_line + (p_fp_hist * c_second_line)))),
    strat_stool = (1 - p_death_bleed) * (c_stool +
                                           (p_fu_stool * p_fp_stool *
                                              (c_first_line + (p_fp_stool * c_second_line)))),
    strat_ubt = (1 - p_death_bleed) * (c_ubt_13 +
                                         (p_fp_ubt * (c_first_line + (p_fp_ubt * c_second_line)))),
    strat_serology = (1 - p_death_bleed) * (c_sero +
                                              (p_fp_sero * (c_first_line + (p_fp_sero * c_second_line)))),
  ),
  c_test_drug_hospital = discount(c_test_drug + c_hospital, 0.03),
  utility = u_bleed
)
```

```{r Bleed State Non Erad}
state_bleed_non_erad <- define_state(
  c_test_drug = dispatch_strategy(
    strat_nh = (1 - p_death_bleed) * c_biopsy,
    strat_rut = (1 - p_death_bleed) * (c_biopsy +
                                         (p_tp_rut * (c_first_line + ((1 - p_erad_first) * c_second_line)))),
    strat_histology = (1 - p_death_bleed) * (c_biopsy +
                                               (p_tp_hist * (c_first_line + ((1 - p_erad_first) *
                                                                               c_second_line)))),
    strat_stool = (1 - p_death_bleed) * (c_stool +
                                           (p_fu_stool * p_tp_stool *
                                              (c_first_line + ((1 - p_erad_first) * c_second_line)))),
    strat_ubt = (1 - p_death_bleed) * (c_ubt_13 +
                                         (p_tp_ubt * (c_first_line + ((1 - p_erad_first) * c_second_line)))),
    strat_serology = (1 - p_death_bleed) * (c_sero +
                                              (p_tp_sero * (c_first_line + ((1 - p_erad_first) * c_second_line))))
  ),
  c_test_drug_hospital = discount(c_test_drug + c_hospital, 0.03),
  utility = u_bleed
)
```

```{r Bleed State Still Negative}
state_bleed_still_neg <- define_state(
  c_test_drug = dispatch_strategy(
    strat_nh = (1 - p_death_bleed) * c_biopsy,
    
    strat_rut = (1 - p_death_bleed) * (c_biopsy +
                                         (p_fp_rut * (c_first_line + (p_fp_rut * c_second_line)))),
    
    strat_histology = (1 - p_death_bleed) * (c_biopsy +
                                               (p_fp_hist * (c_first_line + (p_fp_hist * c_second_line)))),
    
    strat_stool = (1 - p_death_bleed) * (c_stool +
                                           (p_fu_stool * p_fp_stool *
                                              (c_first_line + (p_fp_stool * c_second_line)))),
    
    strat_ubt = (1 - p_death_bleed) * (c_ubt_13 +
                                         (p_fp_ubt * (c_first_line + (p_fp_ubt * c_second_line)))),
    
    strat_serology = (1 - p_death_bleed) * (c_sero +
                                              (p_fp_sero * (c_first_line + (p_fp_sero * c_second_line))))
  ),
  c_test_drug_hospital = discount(c_test_drug + c_hospital, 0.03),
  utility = u_bleed
)
```

```{r Define Erad State}
state_erad <- define_state(
  c_test_drug = dispatch_strategy(
    strat_nh = 0,
    strat_rut = 0,
    strat_histology = 0,
    strat_stool = 0,
    strat_ubt = 0,
    strat_serology = 0
  ),
  c_test_drug_hospital = discount(c_test_drug, 0.03),
  utility = u_neg
)
```

```{r Define Non Erad State}
state_non_erad <- define_state(
  c_test_drug = dispatch_strategy(
    strat_nh = 0,
    strat_rut = 0,
    strat_histology = 0,
    strat_stool = 0,
    strat_ubt = 0,
    strat_serology = 0
  ),
  c_test_drug_hospital = discount(c_test_drug, 0.03),
  utility = u_pos
)
```

```{r Define Positive State}
state_pos <- define_state(
  c_test_drug = dispatch_strategy(
    strat_nh = 0,
    strat_rut = 0,
    strat_histology = 0,
    strat_stool = 0,
    strat_ubt = 0,
    strat_serology = 0,
  ),
  c_test_drug_hospital = discount(c_test_drug, 0.03),
  utility = u_pos
)
```

```{r Define Negative State}
state_neg <- define_state(
  c_test_drug = dispatch_strategy(
    strat_nh = 0,
    strat_rut = 0,
    strat_histology = 0,
    strat_stool = 0,
    strat_ubt = 0,
    strat_serology = 0
  ),
  c_test_drug_hospital = discount(c_test_drug, 0.03),
  utility = u_neg
)
```

```{r Define Death State}
state_death <- define_state(
  c_test_drug = 0,
  c_test_drug_hospital = 0,
  utility = 0
)
```

```{r Define NAT Strategy}
strat_nh <- define_strategy(
  transition = mat_nat,
  state_bleed_pos = state_bleed_pos,
  state_bleed_neg = state_bleed_neg,
  state_bleed_erad = state_bleed_erad,
  state_bleed_non_erad = state_bleed_non_erad,
  state_bleed_still_neg = state_bleed_still_neg,
  state_erad = state_erad,
  state_non_erad = state_non_erad,
  state_pos = state_pos,
  state_neg = state_neg,
  state_death = state_death
)
```

```{r Define RUT Strategy}
strat_rut <- define_strategy(
  transition = mat_rut,
  state_bleed_pos = state_bleed_pos,
  state_bleed_neg = state_bleed_neg,
  state_bleed_erad = state_bleed_erad,
  state_bleed_non_erad = state_bleed_non_erad,
  state_bleed_still_neg = state_bleed_still_neg,
  state_erad = state_erad,
  state_non_erad = state_non_erad,
  state_pos = state_pos,
  state_neg = state_neg,
  state_death = state_death
)
```

```{r Define NAT Strategy}
strat_histology <- define_strategy(
  transition = mat_hist,
  state_bleed_pos = state_bleed_pos,
  state_bleed_neg = state_bleed_neg,
  state_bleed_erad = state_bleed_erad,
  state_bleed_non_erad = state_bleed_non_erad,
  state_bleed_still_neg = state_bleed_still_neg,
  state_erad = state_erad,
  state_non_erad = state_non_erad,
  state_pos = state_pos,
  state_neg = state_neg,
  state_death = state_death
)
```

```{r Define Stool Strategy}
strat_stool <- define_strategy(
  transition = mat_stool,
  state_bleed_pos = state_bleed_pos,
  state_bleed_neg = state_bleed_neg,
  state_bleed_erad = state_bleed_erad,
  state_bleed_non_erad = state_bleed_non_erad,
  state_bleed_still_neg = state_bleed_still_neg,
  state_erad = state_erad,
  state_non_erad = state_non_erad,
  state_pos = state_pos,
  state_neg = state_neg,
  state_death = state_death
)
```

```{r Define UBT Strategy}
strat_ubt <- define_strategy(
  transition = mat_ubt,
  state_bleed_pos = state_bleed_pos,
  state_bleed_neg = state_bleed_neg,
  state_bleed_erad = state_bleed_erad,
  state_bleed_non_erad = state_bleed_non_erad,
  state_bleed_still_neg = state_bleed_still_neg,
  state_erad = state_erad,
  state_non_erad = state_non_erad,
  state_pos = state_pos,
  state_neg = state_neg,
  state_death = state_death
)
```

```{r Define Serology Strategy}
strat_serology <- define_strategy(
  transition = mat_sero,
  state_bleed_pos = state_bleed_pos,
  state_bleed_neg = state_bleed_neg,
  state_bleed_erad = state_bleed_erad,
  state_bleed_non_erad = state_bleed_non_erad,
  state_bleed_still_neg = state_bleed_still_neg,
  state_erad = state_erad,
  state_non_erad = state_non_erad,
  state_pos = state_pos,
  state_neg = state_neg,
  state_death = state_death
)
```

```{r Run Model 17% HP Prevelance}
res_mod_17 <- run_model(
  strat_nh = strat_nh,
  strat_rut = strat_rut,
  strat_histology = strat_histology,
  strat_stool = strat_stool,
  strat_serology = strat_serology,
  strat_ubt = strat_ubt,
  parameters = param,
  cycles = 35,
  cost = c_test_drug_hospital,
  effect = utility,
  state_time_limit = 1,
  init = c(
    n_patients * 0.17, n_patients * (1 - 0.17), 0, 0, 0, 0, 0, 0, 0, 0
  ),
  method = "life-table"
)
```

```{r Summary Model 17% HP Prevelance}
summary(res_mod_17, threshold = 100000)
```

```{r Define Model 10% HP Prevelance}
#######################################
# HP Prevalence in Peptic Ulcers: 10% #
#######################################
res_mod_10 <- run_model(
  strat_nh = strat_nh,
  strat_rut = strat_rut,
  strat_histology = strat_histology,
  strat_stool = strat_stool,
  strat_serology = strat_serology,
  strat_ubt = strat_ubt,
  parameters = param,
  cycles = 35,
  cost = c_test_drug_hospital,
  effect = utility,
  state_time_limit = 1,
  init = c(
    n_patients * 0.1, n_patients * (1 - 0.1), 0, 0, 0, 0, 0, 0, 0, 0
  ),
  method = "life-table"
)
```

```{r Summary Model 10% HP}
summary(res_mod_10, threshold = 100000)
```

```{r Define Model 30% HP}
#######################################
# HP Prevalence in Peptic Ulcers: 30% #
#######################################
res_mod_30 <- run_model(
  strat_nh = strat_nh,
  strat_rut = strat_rut,
  strat_histology = strat_histology,
  strat_stool = strat_stool,
  strat_serology = strat_serology,
  strat_ubt = strat_ubt,
  parameters = param,
  cycles = 35,
  cost = c_test_drug_hospital,
  effect = utility,
  state_time_limit = 1,
  init = c(
    n_patients * 0.3, n_patients * (1 - 0.3), 0, 0, 0, 0, 0, 0, 0, 0
  ),
  method = "life-table"
)
```

```{r Summary Model 30% HP}
summary(res_mod_30, threshold = 100000)
```

```{r Define Model 60% HP}
#######################################
# HP Prevalence in Peptic Ulcers: 60% #
#######################################
res_mod_60 <- run_model(
  strat_nh = strat_nh,
  strat_rut = strat_rut,
  strat_histology = strat_histology,
  strat_stool = strat_stool,
  strat_serology = strat_serology,
  strat_ubt = strat_ubt,
  
  parameters = param,
  cycles = 35,
  cost = c_test_drug_hospital,
  effect = utility,
  state_time_limit = 1,
  init = c(
    n_patients * 0.6, n_patients * (1 - 0.6), 0, 0, 0, 0, 0, 0, 0, 0
  ),
  method = "life-table"
)
```

```{r Summary Model 60% HP}
summary(res_mod_60, threshold = 100000)
```

```{r Define Model 90% HP}
#######################################
# HP Prevalence in Peptic Ulcers: 90% #
#######################################
res_mod_90 <- run_model(
  strat_nh = strat_nh,
  strat_rut = strat_rut,
  strat_histology = strat_histology,
  strat_stool = strat_stool,
  strat_serology = strat_serology,
  strat_ubt = strat_ubt,
  parameters = param,
  cycles = 35,
  cost = c_test_drug_hospital,
  effect = utility,
  state_time_limit = 1,
  init = c(
    n_patients * 0.9, n_patients * (1 - 0.9), 0, 0, 0, 0, 0, 0, 0, 0
  ),
  method = "life-table"
)
```

```{r Summary Model 90% HP}
summary(res_mod_90, threshold = 100000)
```

```{r Define Model 17% HP Non Invasive, SAT vs UBT}
# So the same pattern holds at every level of prevalence. It seems reasonable
# to just discard the biopsy methods (invasive) and focus on the non-invasive
# methods, which are the two closest to each other in cost-effectiveness. Just
# going to keep the same structure but remove the other strategies, then setup
# the DSAs.

#####################################################
# HP Prevalence in Peptic Ulcers, Non-Invasive: 17% #
#####################################################
res_mod_noninv_17 <- run_model(
  strat_stool = strat_stool,
  strat_ubt = strat_ubt,
  parameters = param,
  cycles = 35,
  cost = c_test_drug_hospital,
  effect = utility,
  state_time_limit = 1,
  init = c(
    n_patients * 0.17, n_patients * (1 - 0.17), 0, 0, 0, 0, 0, 0, 0, 0
  ),
  method = "life-table"
)
```

```{r Define Model 17% HP Non Invasive, UBT vs Sero}
res_mod_noninv_17seroubt <- run_model(
  strat_serology = strat_serology,
  strat_ubt = strat_ubt,
  parameters = param,
  cycles = 35,
  cost = c_test_drug_hospital,
  effect = utility,
  state_time_limit = 1,
  init = c(
    n_patients * 0.17, n_patients * (1 - 0.17), 0, 0, 0, 0, 0, 0, 0, 0
  ),
  method = "life-table"
)
```

```{r Define Model 17% HP Non Invasive, SAT vs Sero}
res_mod_noninv_17serosat <- run_model(
  strat_stool = strat_stool,
  strat_serology = strat_serology,
  parameters = param,
  cycles = 35,
  cost = c_test_drug_hospital,
  effect = utility,
  state_time_limit = 1,
  init = c(
    n_patients * 0.17, n_patients * (1 - 0.17), 0, 0, 0, 0, 0, 0, 0, 0
  ),
  method = "life-table"
)
```

```{r Summary Model 17 Non Invasive, UBT vs SAT}
summary(res_mod_noninv_17, threshold = 100000)
```

```{r Define Model 10% HP Non Invasive}
#####################################################
# HP Prevalence in Peptic Ulcers, Non-Invasive: 10% #
#####################################################
res_mod_noninv_10 <- run_model(
  strat_stool = strat_stool,
  strat_serology = strat_serology,
  strat_ubt = strat_ubt,
  parameters = param,
  cycles = 35,
  cost = c_test_drug_hospital,
  effect = utility,
  state_time_limit = 1,
  init = c(
    n_patients * 0.1, n_patients * (1 - 0.1), 0, 0, 0, 0, 0, 0, 0, 0
  ),
  method = "life-table"
)
```

```{r Summary Model 10% HP Non Invasive}
summary(res_mod_noninv_10, threshold = 100000)
```

```{r Define Model 30% HP Non Invasive}
#####################################################
# HP Prevalence in Peptic Ulcers, Non-Invasive: 30% #
#####################################################
res_mod_noninv_30 <- run_model(
  strat_stool = strat_stool,
  strat_serology = strat_serology,
  strat_ubt = strat_ubt,
  parameters = param,
  cycles = 35,
  cost = c_test_drug_hospital,
  effect = utility,
  state_time_limit = 1,
  init = c(
    n_patients * 0.3, n_patients * (1 - 0.3), 0, 0, 0, 0, 0, 0, 0, 0
  ),
  method = "life-table"
)
```

```{r Summary Model 30% HP Non Invasive}
summary(res_mod_noninv_30, threshold = 100000)
```

```{r Define Model 60% HP Non Invasive}
#####################################################
# HP Prevalence in Peptic Ulcers, Non-Invasive: 60% #
#####################################################
res_mod_noninv_60 <- run_model(
  strat_stool = strat_stool,
  strat_serology = strat_serology,
  strat_ubt = strat_ubt,
  parameters = param,
  cycles = 35,
  cost = c_test_drug_hospital,
  effect = utility,
  state_time_limit = 1,
  init = c(
    n_patients * 0.6, n_patients * (1 - 0.6), 0, 0, 0, 0, 0, 0, 0, 0
  ),
  method = "life-table"
)
```

```{r Summary Model 60% HP Non Invasive}
summary(res_mod_noninv_60, threshold = 100000)
```

```{r Define Model 90% HP Non Invasive}
#####################################################
# HP Prevalence in Peptic Ulcers, Non-Invasive: 90% #
#####################################################
res_mod_noninv_90 <- run_model(
  strat_stool = strat_stool,
  strat_serology = strat_serology,
  strat_ubt = strat_ubt,
  parameters = param,
  cycles = 35,
  cost = c_test_drug_hospital,
  effect = utility,
  state_time_limit = 1,
  init = c(
    n_patients * 0.9, n_patients * (1 - 0.9), 0, 0, 0, 0, 0, 0, 0, 0
  ),
  method = "life-table"
)
```

```{r Summary Model 90% HP Non Invasive}
summary(res_mod_noninv_90, threshold = 100000)
```

```{r Count Model 17% HP Non Invasive}
# Now we need to get counts of each of the states after the first cycle to
# understand how many additional bleeds due to HP are being prevented by
# choosing one strategy over another.
counts_mod_noninv_10 <- data.table(get_counts(res_mod_noninv_10))
counts_mod_noninv_17 <- data.table(get_counts(res_mod_noninv_17))
```

```{r Average Life Expectancy Model 17% HP}
# Here we can pull the average life expectancy from the cohort.
counts_mod_17 <- data.table(get_counts(res_mod_17))
deaths_mod_17 <- counts_mod_17[
  .strategy_names == "strat_ubt" & state_names == "state_death",
  count
]
sum((c(deaths_mod_17[1], diff(deaths_mod_17)) / 10000) * seq(65, 99))
```

```{r Model 17 Bleeding Erad vs Non Erad }
# Obviously only the bleeding associated with eradication vs. non-eradication
# will matter because the other causes of bleeding are not affected by our
# diagnostic strategies.
model_erad <- counts_mod_17[
  (state_names == "state_bleed_erad" |
     state_names == "state_bleed_non_erad") &
    model_time > 1,
  .(Sum = sum(count)),
  by = .(state_names, .strategy_names)
]
```

```{r Dataframe with Bleeds Avoided and NNT 10% HP}
get_bleed_avoid <- function(counts_table) {
  round(
    (counts_table[
      (state_names == "state_bleed_erad" |
         state_names == "state_bleed_non_erad") &
        model_time > 1,
      .(Sum = sum(count)),
      by = .(state_names, .strategy_names)
    ][c(1, 3), Sum] |> sum()) -
      (counts_table[
        (state_names == "state_bleed_erad" |
           state_names == "state_bleed_non_erad") &
          model_time > 1,
        .(Sum = sum(count)),
        by = .(state_names, .strategy_names)
      ][c(2, 4), Sum] |> sum())
  )
}
bleed_avoid_10 <- get_bleed_avoid(counts_mod_noninv_10)

result_10 <- model_erad %>%
  group_by(.strategy_names) %>%
  summarise(combined_sum = round(sum(Sum)))
result_10$diff_to_strat_nh <- abs(result_10$combined_sum - result_10$combined_sum[result_10$.strategy_names == "strat_nh"])
result_10$diff_to_strat_nh[result_10$.strategy_names == "strat_nh"] <- 0
result_10 <- result_10 %>%
  mutate(nnt = 1 / (diff_to_strat_nh / 10000))

round(1 / (bleed_avoid_10 / 10000))
print(result_10)
```


```{r Additional Bleeds Avoided UBT over Stool 17% HP}
# So to know how many additional bleeding events per 10000 patients we avoid
# by choosing UBT over stool antigen:
bleed_avoid_17 <- get_bleed_avoid(counts_mod_noninv_17)
bleed_avoid_17

```

```{r Model 17% HP Number to be Treated}
# We should get a number needed to treat by inverting the absolute risk
# reduction, which in this case should be bleed_avoid_17 divided by 10000.
nnt_17 <- round(1 / (bleed_avoid_17 / 10000))
nnt_17
```

```{r Number to be Treated All Models}
# In this case we have a NNT of 189, which is very unimpressive. However it does
# reflect what we see here, which is that UBT is only just barely an improvement
# over stool antigen from an effectiveness perspective. Let's get the rest of
# them:
counts_mod_noninv_10 <- data.table(get_counts(res_mod_noninv_10))
counts_mod_noninv_17 <- data.table(get_counts(res_mod_noninv_17))
counts_mod_noninv_30 <- data.table(get_counts(res_mod_noninv_30))
counts_mod_noninv_60 <- data.table(get_counts(res_mod_noninv_60))
counts_mod_noninv_90 <- data.table(get_counts(res_mod_noninv_90))

bleed_avoid_10 <- get_bleed_avoid(counts_mod_noninv_10)
bleed_avoid_17 <- get_bleed_avoid(counts_mod_noninv_17)
bleed_avoid_30 <- get_bleed_avoid(counts_mod_noninv_30)
bleed_avoid_60 <- get_bleed_avoid(counts_mod_noninv_60)
bleed_avoid_90 <- get_bleed_avoid(counts_mod_noninv_90)

nnt_10 <- round(1 / (bleed_avoid_10 / 10000))
nnt_17 <- round(1 / (bleed_avoid_17 / 10000))
nnt_30 <- round(1 / (bleed_avoid_30 / 10000))
nnt_60 <- round(1 / (bleed_avoid_60 / 10000))
nnt_90 <- round(1 / (bleed_avoid_90 / 10000))


```

```{r Prevalance of HP increases effectiveness}
# As the prevalence of HP in peptic ulcers increases the effectiveness (number
# of bleeding events avoided) gets higher for choosing UBT over stool antigen.
# Totally makes sense, because if more bleeds are due to HP then a more
# effective test will have greater impact.
bleed_avoid_10 / 35
bleed_avoid_17 / 35
bleed_avoid_30 / 35
bleed_avoid_60 / 35
bleed_avoid_90 / 35
```

```{r Tests Should be Seen as Equivalent}
# We get all the way down to a NNT of 36 for 90% prevalence. Unrealistic, also
# not extremely impressive, again emphasizing that really from a
# cost-effectiveness stand-point the tests should be seen as equivalent.
nnt_10
nnt_17
nnt_30
nnt_60
nnt_90
```

```{r Plots for State Counts for UBT and Stool}
# Here we can make plots for the state counts for UBT and stool antigen.
# Not very remarkable because due to the similarity in efficacy the graphs will
# look virtually identical.
plot(res_mod_noninv_17, type = "counts")
plot(res_mod_noninv_30, type = "counts")
plot(res_mod_noninv_60, type = "counts")
plot(res_mod_noninv_90, type = "counts")
```

```{r Results Table for Model Results}
# Here is a function to make a results table for the model results.
get_results_table <- function(res_mod) {
  res_mod_sum <- summary(res_mod, threshold = 100000)
  
  results_table <- data.frame(
    res_mod_sum$res_values[".strategy_names"],
    res_mod_sum$res_values[c(1:3)] / 10000,
    res_mod_sum$res_comp[".icer"]
  )
  
  names(results_table) <- c(
    "Strategy", "Costs - Tests + Drugs",
    "Costs Total - Tests + Drugs + Hospital", "QALYs", "ICER"
  )
  
  results_table$Strategy <- plyr::mapvalues(
    results_table$Strategy,
    from = c(
      "strat_nh", "strat_rut", "strat_histology", "strat_stool", "strat_serology", "strat_ubt"
    ),
    to = c("Natural History", "RUT", "Histology", "Stool Antigen",'Serology', "UBT"))
  
  results_table
}
```

```{r Get Results Table Model 17% HP}
write.csv(get_results_table(res_mod_17), "results_table_17.csv", row.names = FALSE)
```

```{r Get Results Table Model 30% HP}
write.csv(get_results_table(res_mod_30), "results_table_30.csv", row.names = FALSE)
```

```{r Get Results Table Model 60% HP}
write.csv(get_results_table(res_mod_60), "results_table_60.csv", row.names = FALSE)
```

```{r Get Results Table Model 90% HP}
write.csv(get_results_table(res_mod_90), "results_table_90.csv", row.names = FALSE)
```

```{r Get Results Table Model 17% HP NonInvasive}
write.csv(get_results_table(res_mod_noninv_17), "results_table_17_non_inv.csv", row.names = FALSE)
```

```{r Get Results Table Model 30% HP Non Invasive}
write.csv(get_results_table(res_mod_noninv_30), "results_table_30_non_inv.csv", row.names = FALSE)
```

```{r Get Results Table Model 60% HP Non Invasive}
write.csv(get_results_table(res_mod_noninv_60), "results_table_60_non_inv.csv", row.names = FALSE)
```

```{r Get Results Table Model 90% HP Non Invasive}
write.csv(get_results_table(res_mod_noninv_90), "results_table_90_non_inv.csv", row.names = FALSE)
```

```{r DSA Section}
###############
# DSA Section #
###############

# We're only going to pit UBT vs stool antigen because they are clearly better
# than the invasive strategies, as well as natural history. Here it's called
# dsa_ubt because the stool antigen is the reference strategy, so UBT will be
# compared to it. Nevertheless we're including parameters of both UBT and stool
# antigen to see the effects since both strategies are reflected in this DSA.
dsa_ubt <- define_dsa(
  c_hospital, 16772.02, 17614.04,
  c_first_line, 59.59, 178.77,
  c_second_line, 68.31, 204.93,
  c_ubt_13, 52.74, 158.21,
  c_stool, 7.19, 21.57,
  sens_ubt, 0.90, 0.95,
  spec_ubt, 0.87, 0.96,
  sens_sero, 0.85, 0.90,
  spec_sero, 0.62, 0.75,
  sens_stool, 0.82, 0.91,
  spec_stool, 0.62, 0.78,
  u_bleed, 0.4, 0.6,
  u_pos, 0.8, 1,
  p_infect, 0.001, 0.0063,
  p_reinfect, 0.00725, 0.02175,
  p_erad_first, 0.5, 0.9,
  p_erad_second, 0.61, 1,
  p_ulcer_male, 0.2615, 0.7845,
  p_bleed_pud, 0.01, 0.02,
  p_pud_pos_lifetime, 0.075, 0.225,
  p_death_bleed, 0.058, 0.114,
  p_bleed_erad_first_year, 0.01, 0.032,
  p_bleed_erad_sub_year, 0.0005, 0.0036,
  p_bleed_non_erad_first_year, 0.14, 0.25,
  p_bleed_non_erad_sub_year, 0.05, 0.067,
  p_bleed_still_neg_first_year, 0.01, 0.032,
  p_bleed_still_neg_sub_year, 0.0005, 0.0036
)
```

```{r Run DSA Model 17% HP NonInvasive, SAT vs UBT}
res_dsa_ubt_17 <- run_dsa(
  model = res_mod_noninv_17,
  dsa = dsa_ubt
)
```
```{r Run DSA Model 17% HP NonInvasive, SAT vs Sero}
res_dsa_ubt_17serosat <- run_dsa(
  model = res_mod_noninv_17serosat,
  dsa = dsa_ubt
)
```
```{r Run DSA Model 17% HP NonInvasive, UBT vs Sero}
res_dsa_ubt_17seroubt <- run_dsa(
  model = res_mod_noninv_17seroubt,
  dsa = dsa_ubt
)
```

```{r Run DSA Model 30% HP NonInvasive, SAT vs UBT}
res_dsa_ubt_30 <- run_dsa(
  model = res_mod_noninv_30,
  dsa = dsa_ubt
)
```

```{r Run DSA Model 60% HP NonInvasive, SAT vs UBT}
res_dsa_ubt_60 <- run_dsa(
  model = res_mod_noninv_60,
  dsa = dsa_ubt
)
```

```{r Run DSA Model 90% HP NonInvasive, SAT vs UBT}
res_dsa_ubt_90 <- run_dsa(
  model = res_mod_noninv_90,
  dsa = dsa_ubt
)
```

```{r Plot Model 17% HP Strategy, ICER, SAT vs UBT}
# From the tornado diagram we can see that the main driver of increased cost-
# effectiveness of UBT over stool antigen is the reduced sensitivity of stool
# antigen in the presence of gastric bleeding. Indeed, at the upper bound of the
# observed range of sensitivities stool antigen would actually be cost-effective
# over UBT. Of course the probability of follow-up for the stool antigen does
# not change the conclusion because UBT is cost-effective even at a follow-up
# rate of 100%, although it would be interesting to see the impact of follow-up
# at a higher sensitivity level of the stool antigen test.
plot_res_dsa_ubt_17 <- plot(res_dsa_ubt_17,
                            strategy = "strat_ubt",
                            result = "icer",
                            type = "difference"
)
plot_res_dsa_ubt_17
```
```{r Plot Model 17% HP Strategy, ICER, SAT vs Sero}

plot_res_dsa_ubt_17serosat <- plot(res_dsa_ubt_17serosat,
                                  strategy = "strat_serology",
                                  result = "icer",
                                  type = "difference"
)
plot_res_dsa_ubt_17serosat
```
```{r Plot Model 17% HP Strategy, ICER, UBT vs Sero}
plot_res_dsa_ubt_17seroubt <- plot(res_dsa_ubt_17seroubt,
                                  strategy = "strat_ubt",
                                  result = "icer",
                                  type = "difference"
)
plot_res_dsa_ubt_17seroubt
```

```{r Top 10 Variables from Tornado Diagram, UBT vs SAT}
# Get the top 10 variables from the tornado diagram.
dsa_top_ten <- as.character(data.table(
  icer_diff = abs(
    plot_res_dsa_ubt_17$data$.icer[
      seq(nrow(plot_res_dsa_ubt_17$data))[c(TRUE, FALSE)]
    ] - plot_res_dsa_ubt_17$data$.icer[
      seq(nrow(plot_res_dsa_ubt_17$data))[c(FALSE, TRUE)]
    ]
  ),
  variable = plot_res_dsa_ubt_17$data$.par_names[c(TRUE, FALSE)]
)[order(-icer_diff)][1:10, variable])

dsa_top_ten
```
```{r Top 10 Variables from Tornado Diagram, SAT vs Sero}
# Get the top 10 variables from the tornado diagram.
dsa_top_tenserosat <- as.character(data.table(
  icer_diff = abs(
    plot_res_dsa_ubt_17serosat$data$.icer[
      seq(nrow(plot_res_dsa_ubt_17serosat$data))[c(TRUE, FALSE)]
    ] - plot_res_dsa_ubt_17serosat$data$.icer[
      seq(nrow(plot_res_dsa_ubt_17serosat$data))[c(FALSE, TRUE)]
    ]
  ),
  variable = plot_res_dsa_ubt_17serosat$data$.par_names[c(TRUE, FALSE)]
)[order(-icer_diff)][1:10, variable])

dsa_top_tenserosat
```


```{r Top 10 Variables from Tornado Diagram, UBT vs Sero}
# Get the top 10 variables from the tornado diagram.
dsa_top_tenseroubt <- as.character(data.table(
  icer_diff = abs(
    plot_res_dsa_ubt_17seroubt$data$.icer[
      seq(nrow(plot_res_dsa_ubt_17seroubt$data))[c(TRUE, FALSE)]
    ] - plot_res_dsa_ubt_17seroubt$data$.icer[
      seq(nrow(plot_res_dsa_ubt_17seroubt$data))[c(FALSE, TRUE)]
    ]
  ),
  variable = plot_res_dsa_ubt_17seroubt$data$.par_names[c(TRUE, FALSE)]
)[order(-icer_diff)][1:10, variable])

dsa_top_tenseroubt
```

```{r Trim Top 10 Plot, UBT vs SAT}
# Trim our plot to the top ten most impactful variables.
plot_res_dsa_ubt_17$data <- plot_res_dsa_ubt_17$data[
  plot_res_dsa_ubt_17$data$.par_names %in% dsa_top_ten,
]
plot_res_dsa_ubt_17
```

```{r Trim Top 10 Plot, SAT vs Sero}
# Trim our plot to the top ten most impactful variables.
plot_res_dsa_ubt_17serosat$data <- plot_res_dsa_ubt_17serosat$data[
  plot_res_dsa_ubt_17serosat$data$.par_names %in% dsa_top_tenserosat,
]
plot_res_dsa_ubt_17serosat
```
```{r Trim Top 10 Plot UBT vs Sero}
# Trim our plot to the top ten most impactful variables.
plot_res_dsa_ubt_17seroubt$data <- plot_res_dsa_ubt_17seroubt$data[
  plot_res_dsa_ubt_17seroubt$data$.par_names %in% dsa_top_tenseroubt,
]
plot_res_dsa_ubt_17seroubt
```

```{r Prettify Top 10 Plot UBT vs SAT}
# Time to prettify the plot.
plot_res_dsa_ubt_17$data$.par_names <- plot_res_dsa_ubt_17$data$.par_names |>
  dplyr::recode(
    c_first_line = "Cost of first-line eradication therapy",
    c_hospital = "Cost of hospitalization, PUD",
    c_second_line = "Cost of second-line eradication therapy",
    c_stool = "Cost of SAT",
    c_ubt_13 = "Cost of UBT",
    p_bleed_erad_first_year =
      "Risk of rebleeding, HP eradication, first year",
    p_bleed_erad_sub_year =
      "Risk of rebleeding, HP eradication, subsequent years",
    p_bleed_non_erad_first_year =
      "Risk of rebleeding, HP noneradication, first year",
    p_bleed_non_erad_sub_year =
      "Risk of rebleeding, HP noneradication, subsequent years",
    p_bleed_still_neg_first_year =
      "Risk of rebleeding, HP negative, first year",
    p_bleed_still_neg_sub_year =
      "Risk of rebleeding, HP negative, subsequent years",
    p_bleed_pud = "Risk of bleeding, PUD",
    p_death_bleed = "30-day mortality after peptic ulcer bleed",
    p_erad_first = "First-line eradication rate",
    p_erad_second = "Second-line eradication rate",
    p_fu_stool = "Rate of SAT follow-up",
    p_infect = "HP infection rate",
    p_pud_pos_lifetime = "Lifetime risk of PUD, HP positive",
    p_ulcer_male = "Percentage of male PUD hospitalizations",
    sens_stool = "SAT sensitivity",
    spec_stool = "SAT specificity",
    sens_ubt = "UBT sensitivity",
    spec_ubt = "UBT specificity",
    sens_sero = "Serology sensitivity",
    spec_sero = "Serology specificity",
    u_bleed = "Utility of peptic ulcer bleed",
    u_pos = "Utility of HP positive state",
    p_reinfect = "HP reinfection rate"
  )
plot_res_dsa_ubt_17$data$.strategy_names <- ""

plot_res_dsa_ubt_17$labels$y <- "Parameter"

list(plot_res_dsa_ubt_17$labels)

plot_res_dsa_ubt_17 +
  ggplot2::ggtitle(
    "Figure 1. One-way sensitivity analyses: UBT vs. SAT at prevalence of 17%"
  ) +
  ggplot2::theme(plot.title.position = "plot")

plot_res_dsa_ubt_17$data$.hjust[
  plot_res_dsa_ubt_17$data$.hjust == 0
] <- -0.1

plot_res_dsa_ubt_17$data$.hjust[
  plot_res_dsa_ubt_17$data$.hjust == 1
] <- 1.1
```
```{r Prettify Top 10 Plot SAT vs Sero}
# Time to prettify the plot.
plot_res_dsa_ubt_17serosat$data$.par_names <- plot_res_dsa_ubt_17serosat$data$.par_names |>
  dplyr::recode(
    c_first_line = "Cost of first-line eradication therapy",
    c_hospital = "Cost of hospitalization, PUD",
    c_second_line = "Cost of second-line eradication therapy",
    c_stool = "Cost of SAT",
    c_ubt_13 = "Cost of UBT",
    p_bleed_erad_first_year =
      "Risk of rebleeding, HP eradication, first year",
    p_bleed_erad_sub_year =
      "Risk of rebleeding, HP eradication, subsequent years",
    p_bleed_non_erad_first_year =
      "Risk of rebleeding, HP noneradication, first year",
    p_bleed_non_erad_sub_year =
      "Risk of rebleeding, HP noneradication, subsequent years",
    p_bleed_still_neg_first_year =
      "Risk of rebleeding, HP negative, first year",
    p_bleed_still_neg_sub_year =
      "Risk of rebleeding, HP negative, subsequent years",
    p_bleed_pud = "Risk of bleeding, PUD",
    p_death_bleed = "30-day mortality after peptic ulcer bleed",
    p_erad_first = "First-line eradication rate",
    p_erad_second = "Second-line eradication rate",
    p_fu_stool = "Rate of SAT follow-up",
    p_infect = "HP infection rate",
    p_pud_pos_lifetime = "Lifetime risk of PUD, HP positive",
    p_ulcer_male = "Percentage of male PUD hospitalizations",
    sens_stool = "SAT sensitivity",
    spec_stool = "SAT specificity",
    sens_ubt = "UBT sensitivity",
    spec_ubt = "UBT specificity",
    sens_sero = "Serology sensitivity",
    spec_sero = "Serology specificity",
    u_bleed = "Utility of peptic ulcer bleed",
    u_pos = "Utility of HP positive state",
    p_reinfect = "HP reinfection rate"
  )
plot_res_dsa_ubt_17serosat$data$.strategy_names <- ""

plot_res_dsa_ubt_17serosat$labels$y <- "Parameter"

list(plot_res_dsa_ubt_17serosat$labels)

plot_res_dsa_ubt_17serosat +
  ggplot2::ggtitle(
    "Figure 1. One-way sensitivity analyses: IgG vs. SAT at prevalence of 17%"
  ) +
  ggplot2::theme(plot.title.position = "plot")

plot_res_dsa_ubt_17serosat$data$.hjust[
  plot_res_dsa_ubt_17serosat$data$.hjust == 0
] <- -0.1

plot_res_dsa_ubt_17serosat$data$.hjust[
  plot_res_dsa_ubt_17serosat$data$.hjust == 1
] <- 1.1
```
```{r Prettify Top 10 Plot UBT vs Sero}
# Time to prettify the plot.
plot_res_dsa_ubt_17seroubt$data$.par_names <- plot_res_dsa_ubt_17seroubt$data$.par_names |>
  dplyr::recode(
    c_first_line = "Cost of first-line eradication therapy",
    c_hospital = "Cost of hospitalization, PUD",
    c_second_line = "Cost of second-line eradication therapy",
    c_stool = "Cost of SAT",
    c_ubt_13 = "Cost of UBT",
    p_bleed_erad_first_year =
      "Risk of rebleeding, HP eradication, first year",
    p_bleed_erad_sub_year =
      "Risk of rebleeding, HP eradication, subsequent years",
    p_bleed_non_erad_first_year =
      "Risk of rebleeding, HP noneradication, first year",
    p_bleed_non_erad_sub_year =
      "Risk of rebleeding, HP noneradication, subsequent years",
    p_bleed_still_neg_first_year =
      "Risk of rebleeding, HP negative, first year",
    p_bleed_still_neg_sub_year =
      "Risk of rebleeding, HP negative, subsequent years",
    p_bleed_pud = "Risk of bleeding, PUD",
    p_death_bleed = "30-day mortality after peptic ulcer bleed",
    p_erad_first = "First-line eradication rate",
    p_erad_second = "Second-line eradication rate",
    p_fu_stool = "Rate of SAT follow-up",
    p_infect = "HP infection rate",
    p_pud_pos_lifetime = "Lifetime risk of PUD, HP positive",
    p_ulcer_male = "Percentage of male PUD hospitalizations",
    sens_stool = "SAT sensitivity",
    spec_stool = "SAT specificity",
    sens_ubt = "UBT sensitivity",
    spec_ubt = "UBT specificity",
    sens_sero = "Serology sensitivity",
    spec_sero = "Serology specificity",
    u_bleed = "Utility of peptic ulcer bleed",
    u_pos = "Utility of HP positive state",
    p_reinfect = "HP reinfection rate"
  )
plot_res_dsa_ubt_17seroubt$data$.strategy_names <- ""

plot_res_dsa_ubt_17seroubt$labels$y <- "Parameter"

list(plot_res_dsa_ubt_17seroubt$labels)

plot_res_dsa_ubt_17seroubt +
  ggplot2::ggtitle(
    "Figure 1. One-way sensitivity analyses: IgG vs. UBT at prevalence of 17%"
  ) +
  ggplot2::theme(plot.title.position = "plot")

plot_res_dsa_ubt_17seroubt$data$.hjust[
  plot_res_dsa_ubt_17serosat$data$.hjust == 0
] <- -0.1

plot_res_dsa_ubt_17seroubt$data$.hjust[
  plot_res_dsa_ubt_17seroubt$data$.hjust == 1
] <- 1.1
```

```{r UBT remains cost effective across range of stool sensitivity, SAT vs UBT}
# Here it appears that at high prevalence levels of HP in bleeding peptic ulcers
# the UBT will remain cost-effective across the range of sensitivity of the
# stool antigen test.
plot(res_dsa_ubt_60,
     strategy = "strat_ubt",
     result = "icer",
     type = "difference"
)
plot(res_dsa_ubt_90,
     strategy = "strat_ubt",
     result = "icer",
     type = "difference"
)
```

```{r New parameters to test stool follow up}
# To answer the question about whether stool antigen follow-up is relevant at
# the higher end of the stool antigen sensitivity, we can rerun the model with
# new values for parameters. In this case the upper limit of the sensitivity for
# the stool antigen test: 0.91
new_param <- param
new_param$sens_stool <- 0.91
```

```{r Stool less costly, more effective at higher sensitivity, UBT vs SAT 17% HP}
# We can see now that at a higher sensitivity stool antigen is less costly but
# still less effective, with an ICER that is still under our WTP of 100000.
res_mod_noninv_17_new <- run_model(
  strat_stool = strat_stool,
  strat_ubt = strat_ubt,
  parameters = new_param,
  cycles = 35,
  cost = c_test_drug_hospital,
  effect = utility,
  state_time_limit = 1,
  init = c(
    n_patients * 0.17, n_patients * (1 - 0.17), 0, 0, 0, 0, 0, 0, 0, 0
  ),
  method = "life-table"
)
res_mod_noninv_17_new
```

```{r Summary Model 17% HP new parameters, UBT vs SAT}
summary(res_mod_noninv_17_new, threshold = 100000)
```

```{r Model 30% HP New Parameters, UBT vs SAT}
res_mod_noninv_30_new <- run_model(
  strat_stool = strat_stool,
  strat_ubt = strat_ubt,
  parameters = new_param,
  cycles = 35,
  cost = c_test_drug_hospital,
  effect = utility,
  state_time_limit = 1,
  init = c(
    n_patients * 0.3, n_patients * (1 - 0.3), 0, 0, 0, 0, 0, 0, 0, 0
  ),
  method = "life-table"
)
```

```{r Summary Model 30% HP New Parameters}
summary(res_mod_noninv_30_new, threshold = 100000)
```

```{r Model 60% HP New Parameters}
# And at higher prevalence even improved sensitivity does not prevent stool
# antigen from being dominated.
res_mod_noninv_60_new <- run_model(
  strat_stool = strat_stool,
  strat_ubt = strat_ubt,
  parameters = new_param,
  cycles = 35,
  cost = c_test_drug_hospital,
  effect = utility,
  state_time_limit = 1,
  init = c(
    n_patients * 0.6, n_patients * (1 - 0.6), 0, 0, 0, 0, 0, 0, 0, 0
  ),
  method = "life-table"
)
```

```{r Summary Model 60% HP New Parameters}
summary(res_mod_noninv_60_new, threshold = 100000)
```

```{r Model 90% HP New Parameters}
res_mod_noninv_90_new <- run_model(
  strat_stool = strat_stool,
  strat_ubt = strat_ubt,
  parameters = new_param,
  cycles = 35,
  cost = c_test_drug_hospital,
  effect = utility,
  state_time_limit = 1,
  init = c(
    n_patients * 0.9, n_patients * (1 - 0.9), 0, 0, 0, 0, 0, 0, 0, 0
  ),
  method = "life-table"
)
```

```{r Summary Model 90% HP New Parameters}
summary(res_mod_noninv_90_new, threshold = 100000)
```

```{r Remove Stool Sensitivity from DSA}
# Let's look at the DSA again. We'll remove sensitivity of stool from the ranges
# because we want to hold that constant at its upper bound.
dsa_ubt_new <- define_dsa(
  c_hospital, 16772.02, 17614.04,
  c_first_line, 95.34, 143.02,
  c_second_line, 109.30, 163.94,
  c_ubt_13, 84.38, 126.56,
  sens_ubt, 0.90, 0.95,
  spec_ubt, 0.87, 0.96,
  
  sens_sero, 0.85, 0.90,
  spec_sero, 0.62, 0.75,
  spec_stool, 0.62, 0.78,
  p_fu_stool, 0.9, 1,
  u_bleed, 0.4, 0.6,
  u_pos, 0.8, 1,
  p_reinfect, 0.0116, 0.0174,
  p_erad_first, 0.5, 0.9,
  p_erad_second, 0.61, 1,
  p_ulcer_male, 0.4184, 0.6276
)
```

```{r Model 17% HP Stool Removed DSA}
res_dsa_ubt_17_new <- run_dsa(
  model = res_mod_noninv_17_new,
  dsa = dsa_ubt_new
)
#res_dsa_ubt_17_new
```

```{r Model 30% HP Stool Removed DSA}
res_dsa_ubt_30_new <- run_dsa(
  model = res_mod_noninv_30_new,
  dsa = dsa_ubt_new
)
```

```{r Model 60% HP Stool Removed DSA}
res_dsa_ubt_60_new <- run_dsa(
  model = res_mod_noninv_60_new,
  dsa = dsa_ubt_new
)
```

```{r Model 90% HP Stool Removed DSA}
res_dsa_ubt_90_new <- run_dsa(
  model = res_mod_noninv_90_new,
  dsa = dsa_ubt_new
)
```

```{r View Model 17% HP Stool Removoved DSA}
#res_dsa_ubt_17_new
```

```{r Model 17% HP New Plot Stool Removed }
# So after all that we can see that indeed at 90% follow-up after stool antigen
# testing with a prevalence of 17% we can get the ICER to go negative again.
# However this is not relevant as with a WTP of 100000 UBT is still cost
# effective over stool antigen at any rate of stool antigen follow-up.
plot(res_dsa_ubt_17_new,
     strategy = "strat_ubt",
     result = "icer",
     type = "difference"
)
```


```{r Model 30% HP New Plot Stool Removed }
plot(res_dsa_ubt_30_new,
     strategy = "strat_ubt",
     result = "icer",
     type = "difference"
)
```

```{r Model 60% HP New Plot Stool Removed }
plot(res_dsa_ubt_60_new,
     strategy = "strat_ubt",
     result = "icer",
     type = "difference"
)
```
```{r Model 90% HP New Plot Stool Removed }
plot(res_dsa_ubt_90_new,
     strategy = "strat_ubt",
     result = "icer",
     type = "difference"
)
```

```{r PSA Section}
# PSA Section
get_a <- function(mu, sd) {
  mu * (((mu * (1 - mu)) / (sd^2)) - 1)
}

get_b <- function(mu, sd) {
  (1 - mu) * (((mu * (1 - mu)) / (sd^2)) - 1)
}

get_sd <- function(ci, n) {
  (ci / 3.92) * sqrt(n)
}
```

```{r Define PSA}
rsp <- define_psa(
  ###########################
  # Epidemiology parameters #
  ###########################
  p_ulcer_male ~ beta(
    get_a(0.523, 0.1046),
    get_b(0.523, 0.1046)
  ),
  
  p_bleed_erad_first_year ~ beta(
    get_a(0.029, 0.022 / 2),
    get_b(0.029, 0.022 / 2)
  ),
  p_bleed_erad_sub_year ~ beta(
    get_a(0.0015, 0.0031 / 2),
    get_b(0.0015, 0.0031 / 2)
  ),
  
  p_bleed_non_erad_first_year ~ beta(
    get_a(0.2, 0.1 / 2),
    get_b(0.2, 0.1 / 2)
  ),
  p_bleed_non_erad_sub_year ~ beta(
    get_a(0.058, 0.017 / 2),
    get_b(0.058, 0.017 / 2)
  ),
  
  p_bleed_still_neg_first_year ~ beta(
    get_a(0.029, 0.022 / 2),
    get_b(0.029, 0.022 / 2)
  ),
  p_bleed_still_neg_sub_year ~ beta(
    get_a(0.0015, 0.0031 / 2),
    get_b(0.0015, 0.0031 / 2)
  ),
  p_bleed_pud ~ beta(
    get_a(0.015, 0.005),
    get_b(0.015, 0.005)
  ),
  p_pud_pos_lifetime ~ beta(
    get_a(0.15, 0.03),
    get_b(0.15, 0.03)
  ),
  
  ########################################
  # Test sensitivities and specificities #
  ########################################
  sens_stool ~ beta(
    get_a(0.87, 0.09 / 2),
    get_b(0.87, 0.09 / 2)
  ),
  spec_stool ~ beta(
    get_a(0.70, 0.16 / 2),
    get_b(0.70, 0.16 / 2)
  ),
  sens_ubt ~ beta(
    get_a(0.93, 0.05 / 2),
    get_b(0.93, 0.05 / 2)
  ),
  spec_ubt ~ beta(
    get_a(0.92, 0.09 / 2),
    get_b(0.92, 0.09 / 2)
  ),
  spec_sero ~ beta(
    get_a(0.88, 0.05 / 2),
    get_b(0.88, 0.05 / 2)
  ),
  sens_sero ~ beta(
    get_a(0.69, 0.14 / 2),
    get_b(0.69, 0.14 / 2)
  ),
  spec_rut ~ beta(
    get_a(0.93, 0.06 / 2),
    get_b(0.93, 0.06 / 2)
  ),
  sens_rut ~ beta(
    get_a(0.67, 0.06 / 2),
    get_b(0.67, 0.06 / 2)
  ),
  sens_hist ~ beta(
    get_a(0.70, 0.08 / 2),
    get_b(0.70, 0.08 / 2)
  ),
  spec_hist ~ beta(
    get_a(0.90, 0.09 / 2),
    get_b(0.90, 0.09 / 2)
  ),
  
  ####################################
  # Eradication and recurrence rates #
  ####################################
  p_erad_first ~ beta(
    get_a(0.70, 0.20),
    get_b(0.70, 0.20)
  ),
  p_erad_second ~ beta(
    get_a(0.81, 0.19),
    get_b(0.81, 0.19)
  ),
  p_reinfect ~ beta(
    get_a(0.0145, 0.0029),
    get_b(0.0145, 0.0029)
  ),
  p_infect ~ beta(
    get_a(0.0025, 0.00265),
    get_b(0.0025, 0.00265)
  ),
  
  #######################
  # Death probabilities #
  #######################
  p_death_bleed ~ beta(
    get_a(0.086, 0.056 / 2),
    get_b(0.086, 0.056 / 2)
  ),
  
  #########
  # Costs #
  #########
  c_biopsy ~ gamma(mean = 89, sd = 17.8),
  c_stool ~ gamma(mean = 14.38, sd = 2.88),
  c_sero ~ gamma(mean = 16.85, sd = 3.37),
  c_ubt_13 ~ gamma(mean = 105.47, sd = 21.09),
  c_first_line ~ gamma(mean = 119.18, sd = 23.84),
  c_second_line ~ gamma(mean = 136.62, sd = 27.32),
  c_hospital ~ gamma(mean = 17193.03, sd = 421.01),
  
  #############
  # Utilities #
  #############
  u_bleed ~ beta(
    get_a(0.9583, 0.0083),
    get_b(0.9583, 0.0083)
  ),
  u_pos ~ beta(
    get_a(0.9, 0.1),
    get_b(0.9, 0.1)
  )
)
```

```{r Set Cluster, Seed}
use_cluster(10)

set.seed(1)
pm <- run_psa(
  model = res_mod_17,
  psa = rsp,
  N = n_patients
)
print(pm)
```

```{r Plot PSA}
pm_sum_17 <- summary(pm, threshold = 0)
pm$run_model

plot(pm, type = "ce")

plot_ac_17 <- plot(pm, type = "ac", max_wtp = 100000, log_scale = FALSE, labels = c("strat_histology"="J0", "J2", "J3", "J4","J5", "J6"))
plot_ac_17

plot_ac_17$data$.p
as.data.table(plot_ac_17$data)[
  .strategy_names == "strat_histology"
]
```



# ```{r Plot PSA Evpi, Cov}
# plot(pm, type = "evpi", max_wtp = 100000, log_scale = FALSE)
# plot(pm, type = "cov")
# ```
# 
# ```{r Counts by State}
# plot(res_mod_17, type = "counts", panel = "by_state") +
#   xlab("Time") +
#   theme_bw() +
#   scale_color_brewer(
#     name = "Strategy",
#     palette = "Set1"
#   )
# ```

# ```{r BCEA}
# bcea <- run_bcea(pm, Kmax = 50000, wtp = 100000)
# ```
# 
# ```{r Summary BCEA}
# summary(bcea, wtp = 100000)
# ```
# 
# ```{r Plot BCEA}
# BCEA::ceac.plot(bcea)
# BCEA::ceaf.plot(bcea)
# BCEA::ceef.plot(
#     bcea,
#     relative = FALSE,
#     graph = "ggplot2"
# )
# ```
# 
# ```{r Plot BCEA Color Coded}
# BCEA::ceplane.plot(
#     bcea,
#     graph = "ggplot2",
#     point_colors = c(
#         "red", "green", "blue",
#         "yellow", "orange", "brown"
#     ),
#     wtp = 100000
# )
# ```
# 
# ```{r MAT Pretty 2 Define Transitions}
# mat_pretty_two <- define_transition(
#     state_names = c(
#         "state_bleed_pos", "state_bleed_neg", "state_bleed_erad",
#         "state_bleed_non_erad", "state_bleed_still_neg", "state_erad",
#         "state_non_erad", "state_neg", "state_death"
#     ),
#     # state_bleed_pos
#     0, 0, 0, (1 - p_death_bleed) * p_bleed_non_erad, 0, 0, C, 0, p_death_bleed,
#     # state_bleed_neg
#     0, 0, 0, 0, (1 - p_death_bleed) * p_bleed_still_neg, 0, 0, C, p_death_bleed,
#     # state_bleed_erad
#     0, 0, (1 - p_death_bleed) * p_bleed_erad, 0, 0, C, 0, 0, p_death_bleed,
#     # state_bleed_non_erad
#     0, 0, 0, (1 - p_death_bleed) * p_bleed_non_erad, 0, 0, C, 0, p_death_bleed,
#     # state_bleed_still_neg
#     0, 0, 0, 0, (1 - p_death_bleed) * p_bleed_still_neg, 0, 0, C, p_death_bleed,
#     # state_erad
#     0, 0, (1 - rr_death_bleed * mr) * p_bleed_erad, 0, 0, C, 0, 0,
#     rr_death_bleed * mr,
#     # state_non_erad
#     0, 0, 0, (1 - rr_death_bleed * mr) * p_bleed_non_erad, 0, 0, C, 0,
#     rr_death_bleed * mr,
#     # state_neg
#     0, 0, 0, 0, (1 - rr_death_bleed * mr) * p_bleed_still_neg, 0, 0, C,
#     rr_death_bleed * mr,
#     # state_death
#     0, 0, 0, 0, 0, 0, 0, 0, 1
# )
# ```
# 
# ```{r Plot MAT Pretty 2}
# plot(mat_pretty_two, arr.tcol = "transparent",
#     self.shiftx = c(0.1, 0.1, -0.1, -0.1, -0.1),
#     self.shifty = c(0.1, -0.1, -0.1, -0.1, -0.1),
#     self.arrpos = c(1.6, -1.6, -1.6, -1.6, -1.6)
# )
# ```
# 
# ```{r Rate Conversions}
# rate_to_prob(r = 0.15, per = 1, to = 1)
# rate_to_prob(r = 0.155, per = 1, to = 1)
# ```
# 
# ```{r Get Counts Model 17}
# get_counts(res_mod_17) |>
#     filter(.strategy_names == "strat_nh") |>
#     filter(model_time != 1) |>
#     filter(state_names == "state_bleed_non_erad") |>
#     pull(count) |>
#     sum()
# ```
