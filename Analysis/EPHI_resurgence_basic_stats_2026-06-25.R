#################################################################################
#################### Exploratory statistical comparisons ########################
#################################################################################

# Jonathan Parr
# 6/26/2026

# Setting up working directory
setwd("Ethiopia_2024_Resurgence_Investigation_repo") 
# Update to your local path

#### Background
#Exploratory comparisons of resistance marker prevalence estimates at our sites to published estimates from the following studies:
#  - **Brhane et al. Comms Medicine 2025**
#  -- *pfk13* 622I: overall (5 regions) 2019, 15.7% [95% CI 13.2–18.2%]

#  - **Fola et al. Nat Micro 2023**
#  -- *pfk13* 622I: Tigray 2018, 8.4% [95% CI 6.2–10.5%]
#  -- *pfk13* 622I: Amhara 2018, 46.9% [95% CI 40.7–53.2%]

#Using Wald z-test to compare our results to published prevalence estimates, with standard errors derived from published 95% confidence intervals.


#### Importing data

# Current study - overall
outbreak_all_622I_p <- 0.355      # weighted prevalence
outbreak_all_622I_l <- 0.282      # lower 95% CI
outbreak_all_622I_u <- 0.428       # upper 95% CI

# Current study - Tigray only
outbreak_tig_622I_p <- 0.671
outbreak_tig_622I_l <- 0.514      
outbreak_tig_622I_u <- 0.829       

# Current study - Amhara only
outbreak_amh_622I_p <- 0.332
outbreak_amh_622I_l <- 0.236      
outbreak_amh_622I_u <- 0.427    

# Brhane (overall comparison)
brhane_all_622I_p <- 0.157     
brhane_all_622I_l <- 0.132      
brhane_all_622I_u <- 0.182      

# Fola (Tigray comparison)
fola_tig_622I_p <- 0.084     
fola_tig_622I_l <- 0.062      
fola_tig_622I_u <- 0.105   

# Fola (Amhara comparison)
fola_amh_622I_p <- 0.469    
fola_amh_622I_l <- 0.407      
fola_amh_622I_u <- 0.532   

# # Zeleke (Amhara comparison)
# zeleke_amh_622I_p <- 0.443     
# zeleke_amh_622I_l <- 0.409     
# zeleke_amh_622I_u <- 0.477  


#### Function to perform Wald z-test
compare_prevalence <- function(p_ref, l_ref, u_ref, p_comp, l_comp, u_comp) {
  # Convert CIs to SEs
  se_ref <- (u_ref - l_ref) / (2 * 1.96)
  se_comp <- (u_comp - l_comp) / (2 * 1.96)
  
  # Wald test
  diff <- p_ref - p_comp
  se_diff <- sqrt(se_ref^2 + se_comp^2)
  z <- diff / se_diff
  p_value <- 2 * (1 - pnorm(abs(z)))
  
  # Return results
  return(data.frame(
    Difference = round(diff, 4),
    SE = round(se_diff, 4),
    Z = round(z, 3),
    P_value = round(p_value, 8)
  ))
}


#### Wald z-test comparing *pfk13* 622I prevalence in our study to the published estimates above

# Perform three comparisons
cat("Comparison to Brhane 2019, overall:\n")
print(compare_prevalence(outbreak_all_622I_p, outbreak_all_622I_l, outbreak_all_622I_u, brhane_all_622I_p, brhane_all_622I_l, brhane_all_622I_u))

cat("\nComparison to Fola 2018, Tigray sites:\n")
print(compare_prevalence(outbreak_tig_622I_p, outbreak_tig_622I_l, outbreak_tig_622I_u, fola_tig_622I_p, fola_tig_622I_l, fola_tig_622I_u))

cat("\nComparison to Fola 2018, Amhara sites:\n")
print(compare_prevalence(outbreak_amh_622I_p, outbreak_amh_622I_l, outbreak_amh_622I_u, fola_amh_622I_p, fola_amh_622I_l, fola_amh_622I_u))