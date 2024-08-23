# Progression-Free Survival (PFS) Efficacy Table using R
PFS can de defined as the time from randomisation/ initiation of treatament to disease progression or death from any cause. In most solid tumor oncology trials, Progression-Free Survival (PFS) is used as the primary endpoint. Since PFS involves time-to-event anylysis, we use ADTTE (Time-to-Event) ADaM dataset. Basic principles of Survival Analysis such as censoring, hazard and survival functions will play a big role in generating our PFS Efficacy table.

Additionally, we make use of advanced techniques for **Non-parametric estimation** (Kaplan-Meier), **Hypothesis Testing**: Non-parametric approach (Log-Rank Test) and **Regression Analysis**: Cox Proportional Hazards Model (CPHM). Therefore, familiarity with Survival Analysis is required!

## Instructions
We will use the mock shell below to generate our PFS Efficacy table.
- Get the Kaplan-Meier median progression-free survival (PFS) and its 95% C.I based on Brookmeyer-Crowley method. 
- HR together with its 95% CI from a stratified Cox model with HRR status as a strata. The CI will be calculated using a profile likelihood approach. Use EFRON method to control for ties.
- For generation of the p-value, use a stratified log-rank test with HRR status (mutant versus wildtype) as a strata.


 ![image](https://github.com/user-attachments/assets/045e463b-a67d-44a2-86dc-a6f6741ee594)

