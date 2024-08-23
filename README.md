# Progression-Free Survival (PFS) Efficacy Table using R
PFS can de defined as the time from randomisation/ initiation of treatament to disease progression or death from any cause. In most solid tumor oncology trials, Progression-Free Survival (PFS) is used as the primary endpoint. Since PFS involves time-to-event anylysis, we use ADTTE (Time-to-Event) ADaM dataset. Basic principles of Survival Analysis such as censoring, hazard and survival functions will play a big role in generating our PFS Efficacy table.

Additionally, we make use of advanced techniques for **Non-parametric estimation** (Kaplan-Meier), **Hypothesis Testing**: Non-parametric approach (Log-Rank Test) and **Regression Analysis**: Cox Proportional Hazards Model (CPHM). Therefore, familiarity with Survival Analysis is required!

## Overview
 

We need to get the Kaplan-Meier median progression-free survival (PFS) and 95% C.I
