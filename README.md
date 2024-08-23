# Progression-Free Survival (PFS) Efficacy Table in R
PFS can de defined as the time from randomisation/ initiation of treatament to disease progression or death from any cause. In most solid tumor oncology trials, Progression-Free Survival (PFS) is used as the primary endpoint. Since PFS involves time-to-event anylysis, we use ADTTE (Time-to-Event) ADaM dataset. Basic principles of Survival Analysis such as censoring, hazard and survival functions will play a big role in generating our PFS Efficacy table.

Additionally, we make use of advanced techniques for **Non-parametric estimation** (Kaplan-Meier), **Hypothesis Testing**: Non-parametric approach (Log-Rank Test) and **Regression Analysis**: Cox Proportional Hazards Model (CPHM). Therefore, familiarity with Survival Analysis is required!

## Instructions
We will use the mock shell below to generate our PFS Efficacy table.
- Get the Kaplan-Meier median progression-free survival (PFS) and its 95% C.I based on Brookmeyer-Crowley method. 
- HR together with its 95% CI from a stratified Cox model with HRR status as a strata. The CI will be calculated using a profile likelihood approach. Use EFRON method to control for ties.
- For generation of the p-value, use a stratified log-rank test with HRR status (mutant versus wildtype) as a strata.


 ![image](https://github.com/user-attachments/assets/045e463b-a67d-44a2-86dc-a6f6741ee594)

# Data Manipulation
NB:
- `Surv()` function in the `{survival}` package accepts by default T/F, where TRUE is the event and False stands for censored. That is, 1 for Event and 0 for censored. Event in our case (PFS) is Death/ Disease Progression coded as 0. Recode to 1 so that R does not take it for "censored". Take care to ensure the Event is properly coded as 1.
- We also need to convert some of the variables to type "factor". Always the first level is used as the baseline/ reference in R. Use `factor()` function to order levels using 'levels = ' argument. Else use `relevel(, ref = "")` E.g: `trt01pn = relevel(as.factor(trt01pn), ref= "")`. Trt01pn = 1 rep' 'Active drug' and Trt01pn = 2 rep' the 'Placebo'. Failure to use 'levels = " or relevel(), R will by default choose alphabetically or numerically. Consequently, the first occurence which is '1' for Active drug will be the baseline/ ref'. We don't want that! We are comparing Active drug to Placebo.

```r
adtte1<-adtte %>% rename_with(tolower) %>% 
                  select( usubjid, trt01a, trt01an, aval, param, paramcd, cnsr
                         , evntdesc, stratf1a, stratf2a, trt01pn, trt01p, parqual
                         ) %>% 
                  mutate(
                          aval_months = aval/30.4375
                        , status      = ifelse(cnsr == 0, 1, 0)
                        , trt01pn_    = factor(trt01pn, levels = c("2", "1"),
                                               labels = c("Placebo","Active"))
                        , stratf1a    = relevel(as.factor(stratf1a), ref = "HRRm")
                        , check       = regexpr("Progression-free survival",param)
                        ) %>% 
                filter (grepl("Progression-free survival", param) == "TRUE" & paramcd == "TRPROGT")
```
# 1. Median Progression Free Survival
Median PFS is the value of time t where the survival function, S(t), equals 0.5. That is, 50 % of the cohort is event free.

Although Kaplan-Meier survival curves are calculated  independently for each group, meaning the order or arrangement of factor levels do not affect the computation of the survival estimates or their confidence intervals, it is good practice to explicitly order the factor levels using `factor(, levels = c())` or `relevel(, ref = "")`. Additionally, any numeric variable intended for grouping should be converted to a factor, rather than allowing R to automatically interpret it for grouping. That is what we did in the previous step "**data manipulation**" for the numeric variable `trt01pn`.

Argument `type` in `survfit` function is an older argument that combined `stype` and `ctype`, now deprecated. Legal values were "kaplan-meier" which is equivalent to `stype = 1`, `ctype = 1`, `"fleming-harrington"` which is equivalent to `stype = 2`, `ctype = 1`, and `"fh2"` which is equivalent to `stype = 2`and `ctype = 2`.

Load the `survival` package for this step.

### 1.1) Case 1: `0 = "Censored"` and `1 = "Event"`

 ```r
library (survival)

surv<-survfit(Surv(aval_months, status) ~ trt01pn_, data = adtte1,
                stype = 1, ctype = 1, conf.type = "log-log", conf.int = 0.95 )

print(surv, digits = 4)
names(surv)

```
     
### 1.2) **Case 2:** `0 = "Event"` and `1 = "Censored"`
CDISC ADaM ADTTE dataset has variable *cnsr* coded as 0 for "Event" and 1 for "Censored". Recall  R interprets 1 = "Event" and 0 = "Censored". If you do not want to recode, thanks to `Surv_CNSR()` function from `ggsurvfit` package that takes care of that. It uses CDISC ADaM ADTTE coding conventions- censor = 1, status/event = 0.

```r
install.packages("ggsurvfit")
library(ggsurvfit)

survfit(Surv_CNSR(aval_months, cnsr)~ trt01pn_, data = adtte1, conf.int = 0.95,
        stype = 1, ctype = 1, conf.type = "log-log")
```
### 1.3)  Median PFS Dataframe results
If you require your results for further processing use `surv_median()` function from `survminer` package. One can use results from either **case 1** or **case 2**.

```r
install.packages("survminer")

library(survminer)

med_surv<-surv_median(surv)

View(med_surv)

med_surv<- med_surv %>% mutate (
                                strata = gsub("=", "", strata)
                              , median = sprintf("%10.1f", median)
                              , med_ci = paste(sprintf("%5.1f", lower), 
                                               trimws(sprintf("%5.1f", upper)), 
                                               sep = " - " 
                                              )
                               ) %>% select (-c(lower, upper))
                       
View(med_surv)

```
### 1.4) 25% 50% and 75% survival time and CI
```
quantile(surv, probs = c(25,50,75)/100)
```






