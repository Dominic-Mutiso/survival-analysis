# Progression-Free Survival (PFS) Efficacy Table in R
PFS can de defined as the time from randomisation/ initiation of treatament to disease progression or death from any cause. In most solid tumor oncology trials, Progression-Free Survival (PFS) is used as the primary endpoint. Since PFS involves time-to-event anylysis, we use ADTTE (Time-to-Event) ADaM dataset. Basic principles of Survival Analysis such as censoring, hazard and survival functions will play a big role in generating our PFS Efficacy table.

Additionally, we make use of advanced techniques for **Non-parametric estimation** (Kaplan-Meier), **Hypothesis Testing**: Non-parametric approach (Log-Rank Test) and **Regression Analysis**: Cox Proportional Hazards Model (CPHM). Therefore, familiarity with Survival Analysis is required!

## Instructions
Use the mock shell below to generate a PFS Efficacy table.
- Get the Kaplan-Meier median progression-free survival (PFS) and its 95% C.I based on Brookmeyer-Crowley method. 
- HR together with its 95% CI from a stratified Cox model with HRR status as a strata. The CI will be calculated using a profile likelihood approach. Use EFRON method to control for ties.
- For generation of the p-value, use a stratified log-rank test with HRR status (mutant versus wildtype) as a strata.


 ![image](https://github.com/user-attachments/assets/045e463b-a67d-44a2-86dc-a6f6741ee594)

# Data Manipulation
Worth noting:
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

## 1.1) CNSR Coding differences

### 1.1.1) Case 1: `0 = "Censored"` and `1 = "Event"`

 ```r
library (survival)

surv<-survfit(Surv(aval_months, status) ~ trt01pn_, data = adtte1,
                stype = 1, ctype = 1, conf.type = "log-log", conf.int = 0.95 )

print(surv, digits = 4)
names(surv)

```
     
### 1.1.2) **Case 2:** `0 = "Event"` and `1 = "Censored"`
CDISC ADaM ADTTE dataset has variable *cnsr* coded as 0 for "Event" and 1 for "Censored". Recall  R interprets 1 as "Event" and 0 as "Censored". If you do not want to recode, thanks to `Surv_CNSR()` function from `ggsurvfit` package that takes care of that. It uses CDISC ADaM ADTTE coding conventions- censor = 1, status/event = 0.

```r
install.packages("ggsurvfit")
library(ggsurvfit)

survfit(Surv_CNSR(aval_months, cnsr)~ trt01pn_, data = adtte1, conf.int = 0.95,
        stype = 1, ctype = 1, conf.type = "log-log")
```
## 1.2)  Median PFS Dataframe results
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
## 1.3) 25% 50% and 75% survival time and CI
```
quantile(surv, probs = c(25,50,75)/100)
```

# 2. Survival Probability at Specific Timepoint (6 Months)
For specific time points, we use the `summary()` function with the argument `times =`. If the specific time points are of no interest, use `surv_summary()` from the `survminer` package. Ultimately, to create a dataframe for the specific time points, one may use `names()` or `str()` to view the structure.

```r
#with all time points at which a curve has a step
surv_summary(surv)

#more than one specified timepoint- 0, 6, 12, 18, 24, 30
summary(surv, times = c(0, 6*(1:5)))

#one specific timepoint- 6.
tp<-summary(surv, times = 6)

tp

# Inspect the structure of the summary object
names(tp)

str(tp)

# Convert the summary to a data frame
time_pt <- data.frame(
                        time     = tp$time    
                      , n.risk   = tp$n.risk  
                      , n.event  = tp$n.event 
                      , survival = (tp$surv*100)    
                      , std.err  = tp$std.err 
                      , lower    = tp$lower * 100   
                      , upper    = tp$upper * 100  
                      , strata   = gsub("=", "", tp$strata)
)

# View the data frame
View(time_pt)

time_pt<- time_pt %>% mutate(
                            survival = sprintf("%10.1f", survival)
                           ,tp_ci   = paste( sprintf("%5.1f", lower)
                                           , sprintf("%5.1f", upper)
                                           , sep= " -"
                                            )
                            ) %>% select(strata, survival, tp_ci)

```
# 3. Stratified Log-rank test
Log-rank test- tests whether there is a significant difference in survival between two or more independent groups.

We can use either the `survdiff()` function or a combination of `survfit()` and `surv_pvalue()` to generate the **p-value**. A combination of `survfit()` and `surv_pvalue()` works best if you require the results for further processing.

```r
#***********************************************
#Option 1: Using "survdiff" function
#***********************************************
#a.)
survdiff(Surv(aval_months, status)~ trt01pn + strata(stratf1a), adtte1)

#b.)
pval<-survdiff(Surv_CNSR(aval_months, cnsr)~ trt01pn_ + strata(stratf1a), adtte1)
print(pval, digits = 6)

###Further manually;
pval_<- 1 - pchisq(pval$chisq, length(pval$n)-1)
pval_
```
```r
#**************************************************
#Option 2: Using both survfit() and surv_pvalue()
#**************************************************
surv1<-survfit(Surv(aval_months, status) ~ trt01pn_ + strata(stratf1a), data = adtte1, 
              stype = 1, ctype = 1, conf.int = 0.95, conf.type="log-log")

pvalue<-surv_pvalue(surv1)

View(pvalue)

pvalue<-pvalue %>% mutate (
                             pval   = sprintf("%12.3f", pval)
                           , strata = sub("\\+(.+)", "Active", variable)
                           ) %>% select(pval, strata)

View(pvalue)
```
The *p-value* from the `survdiff()` or `surv_pvalue()` function does not indicate which treatment is superior. It only tells us whether there is a significant difference between the survival curves. Therefore, we should use a survival plot to visually assess which treatment is superior.
```r
# Plot the survival curves
plot<-ggsurvplot(surv1, data = adtte1,
           legend.labs = c("placebo HRRm", "placebo HRRwt", "active HRRm", "active HRRwt"),
           xlab = "Time in days", 
           ylab = "Survival probability",
           ggtheme = theme_minimal(base_size = 15),
           size = 0.5,  # Reduce line size
           font.x = c(7),  # Reduce x-axis font size
           font.y = c(7),  # Reduce y-axis font size
           font.legend = c(7),  # Reduce legend font size
           font.tickslab = c(6),
           palette = c("red", "red", "blue", "blue"),
           linetype = c("dotted", "solid", "dotted", "solid"),
           legend = "right"
           )


plot$plot <- plot$plot +
        guides(color = guide_legend(override.aes = list(shape = NA)),  # Remove the censoring symbol
               linetype = guide_legend(override.aes = list(shape = NA)))

plot
```
![image](https://github.com/user-attachments/assets/97bb6116-620f-44d6-b84e-607595d92e09)


From this plot we can observe the following:
 1. Within each strata (HRR status), active drug is superior to placebo.
 2. Subjects with the wild type HRR status have better survival compared to those with mutant HRR status.
    
 # 4. Cox Proportional Hazards Model (CPHM)
Used to fit regression models to sensored survival data. It accounts for the effects of continuous and discrete covariate (independent
variable) measurements when the dependent variable is possibly censored time-to-event data.

 It assumes that the effects of the independent variables upon survival are constant over time and are additive in one scale.
 B&#770;

When an independent variable is categorical, it’s important to choose a baseline/reference group, as hazard ratios are calculated by comparing the other levels to the reference group.

*Hazard function* is the conditional probability that an event will occur at time *t* having survived to that time.

## 4.1) Unadjusted Cox Regression
![image](https://github.com/user-attachments/assets/dcb06943-9a27-42b0-9408-787674693027)

The estimate of the log hazard ratio treatment effect, B&#770; is 0.05214. Since this is positive, higher hazards are associated with the active drug than with the placebo. That is, the active drug appears to reduce survival- quite unfortunate. The value of exp(B&#770;) is 1.054 meaning the risk of death/disease progression is higher on the active drug by about 5.4%. 

## 4.2) Adjusted Cox Regression
### Stratified Cox Model
![image](https://github.com/user-attachments/assets/7cf30773-cee8-4dc4-bd8a-d5a01f9472d4)

The coefficient is now negative. This means within each stratification, the active drug is effective. We can also  explicitly estimate the stratification effect.

![image](https://github.com/user-attachments/assets/cf999813-23ac-4b23-9916-eee2dff75435)

Again the co-efficient associated with the active drug is negative. The stratum HRRwt has negative co-efficient and thus is associated with lower hazard than the reference HRRm.

## Options of viewing results.
  1. `summary ()`- used to view results in the console.
  2. `tbl_regression ()`- comes from `gtsummary` package. Used to view results in a table.
  3. `tidy ()`- comes from the `broom` library. View results as a dataframe. Useful if you need the results for further processing. Use `??tidy.coxph` to explore this fuction further.
     Syntax: *tidy(x, exponentiate = FALSE, conf.int = FALSE, conf.level = 0.95, ...)*



