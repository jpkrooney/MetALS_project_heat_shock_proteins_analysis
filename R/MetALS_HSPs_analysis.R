library(tidyverse)
library(ggbeeswarm)
library(ggcorrplot)
library(ggpubr)
library(tableone)
library(corrplot)
library(lme4)
library(broom)
library(broom.mixed)
library(patchwork)
library(survival)


# Load datasets
dfall <- readRDS("../../Data/Processed_data/final_collated_data_allcohorts.RDS")


# Define variables lists
demo_clin <- c("case", "sex", "height_cm", "cohort", "weight_kg", "bmi", 
               "site_onset", "age_onset", "age_dx", "age_survey", "yr_onset", 
               "dx_delay", "isced1997", "el_escorial", "alsfrs_total")


#### Subset to only those with at least one HSP ####
dfall <- dfall %>% 
    filter(!(is.na(hsp70_ngperml) & is.na(hsp90_ngperml) & is.na(dnajc7_ngperml)))
dfall$outcome <- ifelse(dfall$case == "Patient", 1, 0)



# Extract first measurement for each person (note only Germany and Ireland2 have repeated measures)
dfbase <- dfall %>% 
    arrange(id, time_svy2lab) %>% 
    group_by(id) %>% 
    filter(row_number()==1) %>% 
    ungroup()
dfbase$case <- as.factor(dfbase$case)

# Define survival outcome variable
dfbase <- dfbase %>%
    mutate(survoutcome = case_when(vital_status %in% c("Died", "Tracheotomy") ~ 1,
                                   vital_status == "Alive" ~ 0))


#### Descriptive statistics ####
# Table 1
#### Make Table 1 of descriptive statistics by cohort ####
non_par <- c("onset2censor_months", "count_long", "baseline2censor_months",
             "baseline2long_end", "alsfrs_total", "dx_delay")
cat_vars <- c("case", "sex", "site_onset", 
              "isced1997", "el_escorial")
tab1 <- CreateTableOne(demo_clin[ demo_clin != "cohort"],
                       strata = c("case", "cohort"), factorVars = cat_vars, data = dfbase )
tab1b <- print(tab1, nonnormal = non_par, exact = c("sex", "site_onset", "staging", "isced1997"),
               quote = FALSE, noSpaces = TRUE, includeNA = TRUE,
               printToggle = FALSE, showAllLevels = TRUE)

#### Make Table 1 of descriptive statistics overall ####
tab1_total <- CreateTableOne(demo_clin[ demo_clin != "cohort"], 
                             strata = c("case"), factorVars = cat_vars, data = dfbase)
tab1_total <- print(tab1_total, nonnormal = non_par, exact = c("sex", "site_onset", "staging", "isced1997"),
                    quote = FALSE, noSpaces = TRUE,
                    printToggle = FALSE, showAllLevels = TRUE)

# Merge by cohort and overall summary stats and output to csv file
table1 <- cbind(tab1b, tab1_total)
write.csv(table1, "Results/BaselineDemo_table1.csv")



#### Calculate correlations of biomarkers ####
corrvars <- c("height_cm", "weight_kg",  "age_survey", "hsp70_ngperml",
              "hsp90_ngperml",  "dnajc7_ngperml")
corsControls <- cor(dfall[dfall$case == "Control", corrvars ],
    use = "pairwise.complete.obs",)
corsPatients <- cor(dfall[dfall$case == "Patient", corrvars ],
                    use = "pairwise.complete.obs",)
# Make correlation plots
Fig2a <- ggcorrplot(corsControls, type = "lower", lab = TRUE) + ggtitle("Controls")
Fig2b <- ggcorrplot(corsPatients, type = "lower", lab = TRUE) + ggtitle("ALS Cases")
Fig2a + Fig2b + plot_layout(guides = "collect")


# Supplementary Figure 2
tiff("Graphs/SuppFigure2.tiff", width = 800, height = 460)
    Fig2a + Fig2b + plot_layout(guides = "collect")
dev.off()



#### Calculate descriptive statistics of biochemical vars ####
biovars <- c("hsp70_ngperml", "hsp90_ngperml", "dnajc7_ngperml",
             "creatinine_mgperdl", "cystatinc_mgper_l")
quants <- c(0.05, 0.25, 0.5, 0.75, 0.95)

dfquantiles <- dfall %>% 
        group_by(case,cohort) %>%
        reframe(enframe(quantile(hsp70_ngperml, quants,
                                   na.rm=TRUE), "quantile", "hsp70"),
                enframe(quantile(hsp90_ngperml, quants,
                                   na.rm=TRUE), "quantile", "hsp90"),
                enframe(quantile(dnajc7_ngperml, quants,
                                 na.rm=TRUE), "quantile", "dnajc7"),
                enframe(quantile(creatinine_mgperdl, quants,
                                 na.rm=TRUE), "quantile", "creatinine"),
                enframe(quantile(cystatinc_mgper_l, quants,
                                 na.rm=TRUE), "quantile", "cystatinc")) %>%
    pivot_wider(names_from = quantile,
                values_from = c(hsp70, hsp90, dnajc7, creatinine, cystatinc),
                names_prefix = "P") %>% 
    #arrange(cohort, case) %>% 
    data.frame()

dfquantiles <- dfquantiles %>%
    left_join(dfall %>% 
                  group_by(case,cohort) %>%
                  reframe(count_hsp70 = sum(!is.na(hsp70_ngperml)),
                          count_hsp90 = sum(!is.na(hsp90_ngperml)),
                          count_dnajc7 = sum(!is.na(dnajc7_ngperml)),
                          count_creatinine = sum(!is.na(creatinine_mgperdl)),
                          count_cystatinc = sum(!is.na(cystatinc_mgper_l))))

dfquantiles <- dfquantiles %>% 
    dplyr::select(case, cohort,
                  count_hsp70, hsp70_P5., hsp70_P25., hsp70_P50., hsp70_P75., hsp70_P95.,
                  count_hsp90, hsp90_P5., hsp90_P25., hsp90_P50., hsp90_P75., hsp90_P95.,
                  count_dnajc7, dnajc7_P5., dnajc7_P25., dnajc7_P50., dnajc7_P75., dnajc7_P95.,
                  count_creatinine, creatinine_P5., creatinine_P25., creatinine_P50.,
                  creatinine_P75., creatinine_P95.,
                  count_cystatinc, cystatinc_P5., cystatinc_P25., cystatinc_P50.,
                  cystatinc_P75., cystatinc_P95.,) %>% 
    arrange(cohort, case)

write.csv(dfquantiles, "Results/SuppTable1_biomarkerquantiles.csv")



#### Create a dataframe with Germany for baseline case control comparions as there are no German controls ####
dfbase_noDE <- dfbase %>% 
    filter(cohort != "Germany")


#### Create a longitudinal subjects only dataframe ####
dflong <- dfall %>% 
    group_by(id) %>% 
    filter(n() > 1 & cohort %in% c("Ireland - MetALS", "Germany")) #%>% 
    #mutate(cohort = as.factor(as.character(cohort)))


#### Plot distributions of HSPs at baseline ####
gFig1A <- ggplot(dfbase_noDE,
       aes(x = case, y = hsp70_ngperml)) + geom_beeswarm(alpha = 0.5) +
    stat_boxplot(geom = "errorbar", width = 0.25) +
    geom_boxplot(alpha=0.1, width = 0.5) +
    scale_y_continuous(breaks = seq(0, 40, by = 5), trans = "log") +
    xlab("") + ylab("log( HSP70 in ng/ml )") + ggtitle("A") +
    stat_compare_means(method = "t.test",
                       aes(label = sprintf("p = %5.3f", as.numeric(..p.format..))),
                       label.x = 1.3,
                       size = 5) +
    theme_bw(base_size = 14)
gFig1B <- ggplot(dfbase_noDE,
             aes(x = case, y = hsp90_ngperml)) + geom_beeswarm(alpha = 0.5) +
    stat_boxplot(geom = "errorbar", width = 0.25) +
    geom_boxplot(alpha=0.1, width = 0.5) +
    scale_y_continuous(limits = c(10, 75), breaks = seq(0, 75, by = 15), trans = "log") +
    xlab("") + ylab("log( HSP90 in ng/ml )") + ggtitle("B") +
    stat_compare_means(method = "t.test",
                       aes(label = sprintf("p = %5.3f", as.numeric(..p.format..))),
                       label.x = 1.3,
                       size = 5) +
    theme_bw(base_size = 14)
gFig1C <- ggplot(dfbase_noDE,
             aes(x = case, y = dnajc7_ngperml)) + geom_beeswarm(alpha = 0.5) +
    stat_boxplot(geom = "errorbar", width = 0.25) +
    geom_boxplot(alpha=0.1, width = 0.5) +
    scale_y_continuous(limits = c(1, 20), breaks = c(1, seq(0, 20, by = 5)), trans = "log") +
    xlab("") + ylab("log( DNAJC7 in ng/ml )") + ggtitle("C") +
    stat_compare_means(method = "t.test",
                       aes(label = sprintf("p = %5.3f", as.numeric(..p.format..))),
                       label.x = 1.3,
                       size = 5) +
    theme_bw(base_size = 14)
gFig1A + gFig1B + gFig1C

# Figure 1
tiff("Graphs/Figure1.tiff", width = 800, height = 460)
    gFig1A + gFig1B + gFig1C
dev.off()



gSuppFig1A <- ggplot(dfbase,
                 aes(x = case, y = hsp70_ngperml)) + geom_beeswarm(alpha = 0.5) +
    stat_boxplot(geom = "errorbar", width = 0.25) +
    geom_boxplot(alpha=0.1, width = 0.5) +
    scale_y_continuous(breaks = seq(0, 40, by = 5), trans = "log") +
    xlab("") + ylab("log( HSP70 in ng/ml )") + ggtitle("A") +
    stat_compare_means(method = "t.test",
                       aes(label = sprintf("p = %5.3f", as.numeric(..p.format..))),
                       label.x = 1.3,
                       label.y = log(37.5),
                       size = 5) +
    theme_bw(base_size = 14) + facet_wrap(~cohort, nrow = 1)

gSuppFig1B <- ggplot(dfbase,
                 aes(x = case, y = hsp90_ngperml)) + geom_beeswarm(alpha = 0.5) +
    stat_boxplot(geom = "errorbar", width = 0.25) +
    geom_boxplot(alpha=0.1, width = 0.5) +
    scale_y_continuous(limits = c(10, 75), breaks = seq(0, 75, by = 15), trans = "log") +
    xlab("") + ylab("log( HSP90 in ng/ml )") + ggtitle("B") +
    stat_compare_means(method = "t.test",
                       aes(label = sprintf("p = %5.3f", as.numeric(..p.format..))),
                       label.x = 1.3,
                       label.y = log(72.5),
                       size = 5) +
    theme_bw(base_size = 14) + facet_wrap(~cohort, nrow = 1)

gSuppFig1C <- ggplot(dfbase,
                 aes(x = case, y = dnajc7_ngperml)) + geom_beeswarm(alpha = 0.5) +
    stat_boxplot(geom = "errorbar", width = 0.25) +
    geom_boxplot(alpha=0.1, width = 0.5) +
    scale_y_continuous(limits = c(1, 20), breaks = c(1, seq(0, 20, by = 5)), trans = "log") +
    xlab("") + ylab("log( DNAJC7 in ng/ml )") + ggtitle("C") +
    stat_compare_means(method = "t.test",
                       aes(label = sprintf("p = %5.3f", as.numeric(..p.format..))),
                       label.x = 1.3,
                       label.y = log(19),
                       size = 5) +
    theme_bw(base_size = 14) + facet_wrap(~cohort, nrow = 1)
plot_spacer() + (gSuppFig1A / gSuppFig1B / gSuppFig1C) + plot_spacer()


tiff("Graphs/SuppFigure1.tiff", width = 800, height = 800)
    gSuppFig1A / gSuppFig1B / gSuppFig1C
dev.off()



#### Make scatter plots of HSPs vs one another ####
# fit some linear models to get linear slopes
m1a <- lm(log(dnajc7_ngperml) ~ log(hsp70_ngperml), dfall[ dfall$case == "Control", ])
m1b <- lm(log(dnajc7_ngperml) ~ log(hsp70_ngperml), dfall[ dfall$case == "Patient", ])
m2a <- lm(log(hsp90_ngperml) ~ log(hsp70_ngperml), dfall[ dfall$case == "Control", ])
m2b <- lm(log(hsp90_ngperml) ~ log(hsp70_ngperml), dfall[ dfall$case == "Patient", ])
m3a <- lm(log(dnajc7_ngperml) ~ log(hsp90_ngperml), dfall[ dfall$case == "Control", ])
m3b <- lm(log(dnajc7_ngperml) ~ log(hsp90_ngperml), dfall[ dfall$case == "Patient", ])

# the plots
g_hsp70_dnajc7 <- ggplot(dfall,
                         aes(x = hsp70_ngperml, y = dnajc7_ngperml)) +
    geom_point() + facet_wrap(~case) + geom_point() +
    scale_x_log10() + scale_y_log10() +
    xlab("log( HSP70 in ng/ml )") + ylab("log( DNAJC7 in ng/ml )") + ggtitle("A") +
    theme_bw(base_size = 14) + geom_smooth(method = "lm") + 
    geom_text(data = data.frame(hsp70_ngperml = 15, dnajc7_ngperml = 16,
                                case = c("Control", "Patient"),
                                label = c(paste0("n = ", nobs(m1a), "\nSlope: ", round(tidy(m1a)[2,]$estimate, 2),
                                                 ", P <0.001"),
                                          paste0("n = ", nobs(m1b), "\nSlope: ", round(tidy(m1b)[2,]$estimate, 2),
                                                                  ", P <0.001"))),
              aes(label = label), fontface='bold')
    #annotate(geom="text", x=15, y=10,
    #         label="Started\nMg\n400mg od", fontface="bold")

g_hsp70_hsp90 <- ggplot(dfall,
                        aes(x = hsp70_ngperml, y = hsp90_ngperml)) +
    geom_point() + facet_wrap(~case) + geom_point() +
    scale_x_log10() + scale_y_log10() +
    xlab("log( HSP70 in ng/ml )") + ylab("log( HSP90 in ng/ml )") + ggtitle("B") +
    theme_bw(base_size = 14) + geom_smooth(method = "lm") +
    geom_text(data = data.frame(hsp70_ngperml = 15, hsp90_ngperml = 40,
                                case = c("Control", "Patient"),
                                label = c(paste0("n = ", nobs(m2a), "\nSlope: ", round(tidy(m2a)[2,]$estimate, 2),
                                                 ", P = ", round(tidy(m2a)[2,]$p.value, 3)),
                                          paste0("n = ", nobs(m2b), "\nSlope: ", "Slope: ", round(tidy(m2b)[2,]$estimate, 2),
                                                 ", P <0.001"))),
              aes(label = label), fontface='bold')

g_hsp90_dnajc7 <- ggplot(dfall,
                         aes(x = hsp90_ngperml, y = dnajc7_ngperml)) +
    geom_point() + facet_wrap(~case) + geom_point() +
    scale_x_log10() + scale_y_log10() +
    xlab("log( HSP90 in ng/ml )") + ylab("log( DNAJC7 in ng/ml )") + ggtitle("C") +
    theme_bw(base_size = 14) + geom_smooth(method = "lm") +
    geom_text(data = data.frame(hsp90_ngperml = c(10, 30), dnajc7_ngperml = 15,
                                case = c("Control", "Patient"),
                                label = c(paste0("n = ", nobs(m3a), "\nSlope: ", round(tidy(m3a)[2,]$estimate, 2),
                                                 ", P = ", round(tidy(m3a)[2,]$p.value, 3)),
                                          paste0("n = ", nobs(m3b), "\nSlope: ", round(tidy(m3b)[2,]$estimate, 2),
                                                 ", P = ", round(tidy(m3b)[2,]$p.value, 3)))),
              aes(label = label), fontface='bold')

plot_spacer() + (g_hsp70_dnajc7 / g_hsp70_hsp90 / g_hsp90_dnajc7) + plot_spacer()

tiff("Graphs/Figure2.tiff", width = 800, height = 800)
    g_hsp70_dnajc7 / g_hsp70_hsp90 / g_hsp90_dnajc7
dev.off()


#### Plot ratios of HSPs
gSuppFig3A <- ggplot(dfall,
       aes(x = case, y = log(dnajc7_ngperml)/log(hsp70_ngperml))) +
    geom_boxplot() + geom_beeswarm(alpha = 0.5) +
    xlab("") + ylab("Ratio log(DNAJC7) / log(HSP70)") + ggtitle("A") +
    stat_compare_means(method = "t.test",
                       aes(label = sprintf("p = %5.3f", as.numeric(..p.format..))),
                       label.x = 1.3,
                       size = 5) +
    theme_bw(base_size = 14)
gSuppFig3B <- ggplot(dfall,
       aes(x = case, y = log(hsp90_ngperml)/log(hsp70_ngperml))) +
    geom_boxplot() + geom_beeswarm(alpha = 0.5) +
    xlab("") + ylab("Ratio log(HSP90) / log(HSP70)") + ggtitle("B") +
    stat_compare_means(method = "t.test",
                       aes(label = sprintf("p = %5.3f", as.numeric(..p.format..))),
                       label.x = 1.3,
                       size = 5) +
    theme_bw(base_size = 14)
gSuppFig3C <- ggplot(dfall,
       aes(x = case, y = log(dnajc7_ngperml)/log(hsp90_ngperml))) +
    geom_boxplot() + geom_beeswarm(alpha = 0.5) +
    xlab("") + ylab("Ratio log(DNAJC7) / log(HSP90)") + ggtitle("C") +
    stat_compare_means(method = "t.test",
                       aes(label = sprintf("p = %5.3f", as.numeric(..p.format..))),
                       label.x = 1.3,
                       size = 5) +
    theme_bw(base_size = 14)

gSuppFig3A + gSuppFig3B + gSuppFig3C

tiff("Graphs/SuppFigure3.tiff", width = 800, height = 460)
    gSuppFig3A + gSuppFig3B + gSuppFig3C
dev.off()



#### Plot distributions of hsps vs onset
gFig3A <- ggplot(dfall %>%
           filter(time_onset2lab < 200),
       aes(x = time_onset2lab, y = hsp70_ngperml, group = id)) +
    geom_point() + geom_line() + scale_y_log10() +
    xlab("") + ylab("ng/ml \n(log scaled)") + ggtitle("HSP70") +
    #facet_wrap(~cohort, scales = "free", ncol=1) +
    theme_bw(base_size = 12)

gFig3B <- ggplot(dfall %>%
           filter(time_onset2lab < 200 & cohort != "Italy"),
       aes(x = time_onset2lab, y = hsp90_ngperml, group = id)) +
    geom_point() + geom_line() + scale_y_log10() +
    xlab("") + ylab("ng/ml \n(log scaled)") + ggtitle("HSP90") +
    #facet_wrap(~cohort, scales = "free", ncol=1) +
    theme_bw(base_size = 12)

gFig3C <- ggplot(dfall %>% filter(time_onset2lab < 200),
       aes(x = time_onset2lab, y = dnajc7_ngperml, group = id)) +
    geom_point() + geom_line() + scale_y_log10() +
    xlab("Months from ALS onset") + ylab("ng/ml \n(log scaled)") + ggtitle("DNAJC7") +
    #facet_wrap(~cohort, scales = "free", ncol=1) +
    theme_bw(base_size = 12)

gFig3A / gFig3B / gFig3C

tiff("Graphs/Figure3.tiff", width = 600, height = 800)
    gFig3A / gFig3B / gFig3C
dev.off()


#### Plot HSPs versus ALSFRS score
gSuppFig4a <- ggplot(dfall,
       aes(x = alsfrs_total, y = hsp70_ngperml, group = id)) + 
    geom_point() + geom_line() + scale_y_log10() +
    #expand_limits(x = 0) +
    xlab("Total ALSFRS-r") + ylab("ng/ml") + ggtitle("HSP70") +
    theme_bw(base_size = 12)

gSuppFig4b <- ggplot(dfall,
       aes(x = alsfrs_total, y = log(hsp90_ngperml), group = id)) +
    geom_point() + geom_line() + scale_y_log10() +
    #expand_limits(x = 0) +
    xlab("Total ALSFRS-r") + ylab("ng/ml") + ggtitle("HSP90") +
    theme_bw(base_size = 12)

gSuppFig4c <- ggplot(dfall,
       aes(x = alsfrs_total, y = log(dnajc7_ngperml), group = id)) +
    geom_point() + geom_line() + scale_y_log10() +
    #expand_limits(x = 0) +
    xlab("Total ALSFRS-r") + ylab("ng/ml") + ggtitle("DNAJC7") +
    theme_bw(base_size = 12)

gSuppFig4a / gSuppFig4b / gSuppFig4c

tiff("Graphs/SuppFigure4_alsfrsr.tiff", width = 600, height = 800)
    gSuppFig4a / gSuppFig4b / gSuppFig4c
dev.off()


### Make flag variable for rows without missing variables for main analysis ###
dfbase_noDE <- dfbase_noDE %>% 
    mutate(HSP70_NA = if_else(!is.na(log(hsp70_ngperml)) & !is.na(sex) & !is.na(age_survey), FALSE, TRUE),
           HSP90_NA = if_else(!is.na(log(hsp90_ngperml)) & !is.na(sex) & !is.na(age_survey), FALSE, TRUE),
           DNAJC7_NA = if_else(!is.na(log(dnajc7_ngperml)) & !is.na(sex) & !is.na(age_survey), FALSE, TRUE),
           HSP70DNAJC7_NA = as.logical(HSP70_NA * DNAJC7_NA)) %>% 
    data.frame()


#### Fit main analysis glm models ####
# cohort to be included as random effect.
glm_hsp70 <- glmer(outcome ~ log(hsp70_ngperml) + sex + age_survey + (1|cohort),
                   data = dfbase_noDE %>% 
                       filter(HSP70_NA == FALSE), family = "binomial")
glm_hsp90 <- glmer(outcome ~ log(hsp90_ngperml) + sex + age_survey + (1|cohort),
                   data = dfbase_noDE %>% 
                       filter(HSP90_NA == FALSE), family = "binomial")
glm_dnajc7 <- glmer(outcome ~ log(dnajc7_ngperml) + sex + age_survey + (1|cohort),
                    data = dfbase_noDE %>% 
                        filter(DNAJC7_NA == FALSE), family = "binomial")
glm_hsp70dnajc7 <- glmer(outcome ~ log(dnajc7_ngperml) + log(hsp70_ngperml) + 
                             sex + age_survey + (1|cohort),
                         data = dfbase_noDE %>% 
                             filter(HSP70DNAJC7_NA == FALSE), family = "binomial")


# Identify outliers greater 4 SD on log scale
idxhsp70 <- c(abs(scale(log(dfbase_noDE$hsp70_ngperml))) > 4 & !is.na(log(dfbase_noDE$hsp70_ngperml)) )
idxhsp90 <- c(abs(scale(log(dfbase_noDE$hsp90_ngperml))) > 4 & !is.na(log(dfbase_noDE$hsp90_ngperml)) )
idxdnajc7 <- c(abs(scale(log(dfbase_noDE$dnajc7_ngperml))) > 4 & !is.na(log(dfbase_noDE$dnajc7_ngperml)) )

dfbase_noDE$hsp70_ngperml[ idxhsp70 ]
dfbase_noDE$hsp90_ngperml[ idxhsp90 ] 
dfbase_noDE$dnajc7_ngperml[ idxdnajc7 ]

# Fit sensitivity analysis
glm_hsp70dnajc7_sensitivity <- glmer(outcome ~ log(dnajc7_ngperml) + log(hsp70_ngperml) + 
                             sex + age_survey + (1|cohort),
                         data = dfbase_noDE %>% 
                             filter(HSP70DNAJC7_NA == FALSE & !idxhsp70 & !idxdnajc7),
                         family = "binomial")



summary_glms <- bind_rows(data.frame(model = "hsp70",
                                     nobs = glance(glm_hsp70)$nobs,
                                     ncontrols = sum(glm_hsp70@resp$y ==0),
                                     npatients = sum(glm_hsp70@resp$y ==1),
                                     tidy(glm_hsp70, conf.int = TRUE,
                                          exponentiate = TRUE)),
                          data.frame(model = "hsp90",
                                     nobs = glance(glm_hsp90)$nobs,
                                     ncontrols = sum(glm_hsp90@resp$y ==0),
                                     npatients = sum(glm_hsp90@resp$y ==1),
                                     tidy(glm_hsp90, conf.int = TRUE,
                                          exponentiate = TRUE)),
                          data.frame(model = "dnajc7",
                                     nobs = glance(glm_dnajc7)$nobs,
                                     ncontrols = sum(glm_dnajc7@resp$y ==0),
                                     npatients = sum(glm_dnajc7@resp$y ==1),
                                     tidy(glm_dnajc7, conf.int = TRUE,
                                          exponentiate = TRUE)),
                          data.frame(model = "hsp70 + dnajc7",
                                     nobs = glance(glm_hsp70dnajc7)$nobs,
                                     ncontrols = sum(glm_hsp70dnajc7@resp$y ==0),
                                     npatients = sum(glm_hsp70dnajc7@resp$y ==1),
                                     tidy(glm_hsp70dnajc7, conf.int = TRUE,
                                          exponentiate = TRUE)),
                          data.frame(model = "hsp70 + dnajc7 sensitvitiy",
                                     nobs = glance(glm_hsp70dnajc7_sensitivity)$nobs,
                                     ncontrols = sum(glm_hsp70dnajc7_sensitivity@resp$y ==0),
                                     npatients = sum(glm_hsp70dnajc7_sensitivity@resp$y ==1),
                                     tidy(glm_hsp70dnajc7_sensitivity, conf.int = TRUE,
                                          exponentiate = TRUE)))
# Write results to csv
write.csv(summary_glms, "Results/table2_glm_results.csv")




#### Fit Cox PH models ####
dfbase_noDE$cohort <- as.character(dfbase_noDE$cohort)
cox1 <- coxph(Surv(time_sample2censor_mnths, survoutcome) ~ cohort + log(hsp70_ngperml) +
                  age_onset + dx_delay + site_onset,
              dfbase_noDE)
cox2 <- coxph(Surv(time_sample2censor_mnths, survoutcome) ~ cohort + log(hsp90_ngperml) +
                  age_onset + dx_delay + site_onset,
              dfbase_noDE)
cox3 <- coxph(Surv(time_sample2censor_mnths, survoutcome) ~ cohort + log(dnajc7_ngperml) +
                  age_onset + dx_delay + site_onset,
              dfbase_noDE)
cox4 <- coxph(Surv(time_sample2censor_mnths, survoutcome) ~ cohort + log(hsp70_ngperml) + log(dnajc7_ngperml) +
                  age_onset + dx_delay + site_onset,
              dfbase_noDE)

coxres <- rbind(data.frame(model = 1, nobs = cox1$n, nevents = cox1$nevent, tidy(cox1, conf.int = TRUE, exponentiate = TRUE)),
                data.frame(model = 2, nobs = cox2$n, nevents = cox2$nevent, tidy(cox2, conf.int = TRUE, exponentiate = TRUE)),
                data.frame(model = 3, nobs = cox3$n, nevents = cox3$nevent, tidy(cox3, conf.int = TRUE, exponentiate = TRUE)),
                data.frame(model = 4, nobs = cox4$n, nevents = cox4$nevent, tidy(cox4, conf.int = TRUE, exponentiate = TRUE)))

write.csv(coxres, "Results/table2_coxPH_results.csv", row.names = FALSE)

#### End of analysis ####






