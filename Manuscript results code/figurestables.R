# ------------------------------------------------------------------------------------------------------
#
# Code to reproduce results from primary tables and figures in the manuscript:
# Personalized azithromycin treatment rules for children with watery diarrhea using machine learning
#
#
# Sara Kim
# sara.kim2@emory.edu
# October 27, 2024
#
# --------------------------------------------------------------------------------------------------------

# ----- Load libraries
library(ggplot)
library(patchwork)
library(drotr)


# ----- Load dataset
abcd_data <- read.csv("abcd_otr_deidentified.csv")[,-1] 


# ----- Figure 1

# --- Estimate treatment effects, averaged over five seeds 
# Shigella rule: day 3 diarrhea, day 90 re-hospitalization/death, laz
list_shig_diar <- list(diar_shig150, diar_shig300, diar_shig450, diar_shig600, diar_shig750)
average_across_seeds(list_shig_diar, -0.07)

list_shig_hosp <- list(hosp_shig150, hosp_shig300, hosp_shig450, hosp_shig600, hosp_shig750)
average_across_seeds(list_shig_hosp, -0.02)

list_shig_laz <- list(laz_shig150, laz_shig300, laz_shig450, laz_shig600, laz_shig750)
average_across_seeds(list_shig_laz, 0.06)

# Rotavirus rule
list_rota_diar <- list(diar_rota150, diar_rota300, diar_rota450, diar_rota600, diar_rota750)
average_across_seeds(list_rota_diar, -0.07)

list_rota_hosp <- list(hosp_rota150, hosp_rota300, hosp_rota450, hosp_rota600, hosp_rota750)
average_across_seeds(list_rota_hosp, -0.02)

list_rota_laz <- list(laz_rota150, laz_rota300, laz_rota450, laz_rota600, laz_rota750)
average_across_seeds(list_rota_laz, 0.06)

# pathogen quantities
list_path_diar <- list(diar_path150, diar_path300, diar_path450, diar_path600, diar_path750)
average_across_seeds(list_path_diar, -0.07)

list_path_hosp <- list(hosp_path150, hosp_path300, hosp_path450, hosp_path600, hosp_path750)
average_across_seeds(list_path_hosp, -0.02)

list_path_laz <- list(laz_path150, laz_path300, laz_path450, laz_path600, laz_path750)
average_across_seeds(list_path_laz, 0.06)

# symptoms
list_clin_diar <- list(diar_clin150, diar_clin300, diar_clin450, diar_clin600, diar_clin750)
average_across_seeds(list_clin_diar, -0.07)

list_clin_hosp <- list(hosp_clin150, hosp_clin300, hosp_clin450, hosp_clin600, hosp_clin750)
average_across_seeds(list_clin_hosp, -0.02)

list_clin_laz <- list(laz_clin150, laz_clin300, laz_clin450, laz_clin600, laz_clin750)
average_across_seeds(list_clin_laz, 0.06)

# pathogens and symptoms
list_pc_diar <- list(diar_pc150, diar_pc300, diar_pc450, diar_pc600, diar_pc750)
average_across_seeds(list_pc_diar, -0.07)

list_pc_hosp <- list(hosp_pc150, hosp_pc300, hosp_pc450, hosp_pc600, hosp_pc750)
average_across_seeds(list_pc_hosp, -0.02)

list_pc_laz <- list(laz_pc150, laz_pc300, laz_pc450, laz_pc600, laz_pc750)
average_across_seeds(list_pc_laz, 0.06)

# malnutrition and demographics
list_mal_diar <- list(diar_mal150, diar_mal300, diar_mal450, diar_mal600, diar_mal750)
average_across_seeds(list_mal_diar, -0.07)

list_mal_hosp <- list(hosp_mal150, hosp_mal300, hosp_mal450, hosp_mal600, hosp_mal750)
average_across_seeds(list_mal_hosp, -0.02)

list_mal_laz <- list(laz_mal150, laz_mal300, laz_mal450, laz_mal600, laz_mal750)
average_across_seeds(list_mal_laz, 0.06)

# no pathogen 
list_nop_diar <- list(diar_nop150, diar_nop300, diar_nop450, diar_nop600, diar_nop750)
average_across_seeds(list_nop_diar, -0.07)

list_nop_hosp <- list(hosp_nop150, hosp_nop300, hosp_nop450, hosp_nop600, hosp_nop750)
average_across_seeds(list_nop_hosp, -0.02)

list_nop_laz <- list(laz_nop150, laz_nop300, laz_nop450, laz_nop600, laz_nop750)
average_across_seeds(list_nop_laz, 0.06)

# all information 
list_all_diar <- list(diar_all150, diar_all300, diar_all450, diar_all600, diar_all750)
average_across_seeds(list_all_diar, -0.07)

list_all_hosp <- list(all02_150, all02_300, all02_450, all02_600, all02_750)
average_across_seeds(list_all_hosp, -0.02)

list_all_laz <- list(laz_all150, laz_all300, laz_all450, laz_all600, laz_all750)
average_across_seeds(list_all_laz, 0.06)

# save outcomes 
te_data <- read.csv("te_data.csv")

# -- day 90 hospitalization/death
te_data$Rule_hosp <- as.factor(te_data$Rule_hosp)
te_data$Rule_hosp <- factor(te_data$Rule_hosp, levels=c("Comprehensive (38.0%)","Shigella (8.4%)","Rotavirus (13.6%)","Pathogen Quantities (30.5%)",
                                                        "Symptoms (24.1%)","Pathogen + Symptoms (30.7%)","Host (41.5%)",
                                                        "Host + Symptoms (41.7%)"))

comp_hosp <- ggplot(te_data, aes(x = Subcategory, y = te_hosp, color = Rule_hosp, group = Rule_hosp)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_hosp, ymax = upper_hosp), width = 0.2) +
  geom_line(aes(group = interaction(Rule_hosp, rep(1:3, each = 2))), position = position_dodge(width = 0.5), linetype=2, color="darkgrey") +
  facet_wrap(~ Rule_hosp, scales = "free_x", ncol = 4) +
  theme_minimal() +
  scale_color_manual(values = c("Comprehensive (38.0%)" = "black","Shigella (8.4%)" = "black", "Rotavirus (13.6%)" = "black", 
                                "Pathogen Quantities (30.5%)" = "black",
                                "Symptoms (24.1%)" = "black", "Pathogen + Symptoms (30.7%)" = "black", "Host (41.5%)" = "black",
                                "Host + Symptoms (41.7%)" = "black")) +
  labs(y = "Average Risk Difference with 95% CI") +
  ggtitle("B: Day 90 Re-hospitalization or Death Outcome") +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

# -- day 3 diarrhea
te_data$Rule_diar <- as.factor(te_data$Rule_diar)
te_data$Rule_diar <- factor(te_data$Rule_diar, levels=c("Comprehensive (27.8%)","Shigella (15.9%)","Rotavirus (32.7%)","Pathogen Quantities (25.3%)",
                                                        "Symptoms (18.7%)","Pathogen + Symptoms (27.3%)","Host (23.8%)",
                                                        "Host + Symptoms (25.8%)"))

comp_diar <- ggplot(te_data, aes(x = Subcategory, y = te_diar, color = Rule_diar, group = Rule_diar)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_diar, ymax = upper_diar), width = 0.2) +
  geom_line(aes(group = interaction(Rule_diar, rep(1:3, each = 2))), position = position_dodge(width = 0.5), linetype=2, color="darkgrey") +
  facet_wrap(~ Rule_diar, scales = "free_x", ncol = 4) +
  theme_minimal() +
  scale_color_manual(values = c("Comprehensive (27.8%)" = "black","Shigella (15.9%)" = "black", "Rotavirus (32.7%)" = "black", 
                                "Pathogen Quantities (25.3%)" = "black",
                                "Symptoms (18.7%)" = "black", "Pathogen + Symptoms (27.3%)" = "black", "Host (23.8%)" = "black",
                                "Host + Symptoms (25.8%)" = "black")) +
  labs(y = "Average Risk Difference with 95% CI") +
  ggtitle("A: Day 3 Diarrhea Outcome") +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

# -- lazdiff
te_data$Rule_laz <- as.factor(te_data$Rule_laz)
te_data$Rule_laz <- factor(te_data$Rule_laz, levels=c("Comprehensive (37.6%)","Shigella (12.0%)","Rotavirus (13.0%)","Pathogen Quantities (41.7%)",
                                                      "Symptoms (34.5%)","Pathogen + Symptoms (34.9%)","Host (35.8%)",
                                                      "Host + Symptoms (39.2%)"))

comp_laz <- ggplot(te_data, aes(x = Subcategory, y = te_laz, color = Rule_laz, group = Rule_laz)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_laz, ymax = upper_laz), width = 0.2) +
  geom_line(aes(group = interaction(Rule_laz, rep(1:3, each = 2))), position = position_dodge(width = 0.5), linetype=2, color="darkgrey") +
  facet_wrap(~ Rule_laz, scales = "free_x", ncol = 4) +
  theme_minimal() +
  scale_color_manual(values = c("Comprehensive (37.6%)" = "black", "Shigella (12.0%)" = "black", "Rotavirus (13.0%)" = "black", 
                                "Pathogen Quantities (41.7%)" = "black",
                                "Symptoms (34.5%)" = "black", "Pathogen + Symptoms (34.9%)" = "black", "Host (35.8%)" = "black",
                                "Host + Symptoms (39.2%)" = "black")) +
  labs(y = "Average Change in LAZ with 95% CI") +
  ggtitle("C: Change in Linear Growth Outcome") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  labs(caption = "Abbreviation: Confidence interval (CI); Length-for-Age Z-Score (LAZ)") +
  theme(plot.caption = element_text(hjust = 0, vjust = 1, size = 8),
        plot.caption.position = "plot")

comp_diar + comp_hosp + comp_laz +
  plot_layout(ncol=1)




# ----- Figure 2
# -- day 90 hospitalization/death
# calculate ccc then incorporate into figure
ccc <- epiR::epi.ccc(abcd_data$CATE_pred_shigmed_hosp, abcd_data$CATE_pred_allmed_hosp, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_shig_hosp <- ggplot(abcd_data, aes(x=CATE_pred_shigmed_hosp, y=CATE_pred_allmed_hosp)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Shigella", y = "Comprehensive", subtitle = "CCC: 0.12 (95% CI: 0.10, 0.13)",
       title = "Panel B: Day 90 Re-hospitalization or Death Outcome") +
  scale_color_brewer(palette="Dark2") + theme_minimal() +
  ylim(-0.2, 0.2) + xlim(-0.2,0.2) +
  theme(plot.title = element_text(size = 12))

ccc <- epiR::epi.ccc(abcd_data$CATE_pred_rotamed_hosp, abcd_data$CATE_pred_allmed_hosp, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_rota_hosp <- ggplot(abcd_data, aes(x=CATE_pred_rotamed_hosp, y=CATE_pred_allmed_hosp)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Rotavirus", y = "Comprehensive", subtitle = "CCC: 0.08 (95% CI: 0.07, 0.10)") +
  scale_color_brewer(palette="Dark2") + theme_minimal()+
  ylim(-0.2, 0.2) + xlim(-0.2,0.2)+
  theme(plot.title = element_text(size = 10))

ccc <- epiR::epi.ccc(abcd_data$CATE_pred_pathmed_hosp, abcd_data$CATE_pred_allmed_hosp, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_path_hosp <- ggplot(abcd_data, aes(x=CATE_pred_pathmed_hosp, y=CATE_pred_allmed_hosp)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Pathogen Quantities", y = "Comprehensive", subtitle = "CCC: 0.27 (95% CI: 0.25, 0.29)") +
  scale_color_brewer(palette="Dark2") + theme_minimal()+
  ylim(-0.2, 0.2) + xlim(-0.2,0.2)+
  theme(plot.title = element_text(size = 10))

ccc <- epiR::epi.ccc(abcd_data$CATE_pred_clinmed_hosp, abcd_data$CATE_pred_allmed_hosp, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_clin_hosp <- ggplot(abcd_data, aes(x=CATE_pred_clinmed_hosp, y=CATE_pred_allmed_hosp)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Symptoms", y = "Comprehensive", subtitle="CCC: 0.34 (95% CI: 0.32, 0.36)") +
  scale_color_brewer(palette="Dark2") + theme_minimal()+
  ylim(-0.2, 0.2) + xlim(-0.2,0.2)+
  theme(plot.title = element_text(size = 10))

ccc <- epiR::epi.ccc(abcd_data$CATE_pred_pcmed_hosp, abcd_data$CATE_pred_allmed_hosp, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_pc_hosp <- ggplot(abcd_data, aes(x=CATE_pred_pcmed_hosp, y=CATE_pred_allmed_hosp)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Pathogen + Symptoms", y = "Comprehensive", subtitle = "CCC: 0.43 (95% CI: 0.41, 0.45)") +
  scale_color_brewer(palette="Dark2") + theme_minimal()+
  ylim(-0.2, 0.2) + xlim(-0.2,0.2)+
  theme(plot.title = element_text(size = 10))

ccc <- epiR::epi.ccc(abcd_data$CATE_pred_malmed_hosp, abcd_data$CATE_pred_allmed_hosp, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_mal_hosp <- ggplot(abcd_data, aes(x=CATE_pred_malmed_hosp, y=CATE_pred_allmed_hosp)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Host", y = "Comprehensive", subtitle="CCC: 0.65 (95% CI: 0.63, 0.66)") +
  scale_color_brewer(palette="Dark2") + theme_minimal()+
  ylim(-0.2, 0.2) + xlim(-0.2,0.2)+
  theme(plot.title = element_text(size = 10))

ccc <- epiR::epi.ccc(abcd_data$CATE_pred_nopmed_hosp, abcd_data$CATE_pred_allmed_hosp, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_nop_hosp <- ggplot(abcd_data, aes(x=CATE_pred_nopmed_hosp, y=CATE_pred_allmed_hosp)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Host + Symptoms", y = "Comprehensive", subtitle="CCC: 0.73 (95% CI: 0.72, 0.74)") +
  scale_color_brewer(palette="Dark2") + theme_minimal()+
  ylim(-0.2, 0.2) + xlim(-0.2,0.2)+
  theme(plot.title = element_text(size = 10))

all_shig_hosp + all_rota_hosp + all_path_hosp + all_clin_hosp +
  all_pc_hosp + all_mal_hosp + all_nop_hosp + 
  plot_layout(ncol=4)

# -- day 3 diarrhea
ccc <- epiR::epi.ccc(abcd_data$CATE_pred_shigmed_diar, abcd_data$CATE_pred_allmed_diar, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_shig_diar <- ggplot(abcd_data, aes(x=CATE_pred_shigmed_diar, y=CATE_pred_allmed_diar)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Shigella", y = "Comprehensive", subtitle="CCC: 0.46 (95% CI: 0.44, 0.48)",
       title="Panel A: Day 3 Diarrhea Outcome") +
  scale_color_brewer(palette="Dark2") + theme_minimal() +
  ylim(-0.2, 0.2) + xlim(-0.2,0.2)+
  theme(plot.title = element_text(size = 12))

ccc <- epiR::epi.ccc(abcd_data$CATE_pred_rotamed_diar, abcd_data$CATE_pred_allmed_diar, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_rota_diar <- ggplot(abcd_data, aes(x=CATE_pred_rotamed_diar, y=CATE_pred_allmed_diar)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Rotavirus", y = "Comprehensive", subtitle="CCC: 0.56 (95% CI: 0.54, 0.57)") +
  scale_color_brewer(palette="Dark2") + theme_minimal()+
  ylim(-0.2, 0.2) + xlim(-0.2,0.2)+
  theme(plot.title = element_text(size = 10))

ccc <- epiR::epi.ccc(abcd_data$CATE_pred_pathmed_diar, abcd_data$CATE_pred_allmed_diar, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_path_diar <- ggplot(abcd_data, aes(x=CATE_pred_pathmed_diar, y=CATE_pred_allmed_diar)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Pathogen Quantities", y = "Comprehensive", subtitle="CCC: 0.76 (95% CI: 0.75, 0.77)") +
  scale_color_brewer(palette="Dark2") + theme_minimal()+
  ylim(-0.2, 0.2) + xlim(-0.2,0.2)+
  theme(plot.title = element_text(size = 10))

ccc <- epiR::epi.ccc(abcd_data$CATE_pred_clinmed_diar, abcd_data$CATE_pred_allmed_diar, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_clin_diar <- ggplot(abcd_data, aes(x=CATE_pred_clinmed_diar, y=CATE_pred_allmed_diar)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Symptoms", y = "Comprehensive", subtitle="CCC: 0.30 (95% CI: 0.27, 0.31)") +
  scale_color_brewer(palette="Dark2") + theme_minimal()+
  ylim(-0.2, 0.2) + xlim(-0.2,0.2)+
  theme(plot.title = element_text(size = 10))

ccc <- epiR::epi.ccc(abcd_data$CATE_pred_pcmed_diar, abcd_data$CATE_pred_allmed_diar, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_pc_diar <- ggplot(abcd_data, aes(x=CATE_pred_pcmed_diar, y=CATE_pred_allmed_diar)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Pathogen + Symptoms", y = "Comprehensive", subtitle="CCC: 0.82 (95% CI: 0.81, 0.83)") +
  scale_color_brewer(palette="Dark2") + theme_minimal()+
  ylim(-0.2, 0.2) + xlim(-0.2,0.2)+
  theme(plot.title = element_text(size = 10))

ccc <- epiR::epi.ccc(abcd_data$CATE_pred_malmed_diar, abcd_data$CATE_pred_allmed_diar, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_mal_diar <- ggplot(abcd_data, aes(x=CATE_pred_malmed_diar, y=CATE_pred_allmed_diar)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Host", y = "Comprehensive", subtitle="CCC: 0.29 (95% CI: 0.26, 0.31)") +
  scale_color_brewer(palette="Dark2") + theme_minimal()+
  ylim(-0.2, 0.2) + xlim(-0.2,0.2)+
  theme(plot.title = element_text(size = 10))

ccc <- epiR::epi.ccc(abcd_data$CATE_pred_nopmed_diar, abcd_data$CATE_pred_allmed_diar, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_nop_diar <- ggplot(abcd_data, aes(x=CATE_pred_nopmed_diar, y=CATE_pred_allmed_diar)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Host + Symptoms", y = "Comprehensive", subtitle="CCC: 0.43 (95% CI: 0.41, 0.45)") +
  scale_color_brewer(palette="Dark2") + theme_minimal()+
  ylim(-0.2, 0.2) + xlim(-0.2,0.2)+
  theme(plot.title = element_text(size = 10))

all_shig_diar + all_rota_diar + all_path_diar + all_clin_diar +
  all_pc_diar + all_mal_diar + all_nop_diar + 
  plot_layout(ncol=4)


# -- lazdiff
ccc <- epiR::epi.ccc(abcd_data$CATE_pred_shigmed_laz, abcd_data$CATE_pred_allmed_laz, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_shig_laz <- ggplot(abcd_data, aes(x=CATE_pred_shigmed_laz, y=CATE_pred_allmed_laz)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Shigella", 
       y="Comprehensive", 
       title="Panel C: Change in Linear Growth Outcome",
       subtitle="CCC: 1.2e-07 (95% CI: 3.0e-08, 2.2e-07)") +
  labs(caption = "Abbreviation: Confidence interval (CI); Length-for-Age Z-Score (LAZ)") +
  scale_color_brewer(palette="Dark2") + theme_minimal() +
  ylim(-0.2, 0.2) + xlim(-0.2,0.2)+
  theme(plot.title = element_text(size = 12))

ccc <- epiR::epi.ccc(abcd_data$CATE_pred_rotamed_laz, abcd_data$CATE_pred_allmed_laz, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_rota_laz <- ggplot(abcd_data, aes(x=CATE_pred_rotamed_laz, y=CATE_pred_allmed_laz)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Rotavirus", y = "Comprehensive", subtitle="CCC: -2.7e-08 (95% CI: -1.0e-0.7, 4.7e-0.8)") +
  scale_color_brewer(palette="Dark2") + theme_minimal()+
  ylim(-0.2, 0.2) + xlim(-0.2,0.2)+
  theme(plot.title = element_text(size = 10))

ccc <- epiR::epi.ccc(abcd_data$CATE_pred_pathmed_laz, abcd_data$CATE_pred_allmed_laz, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_path_laz <- ggplot(abcd_data, aes(x=CATE_pred_pathmed_laz, y=CATE_pred_allmed_laz)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Pathogen Quantities", y = "Comprehensive", subtitle="CCC: 8.1e-07 (95% CI: 6.7e-07, 9.5e-07)") +
  scale_color_brewer(palette="Dark2") + theme_minimal()+
  ylim(-0.2, 0.2) + xlim(-0.2,0.2)+
  theme(plot.title = element_text(size = 10))

ccc <- epiR::epi.ccc(abcd_data$CATE_pred_clinmed_laz, abcd_data$CATE_pred_allmed_laz, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_clin_laz <- ggplot(abcd_data, aes(x=CATE_pred_clinmed_laz, y=CATE_pred_allmed_laz)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Symptoms", y = "Comprehensive", subtitle="CCC: -2.2e-07 (95% CI: -3.9e-07, -4.9e-08)") +
  scale_color_brewer(palette="Dark2") + theme_minimal()+
  ylim(-0.2, 0.2) + xlim(-0.2,0.2)+
  theme(plot.title = element_text(size = 10))

ccc <- epiR::epi.ccc(abcd_data$CATE_pred_pcmed_laz, abcd_data$CATE_pred_allmed_laz, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_pc_laz <- ggplot(abcd_data, aes(x=CATE_pred_pcmed_laz, y=CATE_pred_allmed_laz)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Pathogen + Symptoms", y = "Comprehensive", subtitle="CCC: 3.0e-07 (95% CI: -4.8e-07, 1.1e-06)") +
  scale_color_brewer(palette="Dark2") + theme_minimal()+
  ylim(-0.2, 0.2) + xlim(-0.2,0.2)+
  theme(plot.title = element_text(size = 10))

ccc <- epiR::epi.ccc(abcd_data$CATE_pred_malmed_laz, abcd_data$CATE_pred_allmed_laz, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_mal_laz <- ggplot(abcd_data, aes(x=CATE_pred_malmed_laz, y=CATE_pred_allmed_laz)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Host", y = "Comprehensive", subtitle="CCC: 6.7e-09 (95% CI: -1.8e-07, 1.9e-06)") +
  scale_color_brewer(palette="Dark2") + theme_minimal()+
  ylim(-0.2, 0.2) + xlim(-0.2,0.2)+
  theme(plot.title = element_text(size = 10))

ccc <- epiR::epi.ccc(abcd_data$CATE_pred_nopmed_laz, abcd_data$CATE_pred_allmed_laz, ci="z-transform", conf.level=0.95)
ccc$rho.c
all_nop_laz <- ggplot(abcd_data, aes(x=CATE_pred_nopmed_laz, y=CATE_pred_allmed_laz)) +
  geom_point(size=1, alpha=0.3, color="black") +
  geom_abline(color="darkgray", linetype=3) +
  labs(x="Host + Symptoms", y = "Comprehensive", subtitle="CCC: -2.0e-07 (95% CI: -4.1e-07, -7.1e-10)") +
  scale_color_brewer(palette="Dark2") + theme_minimal()+
  ylim(-0.2, 0.2) + xlim(-0.2,0.2)+
  theme(plot.title = element_text(size = 10))

all_shig_diar + all_rota_diar + all_path_diar + all_clin_diar + all_pc_diar + all_mal_diar + all_nop_diar +
  all_shig_hosp + all_rota_hosp + all_path_hosp + all_clin_hosp + all_pc_hosp + all_mal_hosp + all_nop_hosp +
  all_shig_laz + all_rota_laz + all_path_laz + all_clin_laz + all_pc_laz + all_mal_laz + all_nop_laz +
  plot_layout(ncol=7)


# ----- Figure 3
threshold_data <- read.csv("threshold_data.csv")

# -- day 90 hospitalization/death
# plot
thresh.hosp1 <- ggplot(data=threshold_data, aes(x=threshold_hosp, y=prop_hosp)) + 
  geom_point(aes(y=prop_hosp), color = "black") + 
  geom_line(aes(y=prop_hosp), color = "black", group=1) +
  geom_ribbon(aes(x = 1:length(threshold_hosp), ymin=prop_lower_hosp, ymax=prop_upper_hosp), 
              alpha=0.1, fill = "darkgrey",  
              color = "darkgrey", linetype = "dotted") +
  theme_minimal() +
  xlab("Threshold (RD)") + ylab("Proportion Treated") +
  ggtitle("B: Day 90 Re-hospitalization or Death Outcome ") +
  scale_y_continuous(labels = scales::percent)

thresh.hosp2 <- ggplot(data=threshold_data, aes(x=threshold_hosp, y=te_hosp)) + 
  geom_point(aes(y=te_hosp), color = "black") + 
  geom_line(aes(y=te_hosp), color = "black", group=1) +
  geom_ribbon(aes(x = 1:length(threshold_hosp), ymin=te_lower_hosp, ymax=te_upper_hosp), 
              alpha=0.1, fill = "darkgrey",  
              color = "darkgrey", linetype = "dotted") +
  theme_minimal() +
  xlab("Threshold (RD)") + ylab("Average RD among Recommended for Treatment")


# -- day 3 diarrhea
threshold_data$threshold_diar <- as.factor(threshold_data$threshold_diar)
threshold_data$threshold_diar <- factor(threshold_data$threshold_diar, levels=c("-2%","-3%","-4%","-5%","-6%","-7%","-8%","-9%","-10%"))

# plot
thresh.diar1 <- ggplot(data=threshold_data, aes(x=threshold_diar, y=prop_diar)) + 
  geom_point(aes(y=prop_diar), color = "black") + 
  geom_line(aes(y=prop_diar), color = "black", group=1) +
  geom_ribbon(aes(x = 1:length(threshold_diar), ymin=prop_lower_diar, ymax=prop_upper_diar), 
              alpha=0.1, fill = "darkgrey",  
              color = "darkgrey", linetype = "dotted") +
  theme_minimal() +
  xlab("Threshold (RD)") + ylab("Proportion Treated") +
  ggtitle("A: Day 3 Diarrhea Outcome") +
  scale_y_continuous(labels = scales::percent)

thresh.diar2 <- ggplot(data=threshold_data, aes(x=threshold_diar, y=te_diar)) + 
  geom_point(aes(y=te_diar), color = "black") + 
  geom_line(aes(y=te_diar), color = "black", group=1) +
  geom_ribbon(aes(x = 1:length(threshold_diar), ymin=te_lower_diar, ymax=te_upper_diar), 
              alpha=0.1, fill = "darkgrey",  
              color = "darkgrey", linetype = "dotted") +
  theme_minimal() +
  xlab("Threshold (RD)") + ylab("Average RD among Recommended for Treatment")

# -- lazdiff

threshold_data$threshold_laz <- as.factor(threshold_data$threshold_laz)

# plot
thresh.laz1 <- ggplot(data=threshold_data, aes(x=threshold_laz, y=prop_laz)) + 
  geom_point(aes(y=prop_laz), color = "black") + 
  geom_line(aes(y=prop_laz), color = "black", group=1) +
  geom_ribbon(aes(x = 1:length(threshold_laz), ymin=prop_lower_laz, ymax=prop_upper_laz), 
              alpha=0.1, fill = "darkgrey",  
              color = "darkgrey", linetype = "dotted") +
  theme_minimal() +
  xlab("Threshold (Change in LAZ)") + ylab("Proportion Treated") +
  ggtitle("C: Change in Linear Growth Outcome") +
  scale_y_continuous(labels = scales::percent) +
  labs(caption = "Abbreviation: Risk difference (RD); Length-for-Age Z-Score (LAZ)") +
  theme(plot.caption = element_text(hjust = 0, vjust = 1, size = 8),
        plot.caption.position = "plot")

thresh.laz2 <- ggplot(data=threshold_data, aes(x=threshold_laz, y=te_laz)) + 
  geom_point(aes(y=te_laz), color = "black") + 
  geom_line(aes(y=te_laz), color = "black", group=1) +
  geom_ribbon(aes(x = 1:length(threshold_diar), ymin=te_lower_laz, ymax=te_upper_laz), 
              alpha=0.1, fill = "darkgrey",  
              color = "darkgrey", linetype = "dotted") +
  theme_minimal() +
  xlab("Threshold (Change in LAZ)") + ylab("Average Change in LAZ among Recommended for Treatment")

# thresh.laz1 + thresh.laz2 +
#   plot_layout(ncol=2)

thresh.diar1 + thresh.diar2 + 
  thresh.hosp1 + thresh.hosp2 +
  thresh.laz1 + thresh.laz2 +
  plot_layout(ncol=2)



# ----- Table 1
# day 3 diarrhea
subday3diar <- main_data %>%
  mutate(
    d_bacteria = case_when(d_bacteria == 0 ~ "Unlikely",
                           d_bacteria == 1 ~ "Possible",
                           d_bacteria == 2 ~ "Likley"),
    dy1_ant_sex = case_when(dy1_ant_sex == 1 ~ "Male",
                            dy1_ant_sex == 2 ~ "Female"),
    site = case_when(site == 2 ~ "Bangladesh",
                     site == 3 ~ "Kenya",
                     site == 4 ~ "Malawi",
                     site == 5 ~ "Mali",
                     site == 6 ~ "India",
                     site == 7 ~ "Tanzania",
                     site == 8 ~ "Pakistan"),
    dy1_scrn_sstools = case_when(dy1_scrn_sstools == 0 ~ "0",
                                 dy1_scrn_sstools == 1 ~ "1",
                                 dy1_scrn_sstools == 2 ~ "2", 
                                 dy1_scrn_sstools > 2 ~ ">2"),
    dy1_scrn_dehydr = case_when(dy1_scrn_dehydr == 1 ~ "None",
                                dy1_scrn_dehydr == 2 ~ "Some",
                                dy1_scrn_dehydr == 3 ~ "Severe"),
    dy1_scrn_vomitall = case_when (dy1_scrn_vomitall == 1 ~ "Yes",
                                   dy1_scrn_vomitall == 2 ~ "No")) %>%
  select(decision_allmed_diar, day3diar, d_bacteria, shigdquant, stetecdquant, 
         tepecdquant, choleradquant, campdquant, salmondquant, rotadquant, cryptodquant,
         agemchild, dy1_ant_sex, site, an_ses_quintile, an_tothhlt5,
         wfazscore, wflzscore, lfazscore, avemuac,
         dy1_scrn_lstools, dy1_scrn_sstools, dy1_scrn_diardays, dy1_scrn_dehydr, dy1_scrn_vomitall)
subday3diar %>%
  tbl_summary(by=decision_allmed_diar,
              label = list(
                day3diar ~ "Diarrhea on Day 3",
                d_bacteria ~ "Bacteria Attributed Diarrhea",
                shigdquant ~ "Likely Shigella",
                stetecdquant ~ "Likely ST ETEC",
                tepecdquant ~ "Likely T EPEC",
                choleradquant ~ "Likely V Cholerae",
                campdquant ~ "Likely Campylobacter",
                salmondquant ~ "Likely Salmonella",
                rotadquant ~ "Likely Rotavirus",
                cryptodquant ~ "Likely Cryptosporidium",
                agemchild ~ "Age (Months)",
                dy1_ant_sex ~ "Sex",
                site ~ "Study Site",
                an_ses_quintile ~ "SES Quintile",
                an_tothhlt5 ~ "# <5 Years in HH",
                wfazscore ~ "Weight for Age Z-Score",
                wflzscore ~ "Weight for Length Z-Score",
                lfazscore ~ "Length for Age Z-Score",
                avemuac ~ "Average Middle-Upper Arm Circumference",
                dy1_scrn_lstools ~ "# Loose Stools in 24 Hours",
                dy1_scrn_sstools ~ "# Solid Stools in 24 Hours",
                dy1_scrn_diardays ~ "Illness Duration (Days)",
                dy1_scrn_dehydr ~ "Dehydration Status",
                dy1_scrn_vomitall ~ "Vomit"),
              statistic = list(all_continuous()~"{mean} ({sd})",
                               all_categorical() ~ "{n} ({p}%)"),
              digits = all_continuous() ~ 2
  ) %>%
  modify_caption("**Characteristics by Treatment under All Information Rule**") %>%
  bold_labels() 

# hospitalization/death
subhosp <- abcd_data2 %>%
  mutate(d_bacteria = case_when(d_bacteria == 0 ~ "Unlikely",
                                d_bacteria == 1 ~ "Possible",
                                d_bacteria == 2 ~ "Likley"),
         dy1_ant_sex = case_when(dy1_ant_sex == 1 ~ "Male",
                                 dy1_ant_sex == 2 ~ "Female"),
         site = case_when(site == 2 ~ "Bangladesh",
                          site == 3 ~ "Kenya",
                          site == 4 ~ "Malawi",
                          site == 5 ~ "Mali",
                          site == 6 ~ "India",
                          site == 7 ~ "Tanzania",
                          site == 8 ~ "Pakistan"),
         dy1_scrn_sstools = case_when(dy1_scrn_sstools == 0 ~ "0",
                                      dy1_scrn_sstools == 1 ~ "1",
                                      dy1_scrn_sstools == 2 ~ "2", 
                                      dy1_scrn_sstools > 2 ~ ">2"),
         dy1_scrn_dehydr = case_when(dy1_scrn_dehydr == 1 ~ "None",
                                     dy1_scrn_dehydr == 2 ~ "Some",
                                     dy1_scrn_dehydr == 3 ~ "Severe"),
         dy1_scrn_vomitall = case_when (dy1_scrn_vomitall == 1 ~ "Yes",
                                        dy1_scrn_vomitall == 2 ~ "No")) %>%
  select(decision_allmean_hosp, an_hosp90death, d_bacteria, shigella_new, st_etec_new, 
         tepec_new, v_cholerae_new, campylobacter_new, salmonella_new, rotavirus_new, cryptosporidium_new,
         agemchild, dy1_ant_sex, site, an_ses_quintile, an_tothhlt5,
         wfazscore, wflzscore, lfazscore, avemuac,
         dy1_scrn_lstools, dy1_scrn_sstools, dy1_scrn_diardays, dy1_scrn_dehydr, dy1_scrn_vomitall)
subhosp %>%
  tbl_summary(by=decision_allmean_hosp,
              label = list(decision_allmean_hosp ~ "Treatment under All Information Rule",
                           an_hosp90death ~ "Hospitalization or Death by Day 90",
                           d_bacteria ~ "Bacteria Attributed Diarrhea",
                           shigella_new ~ "Shigella Quantity",
                           st_etec_new ~ "ST ETEC Quantity",
                           tepec_new ~ "T EPEC Quantity",
                           v_cholerae_new ~ "V Cholerae Quantity",
                           campylobacter_new ~ "Campylobacter Quantity",
                           salmonella_new ~ "Salmonella Quantity",
                           rotavirus_new ~ "Rotavirus Quantity",
                           cryptosporidium_new ~ "Cryptosporidium Quantity",
                           agemchild ~ "Age (Months)",
                           dy1_ant_sex ~ "Sex",
                           site ~ "Study Site",
                           an_ses_quintile ~ "SES Quintile",
                           an_tothhlt5 ~ "# <5 Years in HH",
                           wfazscore ~ "Weight for Age Z-Score",
                           wflzscore ~ "Weight for Length Z-Score",
                           lfazscore ~ "Length for Age Z-Score",
                           avemuac ~ "Average Middle-Upper Arm Circumference",
                           dy1_scrn_lstools ~ "# Loose Stools in 24 Hours",
                           dy1_scrn_sstools ~ "# Solid Stools in 24 Hours",
                           dy1_scrn_diardays ~ "Illness Duration (Days)",
                           dy1_scrn_dehydr ~ "Dehydration Status",
                           dy1_scrn_vomitall ~ "Vomit"),
              statistic = list(all_continuous()~"{mean} ({sd})",
                               all_categorical() ~ "{n} ({p}%)"),
              digits = all_continuous() ~ 2
  ) %>%
  modify_caption("**Characteristics by Treatment under All Information Rule**") %>%
  bold_labels() 

# laz
sublaz <- main_data %>%
  mutate(
    d_bacteria = case_when(d_bacteria == 0 ~ "Unlikely",
                           d_bacteria == 1 ~ "Possible",
                           d_bacteria == 2 ~ "Likley"),
    dy1_ant_sex = case_when(dy1_ant_sex == 1 ~ "Male",
                            dy1_ant_sex == 2 ~ "Female"),
    site = case_when(site == 2 ~ "Bangladesh",
                     site == 3 ~ "Kenya",
                     site == 4 ~ "Malawi",
                     site == 5 ~ "Mali",
                     site == 6 ~ "India",
                     site == 7 ~ "Tanzania",
                     site == 8 ~ "Pakistan"),
    dy1_scrn_sstools = case_when(dy1_scrn_sstools == 0 ~ "0",
                                 dy1_scrn_sstools == 1 ~ "1",
                                 dy1_scrn_sstools == 2 ~ "2", 
                                 dy1_scrn_sstools > 2 ~ ">2"),
    dy1_scrn_dehydr = case_when(dy1_scrn_dehydr == 1 ~ "None",
                                dy1_scrn_dehydr == 2 ~ "Some",
                                dy1_scrn_dehydr == 3 ~ "Severe"),
    dy1_scrn_vomitall = case_when (dy1_scrn_vomitall == 1 ~ "Yes",
                                   dy1_scrn_vomitall == 2 ~ "No")) %>%
  select(decision_allmed_laz, day3diar, d_bacteria, shigdquant, stetecdquant, 
         tepecdquant, choleradquant, campdquant, salmondquant, rotadquant, cryptodquant,
         agemchild, dy1_ant_sex, site, an_ses_quintile, an_tothhlt5,
         wfazscore, wflzscore, lfazscore, avemuac,
         dy1_scrn_lstools, dy1_scrn_sstools, dy1_scrn_diardays, dy1_scrn_dehydr, dy1_scrn_vomitall)
sublaz %>%
  tbl_summary(by=decision_allmed_laz,
              label = list(
                day3diar ~ "Diarrhea on Day 3",
                d_bacteria ~ "Bacteria Attributed Diarrhea",
                shigdquant ~ "Likely Shigella",
                stetecdquant ~ "Likely ST ETEC",
                tepecdquant ~ "Likely T EPEC",
                choleradquant ~ "Likely V Cholerae",
                campdquant ~ "Likely Campylobacter",
                salmondquant ~ "Likely Salmonella",
                rotadquant ~ "Likely Rotavirus",
                cryptodquant ~ "Likely Cryptosporidium",
                agemchild ~ "Age (Months)",
                dy1_ant_sex ~ "Sex",
                site ~ "Study Site",
                an_ses_quintile ~ "SES Quintile",
                an_tothhlt5 ~ "# <5 Years in HH",
                wfazscore ~ "Weight for Age Z-Score",
                wflzscore ~ "Weight for Length Z-Score",
                lfazscore ~ "Length for Age Z-Score",
                avemuac ~ "Average Middle-Upper Arm Circumference",
                dy1_scrn_lstools ~ "# Loose Stools in 24 Hours",
                dy1_scrn_sstools ~ "# Solid Stools in 24 Hours",
                dy1_scrn_diardays ~ "Illness Duration (Days)",
                dy1_scrn_dehydr ~ "Dehydration Status",
                dy1_scrn_vomitall ~ "Vomit"),
              statistic = list(all_continuous()~"{mean} ({sd})",
                               all_categorical() ~ "{n} ({p}%)"),
              digits = all_continuous() ~ 2
  ) %>%
  modify_caption("**Characteristics by Treatment under All Information Rule**") %>%
  bold_labels() 


# ----- Table 2
# get sensitivity and specificity values day 3 diarrhea
abcd_data %>%
  tbl_cross(row = decision_shigmed_diar, col = decision_allmed_diar, percent = "col",
            label = list(decision_shigmed_diar ~ "Shigella Quantity Rule",
                         decision_allmed_diar ~ "All Information Rule")) %>%
  bold_labels()
abcd_data %>%
  tbl_cross(row = decision_rotamed_diar, col = decision_allmed_diar, percent = "col",
            label = list(decision_rotamed_diar ~ "Rotavirus Quantity Rule",
                         decision_allmed_diar ~ "All Information Rule")) %>%
  bold_labels()
abcd_data %>%
  tbl_cross(row = decision_pathmed_diar, col = decision_allmed_diar, percent = "col",
            label = list(decision_pathmed_diar ~ "Pathogen Quantity Rule",
                         decision_allmed_diar ~ "All Information Rule")) %>%
  bold_labels()
abcd_data %>%
  tbl_cross(row = decision_clinmed_diar, col = decision_allmed_diar, percent = "col",
            label = list(decision_clinmed_diar ~ "Symptoms Rule",
                         decision_allmed_diar ~ "All Information Rule")) %>%
  bold_labels()
abcd_data %>%
  tbl_cross(row = decision_pcmed_diar, col = decision_allmed_diar, percent = "col",
            label = list(decision_pcmed_diar ~ "PC Rule",
                         decision_allmed_diar ~ "All Information Rule")) %>%
  bold_labels()
abcd_data %>%
  tbl_cross(row = decision_malmed_diar, col = decision_allmed_diar, percent = "col",
            label = list(decision_malmed_diar ~ "host Rule",
                         decision_allmed_diar ~ "All Information Rule")) %>%
  bold_labels()
abcd_data %>%
  tbl_cross(row = decision_nopmed_diar, col = decision_allmed_diar, percent = "col",
            label = list(decision_nopmed_diar ~ "no pathogen Rule",
                         decision_allmed_diar ~ "All Information Rule")) %>%
  bold_labels()


# get sensitivity and specificity values day 90 hospitalization/death
abcd_data %>%
  tbl_cross(row = decision_shigmed_hosp, col = decision_allmed_hosp, percent = "col",
            label = list(decision_shigmed_hosp ~ "Shigella Quantity Rule",
                         decision_allmed_hosp ~ "All Information Rule")) %>%
  bold_labels()
abcd_data %>%
  tbl_cross(row = decision_rotamed_hosp, col = decision_allmed_hosp, percent = "col",
            label = list(decision_rotamed_hosp ~ "Rotavirus Quantity Rule",
                         decision_allmed_hosp ~ "All Information Rule")) %>%
  bold_labels()
abcd_data %>%
  tbl_cross(row = decision_pathmed_hosp, col = decision_allmed_hosp, percent = "col",
            label = list(decision_pathmed_hosp ~ "Pathogen Quantity Rule",
                         decision_allmed_hosp ~ "All Information Rule")) %>%
  bold_labels()
abcd_data %>%
  tbl_cross(row = decision_clinmed_hosp, col = decision_allmed_hosp, percent = "col",
            label = list(decision_clinmed_hosp ~ "Symptoms Rule",
                         decision_allmed_hosp ~ "All Information Rule")) %>%
  bold_labels()
abcd_data %>%
  tbl_cross(row = decision_pcmed_hosp, col = decision_allmed_hosp, percent = "col",
            label = list(decision_pcmed_hosp ~ "PC Rule",
                         decision_allmed_hosp ~ "All Information Rule")) %>%
  bold_labels()
abcd_data %>%
  tbl_cross(row = decision_malmed_hosp, col = decision_allmed_hosp, percent = "col",
            label = list(decision_malmed_hosp ~ "host Rule",
                         decision_allmed_hosp ~ "All Information Rule")) %>%
  bold_labels()
abcd_data %>%
  tbl_cross(row = decision_nopmed_hosp, col = decision_allmed_hosp, percent = "col",
            label = list(decision_nopmed_hosp ~ "no pathogen Rule",
                         decision_allmed_hosp ~ "All Information Rule")) %>%
  bold_labels()


# get sensitivity and specificity values lazdiff
abcd_data %>%
  tbl_cross(row = decision_shigmed_laz, col = decision_allmed_laz, percent = "col",
            label = list(decision_shigmed_laz ~ "Shigella Quantity Rule",
                         decision_allmed_laz ~ "All Information Rule")) %>%
  bold_labels()
abcd_data %>%
  tbl_cross(row = decision_rotamed_laz, col = decision_allmed_laz, percent = "col",
            label = list(decision_rotamed_laz ~ "Rotavirus Quantity Rule",
                         decision_allmed_laz ~ "All Information Rule")) %>%
  bold_labels()
abcd_data %>%
  tbl_cross(row = decision_pathmed_laz, col = decision_allmed_laz, percent = "col",
            label = list(decision_pathmed_laz ~ "Pathogen Quantity Rule",
                         decision_allmed_laz ~ "All Information Rule")) %>%
  bold_labels()
abcd_data %>%
  tbl_cross(row = decision_clinmed_laz, col = decision_allmed_laz, percent = "col",
            label = list(decision_clinmed_laz ~ "Symptoms Rule",
                         decision_allmed_laz ~ "All Information Rule")) %>%
  bold_labels()
abcd_data %>%
  tbl_cross(row = decision_pcmed_laz, col = decision_allmed_laz, percent = "col",
            label = list(decision_pcmed_laz ~ "PC Rule",
                         decision_allmed_laz ~ "All Information Rule")) %>%
  bold_labels()
abcd_data %>%
  tbl_cross(row = decision_malmed_laz, col = decision_allmed_laz, percent = "col",
            label = list(decision_malmed_laz ~ "host Rule",
                         decision_allmed_laz ~ "All Information Rule")) %>%
  bold_labels()
abcd_data %>%
  tbl_cross(row = decision_nopmed_laz, col = decision_allmed_laz, percent = "col",
            label = list(decision_nopmed_laz ~ "no pathogen Rule",
                         decision_allmed_laz ~ "All Information Rule")) %>%
  bold_labels()
