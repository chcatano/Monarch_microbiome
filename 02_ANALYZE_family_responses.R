#----
# R script to analyze caterpillar microbiota variation for the Immigration Experiment
# Creator: Christopher P. Catano, chcatano@gmail.com
# Data input: "family_diveristy_data_final.csv"
# Data output:
#----


# Check for and install/load required packages
rm(list = ls())

for (package in c('DHARMa', 'effects', 'emmeans', 'gamlss', 'ggeffects', 'ggpubr',  
                  'glmmTMB', 'grid', 'lme4', 'report', 'sjPlot', 'splines', 
                  'tidyverse')) {
  if (!require(package, character.only = T, quietly = T)) {
    install.packages(package)
    library(package, character.only = T)
  }
}

# Set graphics parameters
theme_set(theme_bw() +  
            theme(legend.position = "bottom", panel.grid.minor = element_blank(), 
                  panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5), 
                  text = element_text(size = 10), axis.text = element_text(size = 8)))

#### 1. Import & prep data ####
# data set created in 01_CALCULATE_community_response.R
my.data <- read.csv("data/family_diveristy_data_final.csv", row.names = 1)
str(my.data)

# Make immigration treatments continuous for modeling and replace with number (0-5)
# corresponding to the microbe density
my.data$Immigration <- my.data$Treatment
my.data <- my.data %>%
  filter(Immigration != "C") %>%
  mutate(Immigration = replace(Immigration, Immigration == "D0", 5)) %>%
  mutate(Immigration = replace(Immigration, Immigration == "D1", 4)) %>%
  mutate(Immigration = replace(Immigration, Immigration == "D2", 3)) %>%
  mutate(Immigration = replace(Immigration, Immigration == "D3", 2)) %>%
  mutate(Immigration = replace(Immigration, Immigration == "D4", 1)) %>%
  mutate(Immigration = replace(Immigration, Immigration == "D5", 0)) %>%
  mutate(across(Immigration, as.numeric)) %>%
  mutate(across(Cohort, ~ factor(., levels = c("3rd_instar", "5th_instar"))))

# standardize immigration variable for interpretability and convergence in GLMMs
my.data$std_Immigration <- scale(my.data$Immigration)


#### 2. Model Richness ####
# fit models with different random effects structures and choice the best fit
# based on model convergence and model support (AIC or BIC).

# Random slopes and intercepts ("boundary (singular)")
mod.alpha0.rs <- lmer(focal.alpha.0 ~ Immigration * Cohort + 
                        (1 + std_Immigration | Plant_ID) + 
                        (1 + std_Immigration | Lineage), 
                      data = my.data)

# Random slopes and intercept for Lineage, random intercept for Plant ("boundary (singular)")
mod.alpha0.rsL <- lmer(focal.alpha.0 ~ Immigration * Cohort + 
                         (1 | Plant_ID) + 
                         (1 + std_Immigration | Lineage), 
                       data = my.data)

# Random slopes and intercept for Plant, random intercept for Lineage ("boundary (singular)")
mod.alpha0.rsP <- lmer(focal.alpha.0 ~ Immigration * Cohort + 
                         (1 + std_Immigration | Plant_ID) + 
                         (1 | Lineage), 
                       data = my.data)

# Random intercepts
mod.alpha0.ri <- lmer(focal.alpha.0 ~ Immigration * Cohort + 
                        (1 | Plant_ID) + 
                        (1 | Lineage), 
                      data = my.data)

# Random intercepts only model fit best
anova(mod.alpha0.rs, mod.alpha0.rsL, mod.alpha0.rsP, mod.alpha0.ri)

# Check diagnostics: residual plots created through simulation-based approach
# in DHARMa package (all diagnostics check out)
simulationOutput <- simulateResiduals(mod.alpha0.ri, plot = F, use.u = T, n = 1000)
(diag.beta <- plot(simulationOutput))

# Get immigration effect slope for each cohort
emt <- emtrends(mod.alpha0.ri, pairwise ~ Cohort, var = "Immigration")
summary(emt, infer = c(TRUE, TRUE))

# model results table
plot_model(mod.alpha0.ri, show.values = TRUE, type = "int", show.data = T, jitter = 0.2)
tab_model(mod.alpha0.ri, show.se = TRUE, show.std = TRUE)

# extract predictions (estimated marginal means) and confidence intervals for parameters
dat <- ggpredict(mod.alpha0.ri, terms = c("Immigration", "Cohort"))
moderator_values <- sort(as.character(unique(dat$group)))

# change facet names
cohort_names <- c('3rd_instar' = "3rd instar", '5th_instar' = "5th instar")

# Plot richness response
(fam_rich_plot <- ggplot(as.data.frame(dat), 
                         aes(x = x, y = predicted, group = group, color = as.character(group)), 
                         fill = as.character(group)) + 
    geom_point(data = attr(dat, "rawdata"), aes(x = jitter(x), y = response, alpha = 0.6)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), color = FALSE, alpha = 0.2) +
    geom_line() +
    facet_wrap(~ group, labeller = as_labeller(cohort_names)) +
    scale_color_discrete(guide = "legend", aesthetics = c("colour", "fill"), 
                         breaks = moderator_values, 
                         labels = paste(format(moderator_values), 
                                        c("3rd instar", "5th instar"))) +
    scale_color_manual(values = c("darkorange", "blue")) +
    scale_fill_manual(values = alpha(c("darkorange", "blue"), .2))+
    ylab(expression(paste(alpha, "-diversity"))) +
    scale_x_continuous(name = expression(paste("Immigration density (1.7 CFU/", mu, "L)")),
                       breaks = c(0, 1, 2, 3, 4, 5),
                       labels = c('0' = '10^0', '1' = '10^1', 
                                  '2' = '10^2', '3' = '10^3', 
                                  '4' = '10^4', '5' = '10^5')) +
    theme(legend.position = "none") 
)


#### 3. Model Community size (J) ####

# Random slopes and intercepts 
mod.J.rs <- glmmTMB(focal.com.size ~ Immigration * Cohort + 
                      (1 + Immigration | Plant_ID) + 
                      (1 + Immigration | Lineage),
                    family = poisson(), data = my.data) 

# Random slopes and intercept for Lineage, random intercept for Plant 
mod.J.rsL <- glmmTMB(focal.com.size ~ Immigration * Cohort + 
                       (1 | Plant_ID) + 
                       (1 + Immigration | Lineage), 
                     family = poisson(), data = my.data) 

# Random slopes and intercept for Plant, random intercept for Lineage 
mod.J.rsP <- glmmTMB(focal.com.size ~ Immigration * Cohort + 
                       (1 + Immigration | Plant_ID) + 
                       (1 | Lineage), 
                     family = poisson(), data = my.data) 

# Random intercepts
mod.J.ri <- glmmTMB(focal.com.size ~ Immigration * Cohort + 
                      (1 | Plant_ID) + 
                      (1 | Lineage), 
                    family = poisson(), data = my.data) 

anova(mod.J.rs, mod.J.rsL, mod.J.rsP, mod.J.ri)
# Random slopes model fits best (lowest AIC & BIC)

# Check diagnostics: residual plots created through simulation-based approach
# in DHARMa package
simulationOutput <- simulateResiduals(mod.J.rs, plot = F, use.u = T, n = 1000)
(diag.beta <- plot(simulationOutput))
# check for potential outlier influence. Test run with and without outlier are consistent

# Get immigration effect slope for each cohort
emt <- emtrends(mod.J.rs, pairwise ~ Cohort, var = "Immigration")
summary(emt, infer = c(TRUE, TRUE))

# Model results table
tab_model(mod.J.rs, show.se = TRUE, show.std = TRUE)
#df.residual(mod.J.rs) # change df to 68 in output table

# Extract predictions (estimated marginal means) and confidence intervals for parameters
dat <- ggpredict(mod.J.rs, terms = c("Immigration", "Cohort"))
moderator_values <- sort(as.character(unique(dat$group)))

# Change facet names
cohort_names <- c('3rd_instar' = "3rd instar", '5th_instar' = "5th instar")

# for plotting only, remove outlier to better visualize predicted function. Move
# plot with outlier to suppliments.
raw.nooutlier <-  attr(dat, "rawdata") %>%
  filter(response < 300)

# Plot community size response
(j_plot <- ggplot(as.data.frame(dat), 
                  aes(x = x, y = predicted, group = group, color = as.character(group)), 
                  fill = as.character(group)) +
    geom_point(data = raw.nooutlier, aes(x = jitter(x), y = response, alpha = 0.6)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), color = FALSE, alpha = 0.2) +
    geom_line() +
    facet_wrap(~ group, labeller = as_labeller(cohort_names)) +
    scale_color_discrete(guide = "legend", aesthetics = c("colour", "fill"), 
                         breaks = moderator_values, 
                         labels = paste(format(moderator_values), 
                                        c("3rd instar", "5th instar"))) +
    scale_color_manual(values = c("darkorange", "blue")) +
    scale_fill_manual(values = alpha(c("darkorange", "blue"), .2))+
    ylab(expression(paste("Community size (copy #/", mu, "L)"))) +
    scale_x_continuous(name = expression(paste("Immigration density (1.7 CFU/", mu, "L)")),
                       breaks = c(0, 1, 2, 3, 4, 5),
                       labels = c('0' = '10^0', '1' = '10^1', 
                                  '2' = '10^2', '3' = '10^3', 
                                  '4' = '10^4', '5' = '10^5'))  +
    theme(legend.position = "none") 
)


#### 4. Model beta-diversity (Bray-Curtis) ####

# remove caterpillars with no establishment 
beta.data <- my.data %>% 
  drop_na(focal.disp.bc.rel)

# Fit zero-inflated beta regression (BEZI) model 
# NOTE: zero inflated beta-regression with random effects caused issues with 
# large z stats and negative values for random variation, causing incorrectly 
# high R2. Because random effects (lineage and plant) captured zero variation, 
# they were dropped and the model was refit using only fixed effects. NOTE, 
# summary() on gamlss object gives change in log-odds of the outcome per unit 
# change of the predictor, whereas tab_model() gives change in odds per unit
# outcomes by exponentiation the coefficients. The latter is generally more 
# interpretable.
mod.BC <- gamlss(focal.disp.bc.rel ~ std_Immigration * Cohort, 
                 family = BEZI, 
                 data = beta.data)


# Model results table
plot_model(mod.BC, type = "int", show.data = TRUE, jitter = 0.2,
           terms = c("srd_Immigration [all]", "Cohort [all]"))
tab_model(mod.BC, show.se = TRUE)
Rsq(mod.BC) # 0.71
summary(mod.BC)

# Extract predictions (estimated marginal means) and confidence intervals for parameters
cohort_names <- c('3rd_instar' = "3rd instar", '5th_instar' = "5th instar")

dat <- ggpredict(mod.BC, terms = c("std_Immigration [all]", "Cohort [all]"))
moderator_values <- sort(as.character(unique(dat$group)))

# Plot beta BR response
(fam_bc_plot <- ggplot(as.data.frame(dat), 
                       aes(x = x, y = predicted, group = group, color = as.character(group)), 
                       fill = as.character(group)) +
    geom_point(data = attr(dat, "rawdata"), aes(x = jitter(x), y = response, alpha = 0.6)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), color = FALSE, alpha = 0.2) +
    geom_line() +
    facet_wrap(~ group, labeller = as_labeller(cohort_names)) +
    scale_color_manual(values = c("darkorange", "blue")) +
    scale_fill_manual(values = alpha(c("darkorange", "blue"), .2))+
    ylab(expression(paste(beta, "-diversity"))) +
    scale_x_continuous(name = expression(paste("Immigration density (1.7 CFU/", mu, "L)")),
                       breaks = c(-1.3732664, -0.8004601, -0.2276538, 
                                  0.3451525, 0.9179588, 1.4907652),
                       labels = c('-1.3732664'= '10^0', '-0.8004601'= '10^1', 
                                "-0.2276538" = "10^2", '0.3451525'= '10^3', 
                                '0.9179588'= '10^4', "1.4907652" = "10^5")) +
    theme(legend.position = "none") 
)


#### 5. Plot model results ####
jpeg("Fig3_ModelResponses.jpg", width = 11, height = 17, units = "cm", res = 1000)
fig <- ggarrange(j_plot + rremove("xlab"), fam_rich_plot + rremove("xlab"), fam_bc_plot + rremove("xlab"), 
                 ncol = 1, nrow = 3, align = "v",
                 labels = c("A", "B", "C"))
annotate_figure(fig, 
                bottom = textGrob(expression(paste("Immigration density (1.7 CFU/", mu, "L)")),
                                  gp = gpar(fontsize = 10)))
dev.off()

