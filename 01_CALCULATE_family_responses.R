#----
# R script to calculate microbe abundance, diversity, and null-model deviation
# for Immigration Experiment
# Creator: Christopher P. Catano, chcatano@gmail.com
# Data input: "all_family_data.csv"
# Data output: "family_diveristy_data_final.csv"
#----


# Set up environment 
rm(list = ls())

# Check for and install/load required packages
for (package in c('tidyverse', 'vegan', 'ggpubr', 'NST')) {
  if (!require(package, character.only = T, quietly = T)) {
    install.packages(package)
    library(package, character.only = T)
  }
}


#### 1. Import & prep data ####
# meta and treatment data
fam.data <- read.csv("data/all_family_data.csv", row.names = 1)
fam.data$Cat_ID <- as.character(fam.data$Cat_ID)

# split into treatment/response data vs community matrix
#treatment data
trt.data <- fam.data %>%
  select(c(Cat_ID, Lineage, Treatment, Plant_ID, Cohort, 
           Microbe_Mix, log_community_size))


#### 2. Plot relative abundances of focal bacteria ####
# set graphics parameters
theme_set(theme_bw() +  
            theme(legend.position = "bottom", panel.grid.minor = element_blank(), 
                  panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5), 
                  text = element_text(size = 10), axis.text = element_text(size = 8)))

# Determine relative abundance of focal taxa and non focal taxa
fam.mat <- fam.data %>%
  select(-c(Cat_ID, Lineage, Treatment, Plant_ID, Cohort, Microbe_Mix, log_community_size))
fam.mat.rel <- as.data.frame(t(apply(fam.mat, 1, function(x) x/sum(x))))

# Create stacked bar plot but first group all non-focal species and get their
# summed relative abundance
fam.mat.rel.nonfocal <- fam.mat.rel %>%
  select(-c(Comamonadaceae, Bacillaceae, Enterobacteriaceae, Xanthomonadaceae, Pseudomonadaceae))

fam.mat.rel.summed.nonfocal <- as.data.frame(apply(fam.mat.rel.nonfocal, 1, function(x) sum(x)))
colnames(fam.mat.rel.summed.nonfocal) <- "non focal"

fam.mat.rel.focal <- fam.mat.rel %>%
  select(c(Comamonadaceae, Bacillaceae, Enterobacteriaceae, Xanthomonadaceae, Pseudomonadaceae))

# Bind together
RA_data_all <- cbind(trt.data, fam.mat.rel.focal, fam.mat.rel.summed.nonfocal)
RA_data_all2 <- RA_data_all %>%
  pivot_longer(cols = -c(Cat_ID, Lineage, Treatment, Plant_ID, Cohort, Microbe_Mix, log_community_size), 
               names_to = "Taxa", 
               values_to = 'r.abund')

# Calculate proportion of all microbes represented by focal taxa for each treatment
prop <- RA_data_all 
prop$rel.prop <- 1 - RA_data_all$`non focal`
grouped <- prop %>%
  group_by(Treatment, Cohort) %>%
  summarise(mean = mean(rel.prop), sd = sd(rel.prop))

# Plot stacked bar plot
treat_names <- c("C" = "Control", 'D0'= '10^5', 'D1'= '10^4', "D2" = "10^3", 
                 'D3'= '10^2', 'D4'= '10^1', "D5" = "10^0")

(taxBar_3rd <- RA_data_all2 %>%
  filter(Cohort == "3rd_instar") %>%
  ggplot(aes(fill = factor(Taxa, levels = c("non focal", "Comamonadaceae", "Pseudomonadaceae", "Xanthomonadaceae",
                                            "Bacillaceae", "Enterobacteriaceae")), 
                           y = r.abund, x = Cat_ID)) + 
  geom_bar(position = "stack", stat = "identity") +
  facet_grid(~ factor(Treatment, levels = c("C", 'D5', 'D4', "D3", 'D2', 'D1', "D0")), 
             scales = "free_x", labeller = as_labeller(treat_names)) +
  scale_fill_manual(values = c("#999999", "#fde725", "#3b528b", "#5ec962", "#440154", "#21918c")) +
  ggtitle(expression(paste('Development stage: ', 3^{'rd'}, ' Instar'))) +
  ylab("Relative abundance") +
  xlab("Individual caterpillars") +
  labs(fill = "Family") +
  theme(axis.text.x = element_blank())
)

(taxBar_5th <- RA_data_all2 %>%
  filter(Cohort == "5th_instar") %>%
  ggplot(aes(fill = factor(Taxa, levels = c("non focal", "Comamonadaceae", "Pseudomonadaceae", "Xanthomonadaceae",
                                            "Bacillaceae", "Enterobacteriaceae")), 
             y = r.abund, x = Cat_ID)) + 
  geom_bar(position = "stack", stat = "identity") +
  facet_grid(~ factor(Treatment, levels = c("C", 'D5', 'D4', "D3", 'D2', 'D1', "D0")), 
             scales = "free_x", labeller = as_labeller(treat_names)) +
  scale_fill_manual(values = c("#999999","#fde725", "#3b528b", "#5ec962", "#440154", "#21918c")) +
  ggtitle(expression(paste('Development stage: ', 5^{'th'}, ' Instar'))) +
  ylab("Relative abundance") +
  xlab("Individual caterpillars") +
  labs(fill = "Family") +
  theme(axis.text.x = element_blank())
)


#jpeg("figures/Taxa_Barplot_fig.jpg", width = 22, height = 16, units = "cm", res = 600)
ggarrange(taxBar_3rd, taxBar_5th, 
          common.legend = TRUE, legend="bottom", 
          ncol = 1, nrow = 2, align = "v")
#dev.off()


#### 3. Plot abundance & relative abund. of focal taxa ####
# Note, relative abundance here is each taxa relative to only the other focal taxa,
# note all bacterial taxa in the gut
fam.mat.focal2 <- fam.data %>%
  select(c(Comamonadaceae, Bacillaceae, Enterobacteriaceae, Xanthomonadaceae, Pseudomonadaceae))
fam.mat.rel2 <- as.data.frame(t(apply(fam.mat.focal2, 1, function(x) x/sum(x))))

# Bind together
RA_data_focal <- cbind(trt.data, fam.mat.rel2)
RA_data_focal_long <- RA_data_focal %>%
  pivot_longer(cols = -c(Cat_ID, Lineage, Treatment, Plant_ID, Cohort, Microbe_Mix, log_community_size), 
               names_to = "Taxa", 
               values_to = 'r.abund')

# Calculate mean relative abundance of Enterobacteriaceae across all treatments in 3rd vs 5th cohorts
RA_data_focal_long %>%
  filter(Taxa == "Enterobacteriaceae") %>%
  filter(Treatment != "C") %>%
  select(Cohort, Taxa, r.abund) %>%
  group_by(Cohort) %>%
  summarize(Mean = mean(r.abund, na.rm = TRUE),
            SD   =   sd(r.abund, na.rm = TRUE)) %>%
  
  ungroup() 

# A tibble: 2 × 3
#  Cohort      Mean    SD
#  <chr>      <dbl> <dbl>
#1 3rd_instar 0.486 0.366
#2 5th_instar 0.816 0.337

# If community size is NA (qPCR not run for a few samples) then assign mean
# value calculated for each cohort by immigration treatment
data_totalAbund2 <- RA_data_focal_long %>%
  group_by(Treatment, Cohort) %>%
  mutate(mean_log_community_size = mean(log_community_size, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(log_community_size = coalesce(log_community_size, mean_log_community_size))

# Multiple relative abundance by qPCR result to estimate abundance of each taxa,
# rounding up to the nearest integer
data_totalAbund2$abund <- ceiling(data_totalAbund2$r.abund * exp(data_totalAbund2$log_community_size))

# Remove control treatment and jitter points so that trends by taxa can be visualized
abund.change.data <- data_totalAbund2 %>% 
  select(-c(mean_log_community_size)) %>%
  filter(Treatment != "C") %>%
  mutate(
    x_pos = +(Cohort == '5th_instar'),
    x_pos = if_else(Taxa == 'Comamonadaceae', x_pos - .08, x_pos + .08),
    x_pos = if_else(Taxa == 'Xanthomonadaceae', x_pos - .02, x_pos + .02),
    x_pos = if_else(Taxa == 'Bacillaceae', x_pos - .04, x_pos + .04),
    x_pos = if_else(Taxa == 'Enterobacteriaceae', x_pos - .06, x_pos + .06),
  )

# order treatments for plotting
treat2_names <- c('D0'= '10^5', 'D1'= '10^4', "D2" = "10^3", 
                  'D3'= '10^2', 'D4'= '10^1', "D5" = "10^0")

# plot change in relative abundance (differences between 3rd and 5th instars)
# Note, caterpillars where no focal taxa established are not included in means/CIs
(plot.rabund.change <- abund.change.data %>%
  ggplot(aes(x = x_pos, y = r.abund, 
             colour = factor(Taxa, levels = c("Comamonadaceae", "Pseudomonadaceae", 
                                              "Xanthomonadaceae", "Bacillaceae", 
                                              "Enterobacteriaceae")))) +
    stat_summary(fun = mean, geom = "line", aes(group = Taxa)) +
    stat_summary(fun = mean, geom = "point", aes(group = Taxa), size = 2) +
    stat_summary(fun.data = mean_cl_normal, geom = "linerange", fun.args = list(mult = 1)) +
    facet_wrap(~ factor(Treatment, levels = c('D5', 'D4', "D3", 'D2', 'D1', "D0")), 
                labeller = as_labeller(treat2_names)) +
    ylab("Relative abundance") +
    scale_x_continuous(breaks = c(0.15, 1.1), labels = c("3rd", "5th")) + 
    xlab("Time (instar)") +
    scale_color_manual(values = c("#fde725", "#3b528b", "#5ec962", "#440154", "#21918c")) +
    labs(color = "Family")
)

#jpeg("figures/RA_change_fig.jpg", width = 7.5, height = 6, units = "in", res = 600)
plot.rabund.change
#dev.off()

# plot change in abundance (differences between 3rd and 5th instars)
# Note, caterpillars where no focal taxa established are not included in means/CIs
(plot.abund.change <- abund.change.data %>%
    ggplot(aes(x = x_pos, y = log(abund + 1), 
               colour = factor(Taxa, levels = c("Comamonadaceae", "Pseudomonadaceae", 
                                                "Xanthomonadaceae", "Bacillaceae", 
                                                "Enterobacteriaceae")))) +
    stat_summary(fun = mean, geom = "line", aes(group = Taxa)) +
    stat_summary(fun = mean, geom = "point", aes(group = Taxa), size = 2) +
    stat_summary(fun.data = mean_cl_normal, geom = "linerange", fun.args = list(mult = 1)) +
    facet_wrap(~ factor(Treatment, levels = c('D5', 'D4', "D3", 'D2', 'D1', "D0")), 
               labeller = as_labeller(treat2_names)) +
    ggtitle("Bacterial abundance") +
    ylab(expression("Log 16S rRNA copies per "*mu*"L")) +
    scale_x_continuous(breaks = c(0.15, 1.1), labels = c("3rd", "5th")) + 
    xlab("Development stage (instar)") +
    scale_color_manual(values = c("#fde725", "#3b528b", "#5ec962", "#440154", "#21918c")) +
    labs(color = "Family")
)

jpeg("Fig5.jpg", width = 7.5, height = 6.5, units = "in", res = 600)
plot.abund.change
dev.off()


#### 4. Calculate alpha and beta diversity ####
# Convert abundance data from long to wide format
data_abund_wide <- data_totalAbund2 %>%
  select(-c(r.abund, mean_log_community_size)) %>%
  pivot_wider(names_from = Taxa, values_from = abund) %>%
  replace_na(list(Comamonadaceae = 0, Bacillaceae = 0, Enterobacteriaceae = 0, 
                  Xanthomonadaceae = 0, Pseudomonadaceae = 0)) 

data_relAbund_wide <- data_totalAbund2 %>%
  select(-c(abund, mean_log_community_size)) %>%
  pivot_wider(names_from = Taxa, values_from = r.abund) %>%
  replace_na(list(Comamonadaceae = 0, Bacillaceae = 0, Enterobacteriaceae = 0, 
                  Xanthomonadaceae = 0, Pseudomonadaceae = 0)) 

# Focal taxa richness (q = 0)
trt.data$focal.alpha.0 <- data_abund_wide %>%
  select(c("Comamonadaceae", "Bacillaceae", "Enterobacteriaceae", 
           "Xanthomonadaceae", "Pseudomonadaceae")) %>%
  specnumber()

trt.data$focal.alpha.2 <- data_abund_wide %>%
  select(c("Comamonadaceae", "Bacillaceae", "Enterobacteriaceae", 
           "Xanthomonadaceae", "Pseudomonadaceae")) %>%
  diversity(index = "invsimpson", MARGIN = 1, base = exp(1))
trt.data[sapply(trt.data, is.infinite)] <- 0

trt.data$focal.com.size <- data_abund_wide %>%
  select(c("Comamonadaceae", "Bacillaceae", "Enterobacteriaceae", 
           "Xanthomonadaceae", "Pseudomonadaceae")) %>%
  rowSums()

# Plot observed richness
(plot.focal.rich <- ggplot(trt.data, aes(x = Treatment, y = focal.alpha.0)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = NA)+
    facet_grid(~ Cohort) +
    geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1, alpha = 0.3,
                 position = position_dodge(1)) +
    ylab("Richness") +
    scale_x_discrete(name = expression(paste("Immigration density (1.7 CFU/", mu, "L)")),
                     labels = c("C" = "C", 'D0'= '10^5', 'D1'= '10^4', "D2" = "10^3", 
                              'D3'= '10^2', 'D4'= '10^1', "D5" = "10^0")) #+
)

# Plot focal diversity (Inverse Simpsons)
(plot.focal.div <- ggplot(trt.data, aes(x = Treatment, y = focal.alpha.2)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = NA)+
    facet_grid(~ Cohort) +
    geom_dotplot(binaxis = 'y', stackdir = 'center', alpha = 0.3,
                 position = position_dodge(1)) +
    ylab("Diversity (ENSpie)") +
    scale_x_discrete(name = expression(paste("Immigration density (1.7 CFU/", mu, "L)")),
                     labels = c("C" = "C", 'D0'= '10^5', 'D1'= '10^4', "D2" = "10^3", 
                              'D3'= '10^2', 'D4'= '10^1', "D5" = "10^0")) #+
)

# plot focal community size
(plot.focal.j <- ggplot(trt.data, aes(x = Treatment, y = focal.com.size)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = NA)+
    facet_grid(~ Cohort) +
    geom_dotplot(binaxis = 'y', stackdir = 'center', alpha = 0.3,
                 position = position_dodge(1)) +
    ylab("Community size") +
    scale_x_discrete(name = expression(paste("Immigration density (1.7 CFU/", mu, "L)")),
                     labels = c("C" = "C", 'D0'= '10^5', 'D1'= '10^4', "D2" = "10^3", 
                              'D3'= '10^2', 'D4'= '10^1', "D5" = "10^0")) 
)

# Observed alpha diversity and community size data
jpeg("figures/observed_localDiversity_data.jpg", width = 5, height = 8, units = "in", res = 600)
ggarrange(plot.focal.j, plot.focal.rich, plot.focal.div,  ncol = 1, nrow = 3,
          widths = c(1, 1),  align = "hv")
dev.off()


#### 5. Calculate beta diversity for focal bacteria ####
# Make relative counts for Bray-Curtis since it is sensitive to total counts
fam.mat.focal.rel <- data_relAbund_wide %>%
  select(c("Comamonadaceae", "Bacillaceae", "Enterobacteriaceae", 
           "Xanthomonadaceae", "Pseudomonadaceae"))

# Merge with treatment data and remove individuals where no focal taxa established 
full.focal.data <- cbind(trt.data, fam.mat.focal.rel)
full.focal.data <- full.focal.data %>%
  select(-c(log_community_size)) %>%
  unite("trt_coh", Treatment, Cohort, remove = FALSE) %>%
  filter(rowSums(.[11:15]) > 0)

# This gives us our relativized community matrix where caterpillars with no focal
# taxa established are removed.
focal.rel.mat <- full.focal.data %>%
  select(c(Comamonadaceae, Bacillaceae, Enterobacteriaceae, Xanthomonadaceae, Pseudomonadaceae))

# Calculate Bray-Curtis dissimilarity
focal.bc <- vegdist(focal.rel.mat, method = "bray")
# Calculate the multivariate dispersion of BC
focal.disp.bc <- betadisper(focal.bc, full.focal.data$trt_coh)
boxplot(focal.disp.bc)
# Save distance-to-centroid for each each caterpillar 
full.focal.data$focal.disp.bc.rel <- focal.disp.bc$distances 

# Plot observed beta diversity
(plot.focal.betaBC <- ggplot(full.focal.data, aes(x = Treatment, y = focal.disp.bc.rel)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = NA)+
    facet_grid(~Cohort) +
    geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1, alpha = 0.3,
                 position = position_dodge(1)) +
    ylab("Distance-to-centroid\n(Bray-Curtis dissimilarity)") +
    scale_x_discrete(name=expression(paste("Immigration density (1.7 CFU/", mu, "L)")),
                     labels=c("C" = "C", 'D0'= '10^5', 'D1'= '10^4', "D2" = "10^3", 
                              'D3'= '10^2', 'D4'= '10^1', "D5" = "10^0")) 
)


#### 6. Finalize and export data set for analysis ####
export_ra_data <- trt.data %>%
  select(-c(log_community_size, focal.com.size)) %>%
  left_join(full.focal.data, 
            by = c("Cat_ID", "Lineage", "Treatment", "Plant_ID", "Cohort",
                   "Microbe_Mix", "focal.alpha.0", "focal.alpha.2")) %>%
  replace_na(list(Comamonadaceae = 0, Bacillaceae = 0, Enterobacteriaceae = 0, 
                  Xanthomonadaceae = 0, Pseudomonadaceae = 0, focal.com.size = 0)) %>%
  select(-c(trt_coh, Microbe_Mix))

write.csv(export_ra_data, "data/family_diveristy_data_final.csv")


#### 7. Null model analysis ####

# with abundance data
export_abund_data <- trt.data %>%
  select(-c(log_community_size)) %>%
  left_join(data_abund_wide, 
            by = c("Cat_ID", "Lineage", "Treatment", "Plant_ID", "Cohort",
                   "Microbe_Mix")) %>%
  replace_na(list(Comamonadaceae = 0, Bacillaceae = 0, Enterobacteriaceae = 0, 
                  Xanthomonadaceae = 0, Pseudomonadaceae = 0, focal.com.size = 0)) %>%
  select(-c(log_community_size))

full.data <- export_abund_data %>%
  unite("trt_coh", Treatment, Cohort, remove = FALSE) %>%
  filter(rowSums(.[10:14]) > 0)

# this gives us our community matrix where caterpillars where no focal
# taxa established are removed.
full.mat <- full.data %>%
  select(c(Comamonadaceae, Bacillaceae, Enterobacteriaceae, Xanthomonadaceae, Pseudomonadaceae))


# Null models and stochasticity quantification ####
# Calculate SES (Kraft et al. 2011) as measure of relative stochasticity 
comm <- full.mat
group <- as.data.frame(full.data$trt_coh)
colnames(group) <- "trt_coh"

groupi = group[, 1, drop = FALSE]
# create pool where each species has the same probability of being sampled
meta.com = matrix(1, nrow = 1, ncol = ncol(comm))
colnames(meta.com) <- colnames((comm))

# Equiprobable occurrence frequency, Proportional richness in each sample
tnst.EP <- tNST(comm = comm, group = groupi, meta.group = NULL, meta.com = meta.com,
                dist.method = "bray", abundance.weighted = TRUE, rand = 1000,
                output.rand = FALSE, nworker = 1, LB = FALSE, null.model = "EP",
                between.group = FALSE, SES = TRUE, RC = FALSE)
tnst.EP.grp <- tnst.EP$index.grp
tnst.EP.grp$null.model <- "EP"
tnst.EP.pair.grp <- tnst.EP$index.pair.grp
tnst.EP.pair.grp$null.model <- "EP"

long.summaries <- tnst.EP.grp %>%
  select(c(group, size, SES.i.bray)) %>%
  separate_wider_delim(group, "_", names = c("Treatment", "Instar", "drop"))

# Plot null model deviations
instar_names <- c('3rd'= '3rd Instar', '5th'= '5th Instar')

(ses.plot <- long.summaries %>%
  filter(Treatment != "C") %>%
  ggplot(aes(x = factor(Treatment, levels = c("D5", "D4", "D3", "D2", "D1", "D0")), 
                        y = SES.i.bray, color = Instar, group = Instar, fill = Instar)) + 
  geom_point(size = 2) +
  geom_line() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -2, linetype = "solid", color = "black") +
  ylab("SES") +
  facet_wrap(~ factor(Instar, levels = c('3rd', '5th')), labeller = as_labeller(instar_names)) +
  scale_color_manual(values = c("darkorange", "blue")) +
  scale_fill_manual(values = alpha(c("darkorange", "blue"), .2))+
  scale_x_discrete(name = expression(paste("Immigration density (1.7 CFU/", mu, "L)")),
                   labels = c( "D5" = "10^0", 'D4'= '10^1', 'D3'= '10^2', 
                             "D2" = "10^3",'D1'= '10^4', 'D0'= '10^5')) +
  theme(legend.position = "none")
  
)

jpeg("figures/beta_SES_working.jpg", width = 7, height = 4, units = "in", res = 600)
ses.plot
dev.off()
   
