library(tidyverse)
library(permuter)

#Replicating analysis of data in Tautges 2019: "Deep soil inventories reveal that impacts of cover crops and compost on soil carbon sequestration differ in surface and subsurface soils"

#read in data
data <- read_csv("../Data/Tautges_LTAR_data.csv") %>%
  filter(depthcode != 6) %>% #this is originally where whole profile data was stored
  mutate(Plot = gsub("_", "-", Plot)) %>%
  filter(Treatment != "CWT") #CWT was removed and not reported

#compute whole profile stock and concentration data
#computation of whole profile results seems off, especially for concentration (see e.g. figure 2...)
#it seems concentration was computed for whole profile based on a straight mean (not a mean weighted by length)
#whole profile stock is not a straight sum of the other quantities?
wp_data <- data %>%
  group_by(year, Plot, Treatment, block) %>%
  summarize(percC = mean(percC), TotC = sum(TotC_Mgha)) %>%
  #summarize(percC = sum(Length / sum(Length) * percC), TotC = sum(TotC_Mgha)) %>%
  mutate(depthcode = 6)

############ initial checks #############

#3 blocks
unique(data$block)

#nine cropping systems (treatments)
unique(data$Treatment)

#from Jessica Chiartas:
#OMT = ORG; CMT = CONV; LMT = CONV+WCC; IWC = IWF; IWF = IWF+N; RWC = RWF; RWF = RWF+N; RWL = RWF+WCC.
#CWT got dropped in the paper
#fertilizer=1 => poulty manure; fertilizer=3 => 168 kg/ha synthetic nitrogen; fertilizer=0 => none; fertilizer=2 => winter cover crop (WCC)
#Irrigation=2 => furrow; Irrigation=0 => None; Irrigation=1 => sprinkler

#each treatment is replicated 6 times (twice per block)
replicates_per_block <- data %>% 
  group_by(block, Treatment) %>% 
  summarize(reps_per_block = n_distinct(Plot))

#table 2: initial SOC concentrations and bulk densities (no pH)
table_2 <- data %>%
  filter(year == 1993) %>%
  group_by(Lower_Depth) %>%
  summarize(C_concentration = mean(percC)*10, C_to_N = mean(percC / percN), BD = mean(`Bulk Density`), clay = mean(clay))

#marginal distributions of concentration and stock
hist(data$percC)
hist(data$TotC_Mgha)


#normality or log-normality within treatments?
change_data <- data %>%
  select(year, depthcode, block, Plot, Treatment, TotC = TotC_Mgha, percC) %>%
  bind_rows(wp_data) %>%
  pivot_wider(values_from = c(TotC, percC), names_from = year) %>%
  mutate(TotC_diff = TotC_2012 - TotC_1993, percC_diff = percC_2012 - percC_1993, block = as_factor(block))

#note that there really aren't enough replicates to empirically assess normality (within blocks, treatments) at all
normality_plot <- ggplot(change_data, aes(percC_diff)) +
  geom_histogram() +
  facet_grid(depthcode ~ .)


#figure 2 SOC changes in maize
change_data_corn <- change_data %>%
  filter(Treatment %in% c("OMT","CMT","LMT")) %>%
  group_by(Treatment, depthcode) %>%
  summarize(mean_concentration = mean(percC_diff), se_concentration = sqrt(var(percC_diff)/n()), mean_total = mean(TotC_diff), se_total = sqrt(var(TotC_diff)/n()), sample_size = n()) %>%
  pivot_longer(c("mean_concentration", "mean_total", "se_concentration", "se_total"), values_to = "value") %>%
  separate(name, into = c("statistic", "type"), sep = "_") %>%
  pivot_wider(names_from = "statistic", values_from = "value")


figure_2 <- ggplot(change_data_corn, aes(x = depthcode, y = mean, group = Treatment, fill = Treatment)) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), position = position_dodge(1), width = 0.2) +
  scale_x_reverse() +
  coord_flip() +
  facet_grid(~ type, scales = "free")


#figure 3 SOC changes in wheat
change_data_wheat <- change_data %>%
  filter(Treatment %in% c("OMT","CMT","LMT")) %>%
  group_by(Treatment, depthcode) %>%
  summarize(mean_concentration = mean(percC_diff), se_concentration = sqrt(var(percC_diff)/n()), mean_total = mean(TotC_diff), se_total = sqrt(var(TotC_diff)/n()), sample_size = n()) %>%
  pivot_longer(c("mean_concentration", "mean_total", "se_concentration", "se_total"), values_to = "value") %>%
  separate(name, into = c("statistic", "type"), sep = "_") %>%
  pivot_wider(names_from = "statistic", values_from = "value")


figure_3 <- ggplot(change_data_wheat, aes(x = depthcode, y = mean, group = Treatment, fill = Treatment)) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), position = position_dodge(1), width = 0.2) +
  scale_x_reverse() +
  coord_flip() +
  facet_grid(~ type, scales = "free")



#plot 1993 versus 2012 percent SOC
carbon_plot <- ggplot(change_data, aes(x = percC_1993, y = percC_2012)) +
  geom_point() +
  facet_grid(depthcode ~ .)


############# paired t tests for change within treatments, depths ##########
#I suspect they used a paired t.test (see esp results for OMT aka ORG and LMT aka CONV+WCC)
#I can't get the exact right p-value for OMT, depthcode = 2 (15-30cm), which should be .006
paired_t_test_concentration_results <- change_data %>%
  group_by(depthcode, Treatment) %>%
  summarize(point_estimate = mean(percC_2012 - percC_1993), test_stat = sqrt(n()) * mean(percC_2012 - percC_1993) / sd(percC_2012 - percC_1993), size = n()) %>%
  mutate(parametric_p_value = 2*(1 - pt(abs(test_stat), df = size - 1)))

# 2.59 mg C / ha reported for 15-30cm ORG at top of page 7 is a typo? I estimate 2.42 Mg C / ha, p-value of .04 (not .01)
paired_t_test_stock_results <- change_data %>%
  group_by(depthcode, Treatment) %>%
  summarize(point_estimate = mean(TotC_2012 - TotC_1993), test_stat = sqrt(n()) * mean(TotC_2012 - TotC_1993) / sd(TotC_2012 - TotC_1993), size = n()) %>%
  mutate(parametric_p_value = 2*(1 - pt(abs(test_stat), df = size - 1)))


############ parametric ANOVA ##############
#parametric ANOVA was actually not used in the original paper; they used a "visual inference" method (examining the barcharts with error bars)
#compare deltas (differences) across treatments with a blocking effect (two-way ANOVA)

#topsoil models
concentration_topsoil_model <- lm(percC_diff ~ Treatment + block, data = change_data %>% filter(depthcode == 1))
stock_topsoil_model <- lm(TotC_diff ~ Treatment + block, data = change_data %>% filter(depthcode == 1))
#results
anova(concentration_topsoil_model)
summary(concentration_topsoil_model)
anova(stock_topsoil_model)
summary(stock_topsoil_model)

############# permutation tests ################
#number of permutations
B <- 10000
#using permuter, check for paired differences
pt_stock_change <- change_data %>%
  group_by(Treatment, depthcode) %>%
  summarize(diff_in_means = mean(TotC_diff), perm_p_value = t2p(tst = mean(TotC_diff), distr = one_sample(x = TotC_diff, reps = B), alternative = c("two-sided")))

pt_concentration_change <- change_data %>%
  group_by(Treatment, depthcode) %>%
  summarize(diff_in_means = mean(percC_diff), perm_p_value = t2p(tst = mean(percC_diff), distr = one_sample(x = percC_diff, reps = B), alternative = c("two-sided")))


#compare permutation p-values to parameteric p-values
joined_concentration_change_pvalues <- pt_concentration_change %>%
  left_join(paired_t_test_concentration_results, by = c("Treatment", "depthcode")) %>%
  select(Treatment, depthcode, diff_in_means, perm_p_value, parametric_p_value)

joined_stock_change_pvalues <- pt_stock_change %>%
  left_join(paired_t_test_stock_results, by = c("Treatment", "depthcode")) %>%
  select(Treatment, depthcode, diff_in_means, perm_p_value, parametric_p_value)

#scatter plot of permutation pvalues against parametric pvalues for whole profile
concentration_change_pvalue_plot <- ggplot(joined_concentration_change_pvalues %>% filter(depthcode == 6), aes(x = parametric_p_value, y = perm_p_value, label = Treatment)) + 
  geom_text(size = 3) +
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") 

stock_change_pvalue_plot <- ggplot(joined_stock_change_pvalues %>% filter(depthcode == 6), aes(x = parametric_p_value, y = perm_p_value, label = Treatment)) + 
  geom_text(size = 3) +
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") 


# FDR control by Benjamini-Hochberg
# Are the tests independent? Probably not. If a treatment changes one plot or depth it is likely to effect another similarly. The B-H procedure may still be valid. 
#the Benjamini-Yekutieli procedure is valid under arbitrary dependence
#assume FDR control at 0.1 level
level <- 0.1
m <- nrow(pt_concentration_change)
pt_concentration_change <- pt_concentration_change %>%
  ungroup() %>%
  arrange(perm_p_value) %>%
  mutate(p_value_order = rank(perm_p_value)) %>%
  mutate(bh_threshold = (1:m * level) / m) %>%
  mutate(by_threshold = (1:m * level) / (m * sum(1/1:m)))
  
#no test rejects under BH or BY
FDR_plot <- ggplot(pt_concentration_change, aes(y = perm_p_value, x = p_value_order)) +
  geom_point() +
  geom_hline(yintercept = level, linetype = "dashed", color = "red") +
  geom_hline(yintercept = level / m, linetype = "dashed", color = "blue") +
  geom_line(aes(x = p_value_order, y = bh_threshold), linetype = "dashed", color = "forestgreen")


#NPC for testing intersection null (treatment has an effect on some outcome)




