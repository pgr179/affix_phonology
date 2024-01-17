# Packages
library(tidyverse)
library(lme4)
library(lmerTest)
library(DescTools)
library(effects)
library(sjPlot)
library(ggeffects)
library(car)
library(data.table)
library(performance)


# Read in the data
data <- read.csv(file.choose(), header=TRUE)


table(data$Type)
# Fill in language names
data$Language[data$Language == ""] = NA
data = data %>% fill(Language)

# Limit to just prefixes and suffixes
data = data[data$Type == "P" | data$Type == "S",]

# Set data types
data$Language = as.factor(data$Language)
data$Type = as.factor(data$Type)
data$Syllables = as.numeric(data$Syllables)
data$Segments = as.numeric(data$Segments)
data$Allomorphs = as.numeric(data$Allomorphs)

# Create binary Allomorphy variable
data$Allomorphy = as.numeric(data$Allomorphs >= 2)


## Exploration of the data in tables and plots


# Number of affixes more than two syllables and more than two segments
nrow(data)

nrow(data[data$Segments > 2,])

nrow(data[data$Syllables > 2,])

# Segment lengths
ggplot(data, aes(x = Segments)) +
  geom_histogram()

# Syllables
ggplot(data, aes(x = Syllables)) +
  geom_histogram()

# Allomorphs
ggplot(data, aes(x = Allomorphs)) +
  geom_histogram()


# Number of prefixes and suffixes for each language
rename(count(data, Language, Type), Freq = n)


# Affix lengths by type and language (using averages and standard deviations)
segment_table = data %>% group_by(Language, Type) %>% 
  summarise(mean_segments = mean(Segments), sd_segments = sd(Segments), .groups = 'drop') %>%
  as.data.frame()
segment_table$sd_max = segment_table$mean_segments + segment_table$sd_segments
segment_table$sd_min = segment_table$mean_segments - segment_table$sd_segments

ggplot(segment_table, aes(x = Language)) +
  geom_pointrange(aes(y = mean_segments, ymin = sd_min, ymax = sd_max, color = Type), position = position_dodge(.5)) +
  labs(y = "Segments (means and standard deviations)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))


# Lengths of prefixes vs. suffixes
ggplot(data, aes(x = Segments, group = Type, color = Type)) +
  geom_histogram(aes(y = ..count..), fill = "darkgrey", bins = 33,
               data = ~ subset(., Type %in% c("P"))) +
  geom_histogram(aes(y = -..count..), fill = "darkgrey", bins = 33,
               data = ~ subset(., Type %in% c("S"))) +
  labs(y = "Count", x = "Number of segments")


# Lengths of prefixes vs. suffixes
ggplot(data, aes(x = Syllables, group = Type, color = Type)) +
  geom_histogram(aes(y = ..count..), fill = "darkgrey", bins = 33,
                 data = ~ subset(., Type %in% c("P"))) +
  geom_histogram(aes(y = -..count..), fill = "darkgrey", bins = 33,
                 data = ~ subset(., Type %in% c("S"))) +
  labs(y = "Count", x = "Number of syllables")


# Relationship of allomorphy to affix type
table(data$Type, data$Allomorphy)
chisq.test(data$Type, data$Allomorphy, correct = FALSE)
# X-squared = 0.41203, df = 1, p-value = 0.5209

ggplot(data, aes(x = Allomorphs, group = Type, color = Type)) +
  geom_histogram(aes(y = ..count..), fill = "darkgrey", bins = 33,
                 data = ~ subset(., Type %in% c("P"))) +
  geom_histogram(aes(y = -..count..), fill = "darkgrey", bins = 33,
                 data = ~ subset(., Type %in% c("S"))) +
  labs(y = "Count", x = "Number of allomorphs")


# Plotting segment length as a function of the presence of allomorphy
ggplot(data, aes(y = Allomorphs, x = Segments)) +
  geom_jitter() +
  geom_smooth(method = lm)

ggplot(data, aes(group = Allomorphs, y = Segments)) +
  geom_boxplot()




## Hypothesis 1

# Independent variables: Monosyllabicity (monosyllabic vs. non-monosyllabic)
# Dependent variables: Count
# Random effects: Languages


# Create binary variable for monosyllabicity
data$Syllabicity = ifelse(data$Syllables == 1, "Monosyllabic", "Non-monosyllabic")

# Quick test to see if monosyllabicity varies based on affix type
chisq.test(table(data$Syllabicity, data$Type))

# Create new dataframe with count data
count_data = reshape2::melt(table(data$Syllabicity, data$Language), value.name = "Count")
colnames(count_data)[1] = "Syllabicity"
colnames(count_data)[2] = "Language"
count_data$Syllabicity = relevel(count_data$Syllabicity, "Non-monosyllabic")


# Regression model


# Fixed effects only
#m200 = glm(Count ~ Language + Syllabicity, count_data, family = poisson(link = "sqrt"))
#m201 = glm(Count ~ Language, count_data, family = poisson(link = "sqrt"))
#m201
#anova(m201, m200)
#plot(m200)


# Initial model
m100 = glmer(Count ~ Syllabicity + (1 + Syllabicity|Language), count_data, family = poisson(link = "sqrt"))
m100


# Alternative model and comparison
m101 = glmer(Count ~ 1 + (1 + Syllabicity|Language), count_data, family = poisson(link = "sqrt"))
anova(m100, m101)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# m101    4 428.75 436.40 -210.38   420.75                         
# m100    5 418.88 428.44 -204.44   408.88 11.869  1  0.0005708 ***


# Get confidence intervals for the model
ci = confint(m100)
ci

pp <- profile(m100)

ggplot(as.data.frame(pp),aes(.focal,.zeta))+
  geom_point()+geom_line()+
  facet_wrap(~.par,scale="free_x")+
  geom_hline(yintercept=0,colour="gray")+
  geom_hline(yintercept=c(-1.96,1.96),linetype=2,
             colour="gray")

# Quick look
plot(allEffects(m200))

# Check residuals
plot(m100)

ggplot(data.frame(eta=predict(m100,type="link"),pearson=residuals(m100,type="pearson")),
  aes(x=eta,y=pearson)) +
  geom_point() +
  theme_bw()

rs = DHARMa::simulateResiduals(m100)
plot(rs)
DHARMa::plotResiduals(rs, form = data$Type, quantreg = T)

# Check for overdispersion
check_overdispersion(m100)

# Plot random effects
plot_model(m100, ,type="pred", terms=c("Syllabicity", "Language"), pred.type="re", ci.lvl = NA)


# Generate predictions for plotting
preds = as.data.frame(predictorEffect("Syllabicity", m100))


# Histogram of syllable lengths per language
tiff("syllable_lengths_by_language.tiff", units="in", width=8, height=8, res=300)
ggplot(data, aes(x = Syllables)) +
  geom_histogram(aes(fill = Syllabicity), bins = 9) +
  facet_wrap(vars(Language), ncol = 5) +
  labs(y = "Number of Affixes") +
  scale_fill_manual(values = c("blue", "grey30")) +
  theme(legend.position = "none")
dev.off()


# Plot predictions of model

tiff("count_by_monosyllabicity.tiff", units="in", width=3, height=5, res=300)
ggplot() +
  geom_pointrange(data = preds, aes(x = Syllabicity, y = fit, ymin = lower, ymax = upper)) +
  geom_errorbar(data = preds, aes(x = Syllabicity, y = fit, ymin = lower, ymax = upper), width = .05) +
  #geom_rug(data = p, sides = "l", position = position_jitter(height = 0.15),
  # aes(x = Type, y = Segments), alpha = 1/20) +
  #geom_rug(data = s, sides = "r", position = position_jitter(height = 0.15),
  # aes(x = Type, y = Segments), alpha = 1/20) +
  scale_y_continuous(limits = c(0, 50)) +
  labs(y = "Number of affixes") +
  scale_x_discrete(labels=c("0" = "Non-monosyllabic", "1" = "Monosyllabic")) +
  theme(axis.title.x=element_blank()) +
  geom_line(data = preds, aes(x = Syllabicity, y = fit, group = 1), linetype = 3)
dev.off()





## Hypothesis 2 - Segments

# Independent variables: Type (binary)
# Dependent variables: Segment length (numeric)
# Random effects: Languages


# Normalize dependent variable
data$Segments_bc <- BoxCox(data$Segments, lambda = BoxCoxLambda(data$Segments))
hist(data$Segments_bc)


# Exploring correlation between intercepts and slopes
ggplot(data, aes(x = Type, y = Segments)) +
  geom_smooth(aes(group = Language), method = "lm", se = FALSE)
# Looks like uncorrelated intercepts and slopes


# Regression model

# REML = FALSE models for model selection
m01_ = lmer(Segments_bc ~ Type + (1 + Type|Language), data, REML = FALSE)

# Alternative model and comparison
m02_ = lmer(Segments_bc ~ 1 + (1 + Type|Language), data, REML = FALSE)
anova(m01_, m02_)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)   
# m02_    5 2475.4 2501.8 -1232.7   2465.4                        
# m01_    6 2469.1 2500.8 -1228.6   2457.1 8.2705  1   0.004029 **

# REML = TRUE model for reporting
m01 = lmer(Segments_bc ~ Type + (1 + Type|Language), data)
m01

# Get confidence intervals for the model
ci = confint(m01)
ci

pp <- profile(m01)

ggplot(as.data.frame(pp),aes(.focal,.zeta))+
  geom_point()+geom_line()+
  facet_wrap(~.par,scale="free_x")+
  geom_hline(yintercept=0,colour="gray")+
  geom_hline(yintercept=c(-1.96,1.96),linetype=2,
             colour="gray")

# Check residuals
plot(m01)
rs = DHARMa::simulateResiduals(m01)
plot(rs)
DHARMa::plotResiduals(rs, form = data$Type, quantreg = T)

# Quick look
plot(allEffects(m01))

# QQplot
qqPlot(residuals(m01))

# Plot random effects
plot_model(m01, ,type="pred", terms=c("Type", "Language"), pred.type="re", ci.lvl = NA)



# Generate predictions for plotting
preds = as.data.frame(predictorEffect("Type", m01))

# Generate predictions for random effect levels
rpreds = as.data.frame(ggpredict(m01, c("Type", "Language"), type = "random"))
rpreds$predicted_nbc = BoxCoxInv(rpreds$predicted, lambda = BoxCoxLambda(data$Segments))
table(rpreds$predicted_nbc, rpreds$group)
rpreds
as.data.frame(rpreds[order(rpreds$group),])

# Un-transform predictions back to segment units
preds$fit_nbc = BoxCoxInv(preds$fit, lambda = BoxCoxLambda(data$Segments))
preds$lower_nbc = BoxCoxInv(preds$lower, lambda = BoxCoxLambda(data$Segments))
preds$upper_nbc = BoxCoxInv(preds$upper, lambda = BoxCoxLambda(data$Segments))

# Create separate dfs for prefixes and suffixes for plotting purposes
p = data[data$Type == "P",]
s = data[data$Type == "S",]

# Plot
tiff("affix_length_by_type.tiff", units="in", width=3, height=5, res=300)
ggplot() +
  geom_pointrange(data = preds, aes(x = Type, y = fit_nbc, ymin = lower_nbc, ymax = upper_nbc)) +
  geom_errorbar(data = preds, aes(x = Type, y = fit_nbc, ymin = lower_nbc, ymax = upper_nbc), width = .05) +
  geom_rug(data = p, sides = "l", position = position_jitter(height = 0.15),
           aes(x = Type, y = Segments), alpha = 1/20) +
  geom_rug(data = s, sides = "r", position = position_jitter(height = 0.15),
           aes(x = Type, y = Segments), alpha = 1/20) +
  scale_y_continuous(limits = c(0.75, 5.25)) +
  labs(y = "Number of segments") +
  scale_x_discrete(labels=c("P" = "Prefixes", "S" = "Suffixes")) +
  theme(axis.title.x=element_blank()) +
  geom_line(data = preds, aes(x = Type, y = fit_nbc, group = 1), linetype = 3)
#  geom_line(data = rpreds, aes(x = x, y = predicted_nbc, group = group), linetype = 3, alpha = 0.5)
#  geom_point(data = rpreds, aes(x = x, y = predicted_nbc, group = group), alpha = 0.5)
dev.off()



## Hypothesis 2 - Syllables

# Independent variables: Type (binary)
# Dependent variables: Segment length (numeric)
# Random effects: Languages


# Normalize dependent variable
data$Syllables_bc <- BoxCox(data$Syllables, lambda = BoxCoxLambda(data$Syllables))
hist(data$Syllables_bc)


# Exploring correlation between intercepts and slopes
ggplot(data, aes(x = Type, y = Syllables)) +
  geom_smooth(aes(group = Language), method = "lm", se = FALSE)
# Looks like uncorrelated intercepts and slopes


# Regression model

# REML = FALSE models for model selection
m51_ = lmer(Syllables_bc ~ Type + (0 + Type|Language), data, REML = FALSE)
m51_

# Alternative model and comparison
m52_ = lmer(Syllables_bc ~ 1 + (0 + Type|Language), data, REML = FALSE)
anova(m51_, m52_)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)   
# m52_    5 2289.0 2314.8 -1139.5   2279.0                        
# m51_    6 2282.5 2313.5 -1135.3   2270.5 8.4659  1   0.003619 **


# REML = TRUE model for reporting
m51 = lmer(Syllables_bc ~ Type + (0 + Type|Language), data)
m51

# Get confidence intervals for the model
ci = confint(m51)

pp <- profile(m51)

ggplot(as.data.frame(pp),aes(.focal,.zeta))+
  geom_point()+geom_line()+
  facet_wrap(~.par,scale="free_x")+
  geom_hline(yintercept=0,colour="gray")+
  geom_hline(yintercept=c(-1.96,1.96),linetype=2,
             colour="gray")

# Check residuals
plot(m51)
rs = DHARMa::simulateResiduals(m51)
plot(rs)
DHARMa::plotResiduals(rs, form = data$Type, quantreg = T)

# Quick look
plot(allEffects(m51))

# QQplot
qqPlot(residuals(m51))

# Plot random effects
plot_model(m51, ,type="pred", terms=c("Type", "Language"), pred.type="re", ci.lvl = NA)



# Generate predictions for plotting
preds = as.data.frame(predictorEffect("Type", m51))

# Generate predictions for random effect levels
rpreds = as.data.frame(ggpredict(m51, c("Type", "Language"), type = "random"))
table(rpreds$predicted_nbc, rpreds$group)
rpreds$predicted_nbc = BoxCoxInv(rpreds$predicted, lambda = BoxCoxLambda(data$Syllables))

# Un-transform predictions back to segment units
preds$fit_nbc = BoxCoxInv(preds$fit, lambda = BoxCoxLambda(data$Syllables))
preds$lower_nbc = BoxCoxInv(preds$lower, lambda = BoxCoxLambda(data$Syllables))
preds$upper_nbc = BoxCoxInv(preds$upper, lambda = BoxCoxLambda(data$Syllables))

# Create separate dfs for prefixes and suffixes for plotting purposes
p = data[data$Type == "P",]
s = data[data$Type == "S",]

# Plot
tiff("syllable_length_by_type.tiff", units="in", width=3, height=5, res=300)
ggplot() +
  geom_pointrange(data = preds, aes(x = Type, y = fit_nbc, ymin = lower_nbc, ymax = upper_nbc)) +
  geom_errorbar(data = preds, aes(x = Type, y = fit_nbc, ymin = lower_nbc, ymax = upper_nbc), width = .05) +
  geom_rug(data = p, sides = "l", position = position_jitter(height = 0.15),
           aes(x = Type, y = Syllables), alpha = 1/20) +
  geom_rug(data = s, sides = "r", position = position_jitter(height = 0.15),
           aes(x = Type, y = Syllables), alpha = 1/20) +
  scale_y_continuous(limits = c(-0.25, 3.25)) +
  labs(y = "Number of syllables") +
  scale_x_discrete(labels=c("P" = "Prefixes", "S" = "Suffixes")) +
  theme(axis.title.x=element_blank()) +
  geom_line(data = preds, aes(x = Type, y = fit_nbc, group = 1), linetype = 3)
#  geom_line(data = rpreds, aes(x = x, y = predicted_nbc, group = group), linetype = 3, alpha = 0.5)
#  geom_point(data = rpreds, aes(x = x, y = predicted_nbc, group = group), alpha = 0.5)
dev.off()




## Hypothesis 3

# Independent variable: Type (prefix vs. suffix)
# Dependent variable: Allomorphy (yes or no)
# Random effect: Language


# Exploring correlation between intercepts and slopes
ggplot(data, aes(x = Type, y = Allomorphy)) +
  geom_smooth(aes(group = Language), method = "lm", se = FALSE)
# Looks like uncorrelated intercepts and slopes


# Regression model

# Initial model
m11 = glmer(Allomorphy ~ Type + (1 + Type|Language), data, family = "binomial")
m11

# Alternative model and comparison
m12 = glmer(Allomorphy ~ 1 + (1 + Type|Language), data, family = "binomial")
anova(m11, m12)
# npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)
# m12    4 1025.2 1046.4 -508.62   1017.2                    
# m11    5 1026.4 1052.8 -508.20   1016.4  0.84  1     0.3594


# Get confidence intervals for the model
ci = confint(m11, method="Wald")
ci

pp <- profile(m11)

ggplot(as.data.frame(pp),aes(.focal,.zeta))+
  geom_point()+geom_line()+
  facet_wrap(~.par,scale="free_x")+
  geom_hline(yintercept=0,colour="gray")+
  geom_hline(yintercept=c(-1.96,1.96),linetype=2,
             colour="gray")

# Quick look
plot(allEffects(m11))

# Check residuals
plot(m11)
rs = DHARMa::simulateResiduals(m11)
plot(rs)
DHARMa::plotResiduals(rs, form = data$Type, quantreg = T)

# Plot random effects
plot_model(m11, ,type="pred", terms=c("Type", "Language"), pred.type="re", ci.lvl = NA)


# Generate predictions for plotting
preds = as.data.frame(predictorEffect("Type", m11))

# Generate predictions for random effect levels
rpreds = as.data.frame(ggpredict(m11, c("Type", "Language"), type = "random"))

# Create separate dfs for prefixes and suffixes for plotting purposes
p = data[data$Type == "P",]
s = data[data$Type == "S",]

# Plot
tiff("allomorphy_by_type.tiff", units="in", width=3, height=5, res=300)
ggplot() +
  geom_pointrange(data = preds, aes(x = Type, y = fit, ymin = lower, ymax = upper)) +
  geom_errorbar(data = preds, aes(x = Type, y = fit, ymin = lower, ymax = upper), width = .05) +
  labs(y = "Odds of Allomorphy") +
  scale_x_discrete(labels=c("P" = "Prefixes", "S" = "Suffixes")) +
  theme(axis.title.x=element_blank()) +
  geom_line(data = preds, aes(x = Type, y = fit, group = 1), linetype = 3) +
  scale_y_log10(breaks = c(0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14))
dev.off()



## Hypothesis 3b

# Independent variables: Allomorphy (binary)
# Dependent variables: Segment length (numeric)
# Random effects: Languages


# Convert Allomorphy variable to factor with appropriate labels
data$Allomorphy2 = data$Allomorphy
data$Allomorphy2[data$Allomorphy2 == 0] <- "No"
data$Allomorphy2[data$Allomorphy2 == 1] <- "Yes"
data$Allomorphy2 = as.factor(data$Allomorphy2)

# Exploring correlation between intercepts and slopes
ggplot(data, aes(x = Allomorphy2, y = Segments_bc)) +
  geom_smooth(aes(group = Language), method = "lm", se = FALSE)
# Looks like uncorrelated intercepts and slopes


# Regression model

# REML = FALSE models for model selection
m21_ = lmer(Segments_bc ~ Allomorphy2 + (0 + Allomorphy2|Language), data, REML = FALSE)

m21_b = glmer(Allomorphy ~ Segments_bc + (0 + Segments_bc|Language), data, family = "binomial")

# Alternative model and comparison
m22_ = lmer(Segments_bc ~ 1 + (0 + Allomorphy2|Language), data, REML = FALSE)
anova(m21_, m22_)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)   
# m02_    5 1767.7 1793.5 -878.85   1757.7                        
# m01_    6 1760.8 1791.8 -874.43   1748.8 8.8372  1   0.002952 **


# REML = TRUE model for reporting
m21 = lmer(Segments_bc ~ Allomorphy2 + (0 + Allomorphy2|Language), data)
m21

# Get confidence intervals for the model
ci = confint(m21)

pp <- profile(m21)

ggplot(as.data.frame(pp),aes(.focal,.zeta))+
  geom_point()+geom_line()+
  facet_wrap(~.par,scale="free_x")+
  geom_hline(yintercept=0,colour="gray")+
  geom_hline(yintercept=c(-1.96,1.96),linetype=2,
             colour="gray")

# Check residuals
plot(m21_b)
rs = DHARMa::simulateResiduals(m21)
plot(rs)
DHARMa::plotResiduals(rs, form = data$Allomorphy, quantreg = T)

# Quick look
plot(allEffects(m21_b))

# QQplot
qqPlot(residuals(m01))

# Plot random effects
plot_model(m21_b, ,type="pred", terms=c("Segments_bc", "Language"), pred.type="re", ci.lvl = NA)



# Generate predictions for plotting
preds = as.data.frame(predictorEffect("Type", m01))

# Generate predictions for random effect levels
rpreds = as.data.frame(ggpredict(m01, c("Type", "Language"), type = "random"))
table(rpreds$predicted_nbc, rpreds$group)
rpreds$predicted_nbc = BoxCoxInv(rpreds$predicted, lambda = BoxCoxLambda(data$Segments))

# Un-transform predictions back to segment units
preds$fit_nbc = BoxCoxInv(preds$fit, lambda = BoxCoxLambda(data$Segments))
preds$lower_nbc = BoxCoxInv(preds$lower, lambda = BoxCoxLambda(data$Segments))
preds$upper_nbc = BoxCoxInv(preds$upper, lambda = BoxCoxLambda(data$Segments))

# Create separate dfs for prefixes and suffixes for plotting purposes
p = data[data$Type == "P",]
s = data[data$Type == "S",]

# Plot
tiff("affix_length_by_type.tiff", units="in", width=3, height=5, res=300)
ggplot() +
  geom_pointrange(data = preds, aes(x = Type, y = fit_nbc, ymin = lower_nbc, ymax = upper_nbc)) +
  geom_errorbar(data = preds, aes(x = Type, y = fit_nbc, ymin = lower_nbc, ymax = upper_nbc), width = .05) +
  geom_rug(data = p, sides = "l", position = position_jitter(height = 0.15),
           aes(x = Type, y = Segments), alpha = 1/20) +
  geom_rug(data = s, sides = "r", position = position_jitter(height = 0.15),
           aes(x = Type, y = Segments), alpha = 1/20) +
  scale_y_continuous(limits = c(0.75, 5.25)) +
  labs(y = "Number of segments") +
  scale_x_discrete(labels=c("P" = "Prefixes", "S" = "Suffixes")) +
  theme(axis.title.x=element_blank()) +
  geom_line(data = preds, aes(x = Type, y = fit_nbc, group = 1), linetype = 3)
#  geom_line(data = rpreds, aes(x = x, y = predicted_nbc, group = group), linetype = 3, alpha = 0.5)
#  geom_point(data = rpreds, aes(x = x, y = predicted_nbc, group = group), alpha = 0.5)
dev.off()




##### Extra stuff


## Trying to troubleshoot random effects structures

# Random intercepts only
m01_ = lmer(Segments_bc ~ Type + (1|Language), data, REML = FALSE)
qqmath(ranef(m01_, condVar = TRUE))
plot_model(m01_, ,type="pred", terms=c("Type", "Language"), pred.type="re", ci.lvl = NA)

# Random slopes only
m02_ = lmer(Segments_bc ~ Type + (0 + Type|Language), data, REML = FALSE)
qqmath(ranef(m02_, condVar = TRUE))
plot_model(m02_, ,type="pred", terms=c("Type", "Language"), pred.type="re", ci.lvl = NA)

# Correlated random intercepts and slopes
m03_ = lmer(Segments_bc ~ Type + (1 + Type|Language), data, REML = FALSE)
qqmath(ranef(m03_, condVar = TRUE))
plot_model(m03_, ,type="pred", terms=c("Type", "Language"), pred.type="re", ci.lvl = NA)

# Uncorrelated random intercepts and slopes
m04_ = lmer(Segments_bc ~ Type + (1|Language) + (0 + Type|Language), data, REML = FALSE)
qqmath(ranef(m04_, condVar = TRUE))
plot_model(m04_, ,type="pred", terms=c("Type", "Language"), pred.type="re", ci.lvl = NA)

anova(m01_, m02_, m03_, m04_)

anova(m03_, m02_)






m21_a = glmer(Allomorphy ~ Segments_bc + (1|Language), data, family = "binomial")
plot(allEffects(m21_a))
plot_model(m21_a, ,type="pred", terms=c("Segments_bc", "Language"), pred.type="re", ci.lvl = NA)

m21_b = glmer(Allomorphy ~ Segments_bc + (0 + Segments_bc|Language), data, family = "binomial")
plot(allEffects(m21_b))
plot_model(m21_b, ,type="pred", terms=c("Segments_bc", "Language"), pred.type="re", ci.lvl = NA)

m21_c = glmer(Allomorphy ~ Segments_bc + (1 + Segments_bc|Language), data, family = "binomial")
plot(allEffects(m21_c))
plot_model(m21_c, ,type="pred", terms=c("Segments_bc", "Language"), pred.type="re", ci.lvl = NA)

m21_d = glmer(Allomorphy ~ Segments_bc + (1|Language) + (0 + Segments_bc|Language), data, family = "binomial")
plot(allEffects(m21_d))
plot_model(m21_d, ,type="pred", terms=c("Segments_bc", "Language"), pred.type="re", ci.lvl = NA)


table(data$Allomorphy2, data$Language)



