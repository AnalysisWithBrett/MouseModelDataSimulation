#Simulation for mice project
rm(list=ls())

#Library
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggeffects)
library(tidyverse)
library(lmerTest)
library(lme4)
library(cowplot)
library(MASS)
library(broom)

#Importing data
mice <- read.table(file.choose(), stringsAsFactors = TRUE, header=TRUE) #use this


#Gives random concentrations to each cohort
mice2 <- mice %>% 
  mutate(iga0 = case_when(cohort == 'wild' ~ c(runif(48, min=70, max=400)),
                          cohort == 'labling' ~ c(runif(48, min=0, max=150)),
                          cohort == 'lab' ~ c(runif(48, min=0, max=150))),
         
         iga7 = case_when(cohort == 'wild' ~ c(runif(48, min=70, max=400)),
                          cohort == 'labling' ~ c(runif(48, min=0, max=150)),
                          cohort == 'lab' ~ c(runif(48, min=0, max=150))),
         
         iga14 = case_when(cohort == 'wild' ~ c(runif(48, min=20, max=200)),
                          cohort == 'labling' ~ c(runif(48, min=0, max=70)),
                          cohort == 'lab' ~ c(runif(48, min=0, max=70))),
         
         iga21 = case_when(cohort == 'wild' ~ c(runif(48, min=1, max=70)),
                          cohort == 'labling' ~ c(runif(48, min=0, max=40)),
                          cohort == 'lab' ~ c(runif(48, min=0, max=40))),
         
         iga28 = case_when(cohort == 'wild' ~ c(runif(48, min=70, max=400)),
                          cohort == 'labling' ~ c(runif(48, min=0, max=150)),
                          cohort == 'lab' ~ c(runif(48, min=0, max=150))))

#Converts short format to long format and removing letters from 'days'
mice_long <- mice2 %>% 
  gather(key = day, value = concentration, c(3:7)) %>% 
  mutate(day = parse_number(as.character(day)))

#converts numeric value into factor for ggplot boxplot
mice_long$day <- factor(mice_long$day, levels=c("0", "7", "14", "21", "28"),
                                     labels = c("0", "7", "14", "21", "28"))

#Data visualisation- creating function for background
theme.my.own <- function(){
  theme_bw()+
    theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 12, angle = 45),
          axis.title.x = element_text(size = 14, face = "plain"),             
          axis.title.y = element_text(size = 14, face = "plain"),             
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size = 20, vjust = 1, hjust = 0.5),
          legend.text = element_text(size = 12, face = "italic"),          
          legend.title = element_blank(),                              
          legend.position = c(0.9, 0.9))
}

mice_comp <- list( c("0", "7"), c("0", "14"), c("0", "21"), c("0", "28"), c("7", "14"), c("7", "21"),
                   c("7", "28"),c("14", "21"), c("14", "28"), c("0", "28"))
#Creating colour palette
magic.palette <- c("#f03b20", "#1f78b4", "#2ca25f")
names(magic.palette) <- levels(mice_long$cohort)  #Linking palette with cohort   

magic.palette2 <- c("#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c")
names(magic.palette2) <- levels(mice_long$day)  #Linking palette with days

#Visualising the whole data with boxplots without comparing between cohorts
ggplot(mice_long, aes(x = day,y = concentration, fill = day)) +
  geom_boxplot(aes(day, concentration)) +
  geom_vline(xintercept=c(3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9), lwd = 10, color= "rosybrown1") +
  geom_boxplot(aes(day, concentration), lwd=0.7) +
  xlab("Time (day)") + ylab("IgA concentration (mg/kg)") +
  theme.my.own() + theme(legend.position = c(0.9, 0.8)) +
  scale_fill_manual(values = magic.palette2) +                
  scale_colour_manual(values = magic.palette2) +
  stat_summary(fun.y = mean)

# Making a histogram to assess data distribution
(mice.hist <- ggplot(mice_long, aes(concentration)) + 
    geom_histogram(aes(x = concentration),
                   alpha = 0.9,
                   fill = 'rosybrown1') +   
    labs(x = 'Day',
         y = 'Concentration (normal)',
         title = 'Concentration of IgA') + 
    theme.my.own())

#Fitting a lm for concentration against day
mod <- lm(concentration ~ day, data = mice_long)
summary(mod)
anova(mod)
plot(mod, which = 1)

#Need to transform the data
# Log transform the data
mice_long <- mice_long %>%
  mutate(logconcentration = log(concentration))

#Plotting log transformed data on histogram
(mice_log <- ggplot(mice_long, aes(logconcentration)) + 
    geom_histogram(aes(x = logconcentration),
                   alpha = 0.9,
                   fill = "rosybrown1") + 
    labs(x = 'Day',
         y = 'Concentration (log transformed)',
         title = 'Concentration of IgA') +
    theme.my.own())

#Fitting lm with log data
logmod <- lm(logconcentration ~ day, data = mice_long)
summary(logmod)
anova(logmod)
plot(logmod)

# Find the optimal lambda for Box-Cox
boxcox <- boxcox(mod)
boxcox

## Extract the optimal lambda value
(lambda <- boxcox$x[which.max(boxcox$y)])


# Transform the data using this lambda value
mice_long <- mice_long %>%
  mutate(bcconcentration = ((concentration^lambda-1)/lambda))


# Plot a histogram of the Box-Cox transformed data
(mice_bc <- ggplot(data = mice_long) +
    geom_histogram(aes(x = bcconcentration),
                   alpha = 0.9,
                   fill = 'rosybrown1') +
    labs(x = 'Day',
         y = 'Concentration (BC transformed)',
         title = 'Concentration of IgA') +
    theme.my.own())

#Fitting lm with boxcox data
bcmod <- lm(bcconcentration ~ day, data = mice_long)
summary(bcmod)
anova(bcmod)
plot(bcmod)

# Verify reverse transformations by creating new columns that should match the original
mice_long <- mice_long %>%
  mutate(back_log = exp(logconcentration),
         back_bc = (bcconcentration*lambda + 1)^(1/lambda))

#Comparing between normal and the two transformations
mice_long %>%
  dplyr::select(concentration, back_log, back_bc)

#Checking the datasheet so far
View(mice_long)


###Testing boxcox with interaction###
#Fitting a lm for concentration against day with an interaction and mouse as random effects
int.conc <- lm(concentration ~ day * cohort, data = mice_long)
anova(int.conc)
summary(int.conc)
plot(int.conc, which = 1)

# Find the optimal lambda for Box-Cox
boxcox_int <- boxcox(int.conc)
boxcox_int

## Extract the optimal lambda value
(lambda_int <- boxcox_int$x[which.max(boxcox_int$y)])


# Transform the data using this lambda value
mice_long <- mice_long %>%
  mutate(bcconcentration_int = ((concentration^lambda_int-1)/lambda_int))


# Plot a histogram of the Box-Cox transformed data
(mice_bc <- ggplot(data = mice_long) +
    geom_histogram(aes(x = bcconcentration_int),
                   alpha = 0.9,
                   fill = 'rosybrown1') +
    labs(x = 'Day',
         y = 'Concentration (BC transformed)',
         title = 'Concentration of IgA') +
    theme.my.own())

#Fitting lm with boxcox data
bcmod_int <- lm(bcconcentration_int ~ day * cohort, data = mice_long)
summary(bcmod_int)
anova(bcmod_int)
plot(bcmod_int, which = 1)

#Visualising with boxcox data
(mice_boxplot_int <- ggplot(mice_long, aes(x = day,y = bcconcentration_int, fill = cohort)) +
    geom_boxplot(aes(day, bcconcentration_int)) +
    geom_vline(xintercept=c(2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.1, 4.2, 4.3), 
               lwd = 5, color= "rosybrown1") +
    geom_boxplot(aes(day, bcconcentration_int)) +
    xlab("Time (day)") + ylab("IgA concentration (mg/kg)") +
    theme.my.own() +
    scale_fill_manual(values = magic.palette))

#Visualising with normal data
(mice_boxplot <- ggplot(mice_long, aes(x = day,y = concentration, fill = cohort)) +
  geom_boxplot(aes(day, concentration)) +
  geom_vline(xintercept=c(2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.1, 4.2, 4.3), 
             lwd = 5, color= "rosybrown1") +
  geom_boxplot(aes(day, concentration)) +
  xlab("Time (day)") + ylab("IgA concentration (mg/kg)") +
  theme.my.own() +
  scale_fill_manual(values = magic.palette))

# Plot a histogram of the Box-Cox transformed data
(mice_bc_int <- ggplot(data = mice_long) +
    geom_histogram(aes(x = bcconcentration_int),
                   alpha = 0.9,
                   fill = 'rosybrown1') +
    labs(x = 'Day',
         y = 'Concentration (BC transformed)',
         title = 'Concentration of IgA') +
    theme.my.own())



# Panel of histograms for different transformations
(mice_dist_panel <- plot_grid(mice.hist + labs(title = 'Original data'),  # original data  
                              mice_log + labs(title = 'Log transformation'),  # logarithmic transformation
                              mice_bc + labs(title = 'Box-Cox transformation'),  # Box-Cox transformation
                              mice_bc_int + labs(title = "Box-cox with interaction"),
                              nrow = 2,  # number of row in the panel
                              ncol = 2))  # number of columns in the panel














#Fitting a lm for concenrtation against day with an interaction and mouse as random effects

### UNDER CONSTRUCTION ###

lmer.conc <- lmer(concentration ~ day*cohort + (1|mouse), data = mice_long)
vcov(lmer.conc) #normal
anova(lmer.conc)
summary(lmer.conc)

lmer.log <- lmer(logconcentration ~ day*cohort + (1|mouse), data = mice_long)
vcov(lmer.log) #log transformed
anova(lmer.log)
mice_long
lmer.bc <- lmer(bcconcentration_int ~ day*cohort + (1|mouse), data = mice_long)
summary(lmer.bc) #boxcox transformed


coef(lmer.bc)

#Checking for normality
par(mfrow = c(2,2))
qqnorm(resid(lmer.conc))
qqline(resid(lmer.conc))

qqnorm(resid(lmer.log))
qqline(resid(lmer.log))

qqnorm(resid(lmer.bc))
qqline(resid(lmer.bc))

#Checking for homodescasticity
par(mfrow = c(2,2))
varconc <- plot(lmer.conc, which = 1)

varlog <- plot(lmer.log, which = 1)

varbc <- plot(lmer.bc, which = 1)

# Panel for homogeneity
(var_panel <- plot_grid(varconc, 
                          varlog,
                          varbc,
                          legend,
                          nrow = 2,  # number of row in the panel
                          ncol = 2))  # number of columns in the panel


#Run multiple analyses
(split_plot_cohort <- ggplot(aes(day, concentration, fill=day), data = mice_long) + 
    geom_boxplot(aes(day,concentration)) +
    scale_fill_manual(values = magic.palette2) +                
    scale_colour_manual(values = magic.palette2) +
    facet_wrap(~ cohort) + # create a facet for each mountain range
    xlab("Day") + 
    ylab("Concentration") + theme.my.own() +
    theme(legend.position = 'below'))

(split_plot_day <- ggplot(aes(cohort, concentration, fill=cohort), data = mice_long) + 
    geom_boxplot(aes(cohort, concentration)) +
    scale_fill_manual(values = magic.palette) +                
    scale_colour_manual(values = magic.palette) +
    facet_wrap(~ day) + # create a facet for each mountain range
    xlab("Cohort") + 
    ylab("Concentration") + theme.my.own() +
    theme(legend.position = 'below'))

#calculate summary statistics including our sample size (n), mean hatching time per temperature level, 
#standard deviation and standard error values
(summary_mice <- mice_long %>%
  group_by(day,cohort) %>%
  summarise(n = n(),  # Calculating sample size n
            average_concentration = mean(concentration),  
            # Calculating mean hatching time
            SD = sd(concentration))%>%  # Calculating standard deviation
  mutate(SE = SD / sqrt(n)))  # Calculating standard error

# Making a barplot
dodge <- position_dodge(width = 0.9)

(mice_barplot <- ggplot(aes(day, average_concentration, fill=cohort), data = summary_mice) +                                  
    geom_bar(aes(x = day, y = average_concentration),
             stat = "identity", position = position_dodge()) +               
    geom_errorbar(aes(x = day, ymin = average_concentration - SE,
                      ymax = average_concentration + SE), width = 0.2,
                  colour="black", alpha=0.9,size=1, position = dodge) +  
    # Adding standard error bars        
    scale_fill_manual(values = magic.palette) +                
    scale_colour_manual(values = magic.palette) +
    labs(x = "Day", y = "IgA Concentration (mg/kg)") +
    theme.my.own())

#changing days back to numeric
mice_long2 <- mice_long %>% mutate(day = as.numeric(as.character(day)))

#with quadratic equation rather than linear for normal data
(mice_scatter <- ggplot(mice_long2, aes(day, concentration, colour = cohort)) +
    geom_vline(xintercept=c(14:21), lwd = 8, color= "rosybrown1") +
    geom_jitter(width = 0.7) +
    scale_fill_manual(values = magic.palette) +                
    scale_colour_manual(values = magic.palette) +
    theme.my.own() + scale_x_continuous(breaks = seq(0, 28, by = 7)) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2) + I(x^3), se = FALSE) +
    theme(legend.position = 'right') + xlab("Day") + ylab("Concentration of Iga (µg/ml)"))

(mice_scatter <- ggplot(mice_long2, aes(day, bcconcentration, colour = cohort)) +
    geom_vline(xintercept=c(14:21), lwd = 8, color= "rosybrown1") +
    geom_jitter(width = 0.7) +
    scale_fill_manual(values = magic.palette) +                
    scale_colour_manual(values = magic.palette) +
    theme.my.own() + scale_x_continuous(breaks = seq(0, 28, by = 7)) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2) + I(x^3), se = FALSE))

# Panel of different graphs
(graph_panel <- plot_grid(mice_boxplot, 
                            mice_barplot,
                            mice_scatter,
                            nrow = 2,  # number of row in the panel
                            ncol = 2))  # number of columns in the panel


#with quadratic equation rather than linear for log plot
(mice_scatter_log <- ggplot(mice_long2, aes(day, logconcentration, colour = cohort)) +
    geom_vline(xintercept=c(14:21), lwd = 8, color= "rosybrown1") +
    geom_jitter(width = 0.7) +
    scale_fill_manual(values = magic.palette) +                
    scale_colour_manual(values = magic.palette) +
    theme.my.own() + scale_x_continuous(breaks = seq(0, 28, by = 7)) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE))

#with quadratic equation rather than linear for boxcox plot
(mice_scatter_bc <- ggplot(mice_long2, aes(day, bcconcentration, colour = cohort)) +
    geom_vline(xintercept=c(14:21), lwd = 8, color= "rosybrown1") +
    geom_jitter(width = 0.7) +
    scale_fill_manual(values = magic.palette) +                
    scale_colour_manual(values = magic.palette) +
    theme.my.own() + scale_x_continuous(breaks = seq(0, 28, by = 7)) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2) + I(x^3), se = FALSE))

# Panel for scatter plot
(graph_panel <- plot_grid(mice_scatter, 
                          mice_scatter_log,
                          mice_scatter_bc,
                          nrow = 2,  # number of row in the panel
                          ncol = 2))  # number of columns in the panel





















###EXTRA###
#Quadratic regression for normal data

#Adding a new column containing day to the power of 2 for quadratic regression
mice_long2 <- mice_long2 %>% mutate(day2 = day^2)

#Performing quadratic regression for normal data
quadraticmodconc <- lm(concentration ~ (day + day2) * cohort, data = mice_long2)
summary(quadraticmodconc)
plot(quadraticmodconc)

#Creating a linear model for comparison
linearmodconc <- lm(concentration ~ day * cohort, data=mice_long2)
anova(linearmodconc)
plot(linearmodconc)

#Visualising the quadratic regression
(mice_scatter <- ggplot(mice_long2, aes(day, concentration, fill = cohort)) +
    geom_vline(xintercept=c(14:21), lwd = 8, color= "rosybrown1") +
    geom_jitter(width = 0.7, aes(color= cohort)) +
    scale_fill_manual(values = magic.palette) +                
    scale_colour_manual(values = magic.palette) +
    theme.my.own() + scale_x_continuous(breaks = seq(0, 28, by = 7)) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, aes(color = cohort)) +
    stat_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = cohort)) +
    stat_summary(fun.y = mean, colour = "black"))


anova(linearmodconc, quadraticmodconc)


#Quadratic regression for log data

#Performing quadratic regression for log data
quadraticmodlog <- lm(logconcentration ~ (day + day2), data = mice_long2)
summary(quadraticmodlog)
plot(quadraticmodlog)

#Creating a linear model for comparison
linearmodlog <- lm(logconcentration ~ day, data=mice_long2)
summary(linearmodlog)
plot(linearmodlog)


#Visualising the quadratic regression
(mice_scatter <- ggplot(mice_long2, aes(day, logconcentration, fill= cohort)) +
    geom_vline(xintercept=c(14:21), lwd = 8, color= "rosybrown1") +
    geom_jitter(width = 0.7, aes(color= cohort)) +
    scale_fill_manual(values = magic.palette) +                
    scale_colour_manual(values = magic.palette) +
    theme.my.own() + scale_x_continuous(breaks = seq(0, 28, by = 7)) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, aes(color = cohort)) +
    stat_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = cohort)) +
    stat_summary(fun.y = mean, colour = "black"))


anova(linearmodlog, quadraticmodlog)


#Quadratic regression for boxcox

#Performing quadratic regression for boxcox data
quadraticmodbc <- lm(bcconcentration ~ (day + day2) * cohort, data = mice_long2)
anova(quadraticmodbc)
plot(quadraticmodbc)

#Creating a linear model for comparison
linearmodbc <- lm(bcconcentration ~ day * cohort, data=mice_long2)
anova(linearmodbc)
plot(linearmodbc)


#Visualising the quadratic regression
(mice_scatter <- ggplot(mice_long2, aes(day, bcconcentration, fill = cohort)) +
    geom_vline(xintercept=c(14:21), lwd = 8, color= "rosybrown1") +
    geom_jitter(width = 0.7, aes(color= cohort)) +
    scale_fill_manual(values = magic.palette) +                
    scale_colour_manual(values = magic.palette) +
    theme.my.own() + scale_x_continuous(breaks = seq(0, 28, by = 7)) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, aes(color = cohort)) +
    stat_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = cohort)) +
    stat_summary(fun.y = mean, colour = "black"))

anova(linearmodbc, quadraticmodbc)

#Comparing interactions between normal and boxcox linear model
anova(linearmodbc)
anova(linearmodconc)

#Comparing interactions between normal and boxcox quadratic model
anova(quadraticmodconc)
anova(quadraticmodbc)

tempfile



