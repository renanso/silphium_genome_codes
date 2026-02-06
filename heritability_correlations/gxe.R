gc()
rm(list = ls());ls ()

##################
##Interaction plot
##################
library("ggplot2") # to customize plots
library("emmeans") # TO draw interaction plots
library("dplyr") # for data manipulation
library("ggpubr")

data <- read.csv("final_phenos2.csv", stringsAsFactors=FALSE, fileEncoding="latin1", na.strings = NA)
head(data)
tail(data)
str(data)

# adjusting the factors and variables
data[,1:21] <- lapply(data[,1:21], factor)
data[,22:27] <- lapply(data[,22:27], as.numeric)
str(data)
summary(data)

attach(data)

fm1<- lm(RFC ~ clone * rep * env, data=data)
summary(fm1)
anova(fm1)

fm2<- lm(SM ~ clone * rep * env, data=data)
summary(fm2)
anova(fm2)

fm3<- lm(SNC ~ clone * rep * env, data=data)
summary(fm3)
anova(fm3)

fm4<- lm(RD ~ clone * rep * env, data=data)
summary(fm4)
anova(fm4)

# Calculate mean values for each trait at each environment level
env_means_rfc <- data %>%
  group_by(env) %>%
  summarise(mean_value = mean(RFC, na.rm = TRUE)) %>%
  mutate(env_num = as.numeric(factor(env)),
         trait = "RFC")

env_means_sm <- data %>%
  group_by(env) %>%
  summarise(mean_value = mean(SM, na.rm = TRUE)) %>%
  mutate(env_num = as.numeric(factor(env)),
         trait = "SM")

env_means_snc <- data %>%
  group_by(env) %>%
  summarise(mean_value = mean(SNC, na.rm = TRUE)) %>%
  mutate(env_num = as.numeric(factor(env)),
         trait = "SNC")
env_means_snc<-env_means_snc[c(2,4),]
env_means_snc$env_num<-c(1,2)

env_means_rd <- data %>%
  group_by(env) %>%
  summarise(mean_value = mean(RD, na.rm = TRUE)) %>%
  mutate(env_num = as.numeric(factor(env)),
         trait = "RD")
env_means_rd<-env_means_rd[c(2,4),]
env_means_rd$env_num<-c(1,2)

# Create plots for each trait
plot_rfc <- emmip(fm1, clone~env) + theme(legend.position='none') + 
  ylab("Ray Floret Count") + xlab(NULL) + 
  theme(text = element_text(size = 18)) + 
  scale_colour_manual(values = rep("lightgray", 300)) +
  geom_line(data = env_means_rfc, aes(x = env_num, y = mean_value), 
            colour = "black", linewidth = 1.2, linetype = "solid", inherit.aes = FALSE)

plot_sm <- emmip(fm2, clone~env) + theme(legend.position='none') + 
  ylab("Seed Mass (mg)") + xlab(NULL) + 
  theme(text = element_text(size = 18)) + 
  scale_colour_manual(values = rep("lightgray", 300)) +
  geom_line(data = env_means_sm, aes(x = env_num, y = mean_value), 
            colour = "black", linewidth = 1.2, linetype = "solid", inherit.aes = FALSE)

plot_snc <- emmip(fm3, clone~env) + theme(legend.position='none') + 
  ylab("Seed Number per Capitulum") + xlab(NULL) + 
  theme(text = element_text(size = 18)) + 
  scale_colour_manual(values = rep("lightgray", 300)) +
  geom_line(data = env_means_snc, aes(x = env_num, y = mean_value), 
            colour = "black", linewidth = 1.2, linetype = "solid", inherit.aes = FALSE)

plot_rd <- emmip(fm4, clone~env) + theme(legend.position='none') + 
  ylab("Receptacle Diameter (mm)") + xlab(NULL) + 
  theme(text = element_text(size = 18)) + 
  scale_colour_manual(values = rep("lightgray", 300)) +
  geom_line(data = env_means_rd, aes(x = env_num, y = mean_value), 
            colour = "black", linewidth = 1.2, linetype = "solid", inherit.aes = FALSE)

# Display plots
plot_rfc
plot_sm
plot_snc
plot_rd

tiff("g_x_e_plots.jpeg", width = 30, height = 30, res = 100, units = "cm")
ggarrange(plot_rfc,plot_sm,plot_snc,plot_rd, ncol = 2, nrow = 2)
dev.off()
