#########################################
# Nematodes in a Niayes Gardens (NIANG) #
# Abundance analysis                    #
# Code by: Stephen Wood                 #
# Last updated: August 2018             #
#########################################

#### LOAD PACKAGES ####
library(tidyverse)
library(readxl)
library(lme4)

#### READ IN DATA ####
site_id <- read_excel("sn_nematodes.xls", sheet = "site_id")
abundance <- read_excel("sn_nematodes.xls", sheet = "abundance")

#### MANIPULATE DATA ####
abun_and_id <- full_join(abundance, site_id) %>% select(`nb/100g sol`:Input_level)
names(abun_and_id) <- c('abund','site','block','soil','treatment','treatment_1','treatment_2','input_level')

# Drop observations with double addition because of lack of replication
abun_and_id <- filter(abun_and_id, is.na(input_level) == TRUE) %>% select(-input_level)

# Relevel factors for analysis
abun_and_id$soil <- factor(abun_and_id$soil, levels=c('Bulk Soil', 'Rhizosphere'))
abun_and_id$treatment <- factor(abun_and_id$treatment, levels=c('Témoin absolu', 'Témoin NPK', 'Fumier Cheval', 'Fumier Cheval + NPK', 'Fumier Cheval double', 
                                                                'Fumier Cheval double +NPK', 'Compost', 'Compost + NPK', 'Compost double', 'Compost double+NPK'))

# Convert abundance to per gram of soil
abun_and_id$abund <- abun_and_id$abund / 100

#### VISUALIZE DATA ####
boxplot(abund ~ soil, data=abun_and_id)
boxplot(abund ~ treatment, data=abun_and_id)
boxplot(abund ~ treatment_1, data=abun_and_id)
boxplot(abund ~ treatment_2, data=abun_and_id)

ggplot(abun_and_id,aes(x=soil, y=abund)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.7) +
  ylab("Nematode abundance\n(number per g soil)\n") + scale_y_continuous(breaks = pretty(abun_and_id$abund, n = 8)) + 
  xlab('') + theme_bw() + 
  theme(
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 11)
  )

ggplot(abun_and_id,aes(x=treatment_2, y=abund)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.7) +
  ylab("Nematode abundance\n(number per g soil)\n") + scale_y_continuous(breaks = pretty(abun_and_id$abund, n = 8)) + 
  xlab('') + theme_bw() + 
  theme(
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 11)
  )


#### CONDUCT ANALYSES ####
# Does abundance vary by treatments? #
summary(glmer(abund ~ soil + treatment_2 + (1|site), family=gaussian(link="log"), data=abun_and_id))

