#########################################
# Nematodes in a Niayes Gardens (NIANG) #
# Functional group analysis              #
# Code by: Stephen Wood                 #
# Last updated: August 2018             #
#########################################

#### LOAD PACKAGES ####
library(tidyverse)
library(readxl)
library(lme4)


#### DEFINE CUSTOM FUNCTIONS ####
st.err <- function(x) {
  sd(x)/sqrt(length(x))
}

#### READ IN DATA ####
site_id <- read_excel("sn_nematodes.xls", sheet = "site_id")
taxa_id <- read_excel("sn_nematodes.xls", sheet = "taxa_id")
id <- read_excel("sn_nematodes.xls", sheet = "identification")

#### MANIPULATE DATA ####
id_trans <- id %>% 
  gather(var, val, 2:ncol(id)) %>% 
  spread(Taxons, val) %>%
  select(-c(TOTAL,Actinolaimus))
names(id_trans)[1] <- 'SN'
id_trans$SN <- as.numeric(id_trans$SN)

# Join database with site information level
id_merge <- full_join(id_trans, site_id) %>% select(SN:Lieu, SOL:Input_level)
names(id_merge)[51:56] <- c('site','soil','treatment','treatment_1','treatment_2','input_level')

# Drop observations with double addition because of lack of replication
# Drop observations with 'procep' taxon because of lack of functional classification
id_merge <- filter(id_merge, is.na(input_level) == TRUE) %>% select(-c(input_level, procep))

# Convert from wide to long
id_merge_long <- gather(id_merge, var, val, Achromodora:Zeldia)

# Relevel factors for analysis
id_merge_long$soil <- factor(id_merge_long$soil, levels=c('Bulk Soil', 'Rhizosphere'))

# Add in taxonomic functional groupings and remove NA values
names(id_merge_long)[7] <- 'Taxons'
id_all <- full_join(id_merge_long, taxa_id) %>% filter(is.na(val) != TRUE)


#### VISUALIZE DATA ####
sum <- aggregate(id_all$val, by = list(id_all$SN, id_all$soil, id_all$Group), FUN = sum)
mean <- aggregate(sum$x, by = list(sum$Group.2, sum$Group.3), FUN = mean)
se <- aggregate(sum$x, by = list(sum$Group.2, sum$Group.3), FUN = st.err)

plot <- mean
plot$mean <- mean$x
plot$upper <- mean$x + se$x
plot$lower <- mean$x - se$x
plot <- select(plot, Group.1, Group.2, mean, upper, lower)
names(plot)[1:2] <- c("soil","Group")

ggplot(plot, aes(x = Group, y = mean/100, fill=soil)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=lower/100, ymax=upper/100), position=position_dodge(width=0.9), width=0.25) +
  ylab("Nematode abundance\n(number per g soil)\n") +  
  xlab('') + theme_bw() + 
  scale_fill_manual(values=c('#1b9e77','#d95f02')) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 11)
  )

sum <- aggregate(id_all$val, by = list(id_all$SN, id_all$treatment_2, id_all$Group), FUN = sum)
mean <- aggregate(sum$x, by = list(sum$Group.2, sum$Group.3), FUN = mean)
se <- aggregate(sum$x, by = list(sum$Group.2, sum$Group.3), FUN = st.err)

plot <- mean
plot$mean <- mean$x
plot$upper <- mean$x + se$x
plot$lower <- mean$x - se$x
plot <- select(plot, Group.1, Group.2, mean, upper, lower)
names(plot)[1:2] <- c("Treatment","Group")

ggplot(plot, aes(x = Group, y = mean/100, fill=Treatment)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=lower/100, ymax=upper/100), position=position_dodge(width=0.9), width=0.25) +
  ylab("Nematode abundance\n(number per g soil)\n") +  
  xlab('') + theme_bw() + 
  scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c')) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 11)
  )


#### CONDUCT ANALYSES ####
# Does abundance of functional groups vary by treatments? #
summary(glmer(abund ~ soil + treatment_2 + (1|site), family=gaussian(link="log"), data=abun_and_id))

