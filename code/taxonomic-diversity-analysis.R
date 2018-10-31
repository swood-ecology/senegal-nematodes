#########################################
# Nematodes in a Niayes Gardens (NIANG) #
# Taxonomic group analysis              #
# Code by: Stephen Wood                 #
# Last updated: September 2018          #
#########################################

#### LOAD PACKAGES ####
library(tidyverse)
library(readxl)
library(lme4)


#### DEFINE CUSTOM FUNCTIONS ####
st.err <- function(x) {
  sd(x,na.rm=T)/sqrt(length(x[!is.na(x)]))
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
rm(id_trans)

# Drop observations with double addition because of lack of replication
id_merge <- filter(id_merge, is.na(input_level) == TRUE) 

# Convert from wide to long
id_merge_long <- gather(id_merge, var, val, Achromodora:Zeldia)
rm(id_merge)

# Relevel factors for analysis
id_merge_long$soil <- factor(id_merge_long$soil, levels=c('Bulk Soil', 'Rhizosphere'))

# Add in taxonomic functional groupings and remove NA values
names(id_merge_long)[8] <- 'Taxons'
id_all <- filter(id_merge_long, is.na(val) != TRUE)
rm(id_merge_long)

# Add total abundance variable per plot level
agg <- aggregate(id_all$val, by=list(id_all$site,id_all$soil,id_all$treatment),FUN=sum) 
names(agg) <- c('site','soil','treatment','total_abun')
id_all <- full_join(id_all,agg)
rm(agg)

# Subset top taxa
id_all$ra <- id_all$val / id_all$total_abun
test <- aggregate(id_all$val, by = list(id_all$Taxons), FUN = mean)
test[order(test$x),]
taxa <- c('Rotylenchulus','Aphelenchus','Tylenchidea','Acrobeloides','Meloidogyne','Ditylenchus','Tylenchorynchus',
          'Acrobeles','Prodorylaimus','Zeldia','Cephalobus','Chiloplacus')
id_all <- id_all[which(id_all$Taxons %in% taxa),]
rm(test, taxa)

#### VISUALIZE DATA ####
sum <- aggregate(id_all$val, by = list(id_all$SN, id_all$soil, id_all$Taxons), FUN = sum)
mean <- aggregate(sum$x, by = list(sum$Group.2, sum$Group.3), FUN = mean)
se <- aggregate(sum$x, by = list(sum$Group.2, sum$Group.3), FUN = st.err)

plot <- mean
plot$mean <- mean$x
plot$upper <- mean$x + se$x
plot$lower <- mean$x - se$x
plot <- select(plot, Group.1, Group.2, mean, upper, lower)
names(plot)[1:2] <- c("soil","Taxa")

ggplot(plot, aes(x = Taxa, y = mean/100, fill=soil)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=lower/100, ymax=upper/100), position=position_dodge(width=0.9), width=0.25) +
  ylab("Nematode total abundance\n(number per g soil)\n") +  
  xlab('') + theme_bw() + 
  scale_fill_manual(values=c('#1b9e77','#d95f02')) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 13, angle=45, hjust=1),
    axis.text.y = element_text(size = 11)
  )

sum <- aggregate(id_all$val, by = list(id_all$SN, id_all$treatment_2, id_all$Taxons), FUN = sum)
mean <- aggregate(sum$x, by = list(sum$Group.2, sum$Group.3), FUN = mean)
se <- aggregate(sum$x, by = list(sum$Group.2, sum$Group.3), FUN = st.err)

plot <- mean
plot$mean <- mean$x
plot$upper <- mean$x + se$x
plot$lower <- mean$x - se$x
plot <- select(plot, Group.1, Group.2, mean, upper, lower)
names(plot)[1:2] <- c("Treatment","Taxa")

ggplot(plot, aes(x = Taxa, y = mean/100, fill=Treatment)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=lower/100, ymax=upper/100), position=position_dodge(width=0.9), width=0.25) +
  ylab("Nematode total abundance\n(number per g soil)\n") +  
  xlab('') + theme_bw() + 
  scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c')) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 13, angle=45, hjust=1),
    axis.text.y = element_text(size = 11)
  )


#### CONDUCT ANALYSES ####
# Does abundance of functional groups vary by treatments? #
summary(glmer(abund ~ soil + treatment_2 + (1|site), family=gaussian(link="log"), data=abun_and_id))

