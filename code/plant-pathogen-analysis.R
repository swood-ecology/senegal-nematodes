#########################################
# Nematodes in a Niayes Gardens (NIANG) #
# Plant pathogen analysis               #
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
taxa_id$troph_guild <- paste(taxa_id$Group, taxa_id$`life-history-degree`)

joined <- full_join(id, taxa_id) %>% 
  gather(var, val, 2:ncol(id)) %>%
  filter(is.na(val)!=TRUE)

join_groups <- aggregate(joined$val, by=list(joined$troph_guild, joined$var), FUN=sum) %>%
  spread(Group.1,x)
names(join_groups)[1] <- 'SN'

join_groups$c.p1 <- join_groups$`Bacteriovore 1`
join_groups$c.p2 <- join_groups$`Bacteriovore 2`+join_groups$`Fungivore 2`
join_groups$c.p3 <- join_groups$`Bacteriovore 3`
join_groups$c.p4 <- join_groups$`Bacteriovore 4`+join_groups$`Fungivore 4`+join_groups$`Carnivore 4`+join_groups$`Omnivore 4`
join_groups$c.p5 <- join_groups$`Carnivore 5`+join_groups$`Omnivore 5`
join_groups$total_free <- join_groups$`Bacteriovore 1`+join_groups$`Bacteriovore 2`+join_groups$`Bacteriovore 3`+join_groups$`Bacteriovore 4`+join_groups$`Fungivore 2`+join_groups$`Fungivore 4`+join_groups$`Omnivore 4`+join_groups$`Omnivore 5`+join_groups$`Carnivore 4`+join_groups$`Carnivore 5`
join_groups$Cp1 <- join_groups$c.p1/join_groups$total_free
join_groups$Cp2 <- join_groups$c.p2/join_groups$total_free
join_groups$Cp3 <- join_groups$c.p3/join_groups$total_free
join_groups$Cp4 <- join_groups$c.p4/join_groups$total_free
join_groups$Cp5 <- join_groups$c.p5/join_groups$total_free

join_groups$SN <- as.numeric(join_groups$SN)

paths <- full_join(join_groups, site_id) %>% 
  select(SN, c.p1:c.p5,Lieu:Input_level) %>% 
  gather(type, value, c.p1:c.p5) %>% 
  filter(is.na(Input_level) == TRUE) %>% 
  select(-Input_level)


#### VISUALIZE DATA ####
sum <- aggregate(paths$value, by = list(paths$SN, paths$SOL, paths$type), FUN = sum)
mean <- aggregate(sum$x, by = list(sum$Group.2, sum$Group.3), FUN = mean, na.rm=T)
se <- aggregate(sum$x, by = list(sum$Group.2, sum$Group.3), FUN = st.err)

plot <- mean
plot$mean <- mean$x
plot$upper <- mean$x + se$x
plot$lower <- mean$x - se$x
plot <- select(plot, Group.1, Group.2, mean, upper, lower)
names(plot)[1:2] <- c("Treatment","Type")

ggplot(plot, aes(x = Type, y = mean, fill=Treatment)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), position=position_dodge(width=0.9), width=0.25) +
  ylab("Plant pathogen group abundance\n(number per g soil)\n") +  
  xlab('') + theme_bw() + 
  scale_fill_manual(values=c('#1b9e77','#d95f02')) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 13, angle=45, hjust=1),
    axis.text.y = element_text(size = 11)
  )

sum <- aggregate(paths$value, by = list(paths$SN, paths$Treatment_2, paths$type), FUN = sum)
mean <- aggregate(sum$x, by = list(sum$Group.2, sum$Group.3), FUN = mean, na.rm=T)
se <- aggregate(sum$x, by = list(sum$Group.2, sum$Group.3), FUN = st.err)

plot <- mean
plot$mean <- mean$x
plot$upper <- mean$x + se$x
plot$lower <- mean$x - se$x
plot <- select(plot, Group.1, Group.2, mean, upper, lower)
names(plot)[1:2] <- c("Treatment","Type")

ggplot(plot, aes(x = Type, y = mean, fill=Treatment)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), position=position_dodge(width=0.9), width=0.25) +
  ylab("Plant pathogen group abundance\n(number per g soil)\n") +  
  xlab('') + theme_bw() + 
  scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c')) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 13, angle=45, hjust=1),
    axis.text.y = element_text(size = 11)
  )
