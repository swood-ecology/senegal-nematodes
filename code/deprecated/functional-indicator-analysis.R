#########################################
# Nematodes in a Niayes Gardens (NIANG) #
# Functional indicator analysis         #
# Code by: Stephen Wood                 #
# Last updated: August 2018             #
#########################################

#### LOAD PACKAGES ####
library(tidyverse)
library(readxl)
library(lme4)


#### READ IN DATA ####
site_id <- read_excel("sn_nematodes.xls", sheet = "site_id")
taxa_id <- read_excel("sn_nematodes.xls", sheet = "taxa_id")
id <- read_excel("sn_nematodes.xls", sheet = "identification")


#### DEFINE CUSTOM FUNCTIONS ####
st.err <- function(x) {
  sd(x)/sqrt(length(x))
}


#### MANIPULATE DATA ####
taxa_id$troph_guild <- paste(taxa_id$Group, taxa_id$`life-history-degree`)

joined <- full_join(id, taxa_id) %>% 
  gather(var, val, 2:ncol(id)) %>%
  filter(is.na(val)!=TRUE)

join_groups <- aggregate(joined$val, by=list(joined$troph_guild, joined$var), FUN=sum) %>%
  spread(Group.1,x)
names(join_groups)[1] <- 'SN'

join_groups$Ba1 <- join_groups$`Bacteriovore 1` * 3.2
join_groups$Ba2 <- join_groups$`Bacteriovore 2` * 0.8
join_groups$Ba3 <- join_groups$`Bacteriovore 3` * 1.8
join_groups$Ba4 <- join_groups$`Bacteriovore 4` * 3.2
join_groups$Fu2 <- join_groups$`Fungivore 2` * 0.8
join_groups$Fu4 <- join_groups$`Fungivore 4` * 3.2
join_groups$Om4 <- join_groups$`Omnivore 4` * 3.2
join_groups$Om5 <- join_groups$`Omnivore 5` * 5
join_groups$Ca4 <- join_groups$`Carnivore 4` * 3.2
join_groups$Ca5 <- join_groups$`Carnivore 5` * 5
join_groups$PP2 <- join_groups$`Phytoparasite 2` * 0.8
join_groups$PP3 <- join_groups$`Phytoparasite 3` * 1.8
join_groups$PP5 <- join_groups$`Phytoparasite 5` * 5

join_groups$b <- join_groups$Ba2 + join_groups$Fu2
join_groups$e <- join_groups$Ba1 + join_groups$Fu2
join_groups$s <- join_groups$Ba3 + join_groups$Ba4 + join_groups$Fu4 + 
  join_groups$Om4 + join_groups$Om5 + join_groups$Ca4 + join_groups$Ca5

join_groups$total_free <- join_groups$`Bacteriovore 1`+join_groups$`Bacteriovore 2`+join_groups$`Bacteriovore 3`+join_groups$`Bacteriovore 4`+join_groups$`Fungivore 2`+join_groups$`Fungivore 4`+join_groups$`Omnivore 4`+join_groups$`Omnivore 5`+join_groups$`Carnivore 4`+join_groups$`Carnivore 5`

join_groups$IE <- 100*(join_groups$e/(join_groups$e+join_groups$b))
join_groups$IS <- 100*(join_groups$s/(join_groups$s+join_groups$b))
join_groups$NRC <- 100*((join_groups$`Bacteriovore 1`+join_groups$`Bacteriovore 2`+join_groups$`Bacteriovore 3`+join_groups$`Bacteriovore 4`)/((join_groups$`Bacteriovore 1`+join_groups$`Bacteriovore 2`+join_groups$`Bacteriovore 3`+join_groups$`Bacteriovore 4`)+(join_groups$`Fungivore 2`+join_groups$`Fungivore 4`)))
join_groups$IM <- ((join_groups$`Bacteriovore 1`/join_groups$total_free)*1) + (((join_groups$`Bacteriovore 2`+join_groups$`Fungivore 2`)/join_groups$total_free)*2) + (((join_groups$`Bacteriovore 3`)/join_groups$total_free)*3) + (((join_groups$`Bacteriovore 4`+join_groups$`Fungivore 4`+join_groups$`Omnivore 4`+join_groups$`Carnivore 4`)/join_groups$total_free)*4) + (((join_groups$`Omnivore 5`+join_groups$`Carnivore 5`)/join_groups$total_free)*5)
join_groups$IPP <- ((join_groups$PP2/(join_groups$`Phytoparasite 2`+join_groups$`Phytoparasite 3`+join_groups$`Phytoparasite 5`))*2) + ((join_groups$PP3/(join_groups$`Phytoparasite 2`+join_groups$`Phytoparasite 3`+join_groups$`Phytoparasite 5`))*3) + ((join_groups$PP5/(join_groups$`Phytoparasite 2`+join_groups$`Phytoparasite 3`+join_groups$`Phytoparasite 5`))*5)

fn_vars <- select(join_groups,SN,IE:IPP) 
fn_vars$SN <- as.numeric(fn_vars$SN)
fn_vars <- full_join(site_id,fn_vars) %>% 
  filter(is.na(Input_level)==TRUE) %>% 
  select(-c(Bloc,Input_level))


#### VISUALIZE DATA ####
boxplot(IE ~ SOL, data=fn_vars)
boxplot(IE ~ Treatment_2, data=fn_vars)
boxplot(IS ~ SOL, data=fn_vars)
boxplot(IS ~ Treatment_2, data=fn_vars)
boxplot(NRC ~ SOL, data=fn_vars)
boxplot(NRC ~ Treatment_2, data=fn_vars)
boxplot(IM ~ SOL, data=fn_vars)
boxplot(IM ~ Treatment_2, data=fn_vars)
boxplot(IPP ~ SOL, data=fn_vars)
boxplot(IPP ~ Treatment_2, data=fn_vars)

mean <- aggregate(fn_vars$IE, by = list(fn_vars$SOL, fn_vars$Treatment_2), FUN = mean, na.rm=TRUE)
se <- aggregate(fn_vars$IE, by = list(fn_vars$SOL, fn_vars$Treatment_2), FUN = sd, na.rm=TRUE)
plot <- mean
plot$mean <- mean$x
plot$upper <- mean$x + se$x
plot$lower <- mean$x - se$x
plot <- select(plot, Group.1, Group.2, mean, upper, lower)
names(plot)[1:2] <- c("soil","Group")

ggplot(plot, aes(x = Group, y = mean, fill=soil)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), position=position_dodge(width=0.9), width=0.25) +
  ylab("IE\n") +  
  xlab('') + theme_bw() + 
  scale_fill_manual(values=c('#1b9e77','#d95f02')) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 11)
  )

mean <- aggregate(fn_vars$IS, by = list(fn_vars$SOL, fn_vars$Treatment_2), FUN = mean, na.rm=TRUE)
se <- aggregate(fn_vars$IS, by = list(fn_vars$SOL, fn_vars$Treatment_2), FUN = sd, na.rm=TRUE)
plot <- mean
plot$mean <- mean$x
plot$upper <- mean$x + se$x
plot$lower <- mean$x - se$x
plot <- select(plot, Group.1, Group.2, mean, upper, lower)
names(plot)[1:2] <- c("soil","Group")

ggplot(plot, aes(x = Group, y = mean, fill=soil)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), position=position_dodge(width=0.9), width=0.25) +
  ylab("IS\n") +  
  xlab('') + theme_bw() + 
  scale_fill_manual(values=c('#1b9e77','#d95f02')) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 11)
  )

mean <- aggregate(fn_vars$NRC, by = list(fn_vars$SOL, fn_vars$Treatment_2), FUN = mean, na.rm=TRUE)
se <- aggregate(fn_vars$NRC, by = list(fn_vars$SOL, fn_vars$Treatment_2), FUN = sd, na.rm=TRUE)
plot <- mean
plot$mean <- mean$x
plot$upper <- mean$x + se$x
plot$lower <- mean$x - se$x
plot <- select(plot, Group.1, Group.2, mean, upper, lower)
names(plot)[1:2] <- c("soil","Group")

ggplot(plot, aes(x = Group, y = mean, fill=soil)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), position=position_dodge(width=0.9), width=0.25) +
  ylab("NRC\n") +  
  xlab('') + theme_bw() + 
  scale_fill_manual(values=c('#1b9e77','#d95f02')) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 11)
  )

mean <- aggregate(fn_vars$IM, by = list(fn_vars$SOL, fn_vars$Treatment_2), FUN = mean, na.rm=TRUE)
se <- aggregate(fn_vars$IM, by = list(fn_vars$SOL, fn_vars$Treatment_2), FUN = sd, na.rm=TRUE)
plot <- mean
plot$mean <- mean$x
plot$upper <- mean$x + se$x
plot$lower <- mean$x - se$x
plot <- select(plot, Group.1, Group.2, mean, upper, lower)
names(plot)[1:2] <- c("soil","Group")

ggplot(plot, aes(x = Group, y = mean, fill=soil)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), position=position_dodge(width=0.9), width=0.25) +
  ylab("IM\n") +  
  xlab('') + theme_bw() + 
  scale_fill_manual(values=c('#1b9e77','#d95f02')) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 11)
  )

mean <- aggregate(fn_vars$IPP, by = list(fn_vars$SOL, fn_vars$Treatment_2), FUN = mean, na.rm=TRUE)
se <- aggregate(fn_vars$IPP, by = list(fn_vars$SOL, fn_vars$Treatment_2), FUN = sd, na.rm=TRUE)
plot <- mean
plot$mean <- mean$x
plot$upper <- mean$x + se$x
plot$lower <- mean$x - se$x
plot <- select(plot, Group.1, Group.2, mean, upper, lower)
names(plot)[1:2] <- c("soil","Group")

ggplot(plot, aes(x = Group, y = mean, fill=soil)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), position=position_dodge(width=0.9), width=0.25) +
  ylab("IPP\n") +  
  xlab('') + theme_bw() + 
  scale_fill_manual(values=c('#1b9e77','#d95f02')) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 11)
  )


