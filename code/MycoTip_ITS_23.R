# Project: MycoTip year 1 (2023) ----
# Analysis of root tip HTS data.
# Author and date: Kelsey Patrick 24.09.02
# Summary: Root tips were sampled from Carya ovata and Pinus abies during 
# six time-points from April-December of 2023. DNA was extracted and sequenced
# in the spring of 2024. ASVs were quality controlled and rarified by Peter Kennedy.
library(dplyr)
library(vegan)
library(tidyverse)
library(stringr)
library(pals)

# Tidy data for initial analysis ----
# use vroom::vroom to read csv table into a tibble 
ASVtbl <- vroom::vroom("data/processed_data/ITS_data6.fungi.seq.tax.rar.csv")

# 1510 TOTAL ASVs (species)
# use glimpse function (returns columns vertically and data horizontally), 
# str returns object structure, and dim returns dimensions as rowsXcolumns
# glimpse(ASVtbl); str(ASVtbl); dim(ASVtbl)

# produce pivot table using pivot_longer function and selecting rows 2:69
# rows 2:69 are selected because they contain the sample iDs and associated metadata
# we now have abundance of each ASV in each sample 
ASVtbl <- ASVtbl %>% pivot_longer(2:69, names_to = 'sample' , values_to = 'abundance')

# filter out any ASVs with 0 abundance, 180,404 -> 4,399
ASVtbl <- ASVtbl %>% filter(abundance > 0)

# select relevant columns and produce a new tibble
# extract session, plot, and subplot from 'sample'
mod_ASV <- ASVtbl %>% select(asv_id, sample, abundance, Ectomycorrhiza_lineage_template, 
                             Family, Genus, Species, primary_lifestyle,
                             Ectomycorrhiza_exploration_type_template)
mod_ASV$session <- substr(ASVtbl$sample, 1,1)
mod_ASV$plot <- substr(ASVtbl$sample, 2,3)
mod_ASV$subplot <- substr(ASVtbl$sample, 4,4)
# also add date and change plot abbreviations to species names
mod_ASV$date <- NA
mod_ASV$date[mod_ASV$session==4] <- '23.6.29'
mod_ASV$date[mod_ASV$session==5] <- '23.7.20'
mod_ASV$date[mod_ASV$session==6] <- '23.8.23'
mod_ASV$date[mod_ASV$session==7] <- '23.9.28'
mod_ASV$date[mod_ASV$session==8] <- '23.11.1'
mod_ASV$date[mod_ASV$session==9] <- '23.12.6'
mod_ASV$plot[mod_ASV$plot=='CA'] <- 'Carya_ovata'
mod_ASV$plot[mod_ASV$plot=='PI'] <- 'Picea_abies'

# reorder columns by column index
mod_ASV <- mod_ASV[, c(2,1,3,10,13,11,12,4,5,6,7,8,9)]

# change chr to factors
mod_ASV$asv_id <- as.factor(mod_ASV$asv_id)
mod_ASV$session <- as.factor(mod_ASV$session)
mod_ASV$plot <- as.factor(mod_ASV$plot)
mod_ASV$subplot <- as.factor(mod_ASV$subplot)
mod_ASV$Genus <- as.factor(mod_ASV$Genus)
mod_ASV$primary_lifestyle <- as.factor(mod_ASV$primary_lifestyle)
mod_ASV$Ectomycorrhiza_exploration_type_template <- as.factor(mod_ASV$Ectomycorrhiza_exploration_type_template)

# shift Hyaloscypha to ectomycorrhizal primary lifestyle and Cladophialophore
# and Phialocephala to endophytes ... delete old rows and join new
Hyal <- mod_ASV %>% filter(Genus == 'Hyaloscypha') %>%
  mutate(primary_lifestyle = recode(primary_lifestyle, 'litter_saprotroph' = 'ectomycorrhizal')) 
Phia <- mod_ASV %>% filter(Genus == 'Phialocephala') %>%
  mutate(primary_lifestyle = recode(primary_lifestyle, 'soil_saprotroph' = 'root_endophyte')) 
Clad <- mod_ASV %>% filter(Genus == 'Cladophialophora') %>%
  mutate(primary_lifestyle = recode(primary_lifestyle, 'soil_saprotroph' = 'root_endophyte'))
mod_ASV <- mod_ASV[!(mod_ASV$Genus %in% c('Hyaloscypha', 'Phialocephala', 'Cladophialophora')),]
mod_ASV <- rbind(mod_ASV, Hyal, Phia, Clad)
rm(Hyal, Phia, Clad, ASVtbl)

# Lifestyle type analysis ----

# get occurrences within primary_lifestyle column 
# mod_ASV$primary_lifestyle %>% table()
# shows 19 lifestyle types, EMF dominate at 34,816, litter and soil SAPs both around 17k

# group by primary lifestyle
lifestyle_sum <- mod_ASV %>%
  group_by(primary_lifestyle) %>%
  summarise(occurrence = n(), 'abundance_sum' = sum(abundance)) %>%
  arrange(desc(abundance_sum)) %>% 
  ungroup() %>% 
  mutate('abundance_sum_prop' = abundance_sum / sum(abundance_sum))

# plot total lifestyle type occurrence
ggplot(lifestyle_sum, aes(fct_rev(fct_reorder(primary_lifestyle, occurrence)), occurrence))+ 
  geom_col()+
  coord_flip()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = 'Occurrence of lifestyle types',
       x = element_blank(), y = 'sequence occurrence')

# plot total lifestyle type abundance
ggplot(lifestyle_sum, aes(fct_rev(fct_reorder(primary_lifestyle, occurrence)), abundance_sum))+ 
  geom_col()+
  coord_flip()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = 'Abundance of lifestyle types',
       x = element_blank(), y = 'sequence abundance')

## Main guild mean, median, and proportion, over time ----

# fix below: filter out desired guilds first, combine sap, the group by to 
# get mean, median, and sd (look at EMF analysis)

# group by lifestyle, session, and plot to get abundance and proportion
lifestyle_sum <- mod_ASV %>%
  group_by(primary_lifestyle, session, plot) %>%
  summarise(occurrence = n(), 'abundance_sum' = sum(abundance)) %>%
  arrange(desc(abundance_sum)) %>% 
  ungroup(primary_lifestyle) %>% 
  mutate('abundance_sum_prop' = abundance_sum / sum(abundance_sum)) %>% 
  ungroup()

# combine SAP types into one category
lifestyle_sum$primary_lifestyle <- fct_collapse(lifestyle_sum$primary_lifestyle, 
                                                combined_sap = 
                                                  c('unspecified_saprotroph', 'nectar/tap_saprotroph',
                                                    'soil_saprotroph', 'wood_saprotroph', 'dung_saprotroph',
                                                    'litter_saprotroph'))

# get sum of all 'combined_sap' rows per each session and plot
lifestyle_sum <- lifestyle_sum %>%
  group_by(primary_lifestyle, session, plot) %>%
  summarise(n_sap = n(), 'abundance_sum' = sum(abundance_sum)) %>%
  arrange(desc(abundance_sum)) %>% 
  ungroup(primary_lifestyle) %>% 
  mutate('abundance_sum_prop' = abundance_sum / sum(abundance_sum)) %>% 
  ungroup()

# subset dataframe to only include EMF, SAP, PATH, and root endos
lifestyle_sum <- lifestyle_sum[lifestyle_sum$primary_lifestyle %in% 
                                 c('ectomycorrhizal', 'combined_sap',
                                   'root_endophyte', 'plant_pathogen'), ]

# convert "session" from chr to numeric
lifestyle_sum$session <- as.numeric(as.character(lifestyle_sum$session))

# plot the abundance of lifestyle types over time
ggplot(lifestyle_sum, aes(session, abundance_sum, color = primary_lifestyle))+
  geom_point()+
  geom_line(aes(linetype = primary_lifestyle))+
  facet_grid(~ plot)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(colour="white", fill="white"))+
  scale_x_continuous(labels=c('June', 'July', 'August', 'September', 'November', 'December'))+
  labs(title = 'Abundance of lifestyle types across time',
       x = element_blank(), y = 'sequence abundance')+
  scale_colour_grey()







## Main guild top species (SAP,ENDO,PATH) in each plot ----
# EMF analysis ----

# filter to include only ectomycorrhizal lifestyle 
EMF <- mod_ASV %>% filter(primary_lifestyle == 'ectomycorrhizal')

# group by Genus to get abundance, occurrence, and proportion (overall)
EMF_gen <- EMF %>%
  group_by(Genus) %>%
  summarise(occurrence = n(), 'abundance_sum' = sum(abundance)) %>%
  arrange(desc(abundance_sum)) %>% 
  ungroup() %>% 
  mutate('abundance_sum_prop' = abundance_sum / sum(abundance_sum))

# plot total genera occurrence
ggplot(EMF_gen, aes(fct_rev(fct_reorder(Genus, occurrence)), occurrence))+ 
  geom_col()+
  coord_flip()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = 'Occurrence of genera',
       x = element_blank(), y = 'sequence occurrence')

# plot total genera abundance
ggplot(EMF_gen, aes(fct_rev(fct_reorder(Genus, occurrence)), abundance_sum))+ 
  geom_col()+
  coord_flip()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = 'Abundance of genera',
       x = element_blank(), y = 'sequence abundance')

# subset most abundant general (top 10)
EMF_gen <- EMF_gen[EMF_gen$Genus %in% 
                     c('Amphinema', 'Tomentella', 'Tuber', 'Cenococcum',
                       'Russula', 'Tylospora', 'Hygrophorus', 'Inocybe', 
                       'Wilcoxina', 'Sebacina'),]

# plot total abundance by genus (top 10)
ggplot(EMF_gen, aes(fct_rev(fct_reorder(Genus, abundance_sum)), abundance_sum))+ 
  geom_col()+
  coord_flip()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = 'Abundance of top 10 genera',
       x = element_blank(), y = 'sequence abundance')

## Top 10 mean, median, and proportion over time ----

# filter to include only ectomycorrhizal lifestyle 
EMF <- mod_ASV %>% filter(primary_lifestyle == 'ectomycorrhizal')

# group by genus, session, and plot to get session sum of abundance and proportion
EMF_avg <- EMF %>%
  group_by(Genus, session, plot) %>%
  summarise('session_abun' = sum(abundance), 'session_median' = median(abundance),
            'sd' = sd(abundance)) %>%
  arrange(desc(session_abun)) %>%
  ungroup(Genus) %>% 
  mutate('session_prop' = session_abun / sum(session_abun)) %>% 
  ungroup()

# check that median and sd calculations are correct (SD is between samples, not population)
# checkSD <- EMF%>% filter(Genus == 'Elaphomyces', plot == 'Carya_ovata', session == 4)

# use summarise to see how many subplots per session
subplot_occurrence <- EMF %>%
  group_by(Genus, session, plot, subplot) %>%
  summarise('session_abun' = sum(abundance)) %>%
  arrange(desc(session_abun)) %>% 
  ungroup(subplot) %>% 
  summarise(n = n()) %>% 
  ungroup() 

# only include top 10 genera
EMF_avg <- EMF_avg[EMF_avg$Genus %in% 
                     c('Amphinema', 'Tomentella', 'Tuber', 'Cenococcum',
                       'Russula', 'Tylospora', 'Hygrophorus', 'Inocybe', 
                       'Wilcoxina', 'Sebacina'),]

subplot_occurrence <- subplot_occurrence[subplot_occurrence$Genus %in% 
                                           c('Amphinema', 'Tomentella', 'Tuber', 'Cenococcum',
                                             'Russula', 'Tylospora', 'Hygrophorus', 'Inocybe', 
                                             'Wilcoxina', 'Sebacina'),]

# merge avg and subplot data frames
EMF_avg <- merge(EMF_avg, subplot_occurrence,
                 by=c('Genus', 'plot','session'))

# mutate to get average abundance per session and plot
EMF_avg <- mutate(EMF_avg, 'avg' = session_abun/n)

# convert "session" from chr to number
EMF_avg$session <- as.numeric(as.character(EMF_avg$session))

# plot mean abundance per plot over time (RETURN HERE TO ADD SD)
ggplot(EMF_avg, aes(session, avg, color = Genus))+
  geom_point()+
  geom_line()+
  facet_grid(~ plot)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(colour="white", fill="white"))+
  scale_x_continuous(labels=c('June', 'July', 'August', 'September', 'November', 'December'))+
  labs(title = 'Mean abundance of genera across time',
       x = element_blank(), y = 'Mean sequence abundance')+
  scale_colour_manual(values=unname(cols25()))

## Exploration types ----

# occurrences of exploration types
# returns short-distance-delicate > medium-distance-smooth > short-distance-coarse
EMF <- mod_ASV %>% filter(primary_lifestyle == 'ectomycorrhizal')
EMF$Ectomycorrhiza_exploration_type_template %>% table()

# sum of ET abundance, plus occurrences and proportions 
ET_sum <- EMF %>%
  group_by(Ectomycorrhiza_exploration_type_template) %>%
  summarise(occurrence = n(), 'abundance_sum' = sum(abundance)) %>%
  arrange(desc(abundance_sum)) %>% 
  ungroup() %>% 
  mutate('abundance_sum_prop' = abundance_sum / sum(abundance_sum))

# plot total ET occurrence
ggplot(ET_sum, aes(fct_rev(fct_reorder(Ectomycorrhiza_exploration_type_template,
                                       occurrence)), occurrence))+ 
  geom_col()+
  coord_flip()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = 'Occurrence of exploration types',
       x = element_blank(), y = 'sequence occurrence')

# plot total ET abundance
ggplot(ET_sum, aes(fct_rev(fct_reorder(Ectomycorrhiza_exploration_type_template,
                                       occurrence)), abundance_sum))+ 
  geom_col()+
  coord_flip()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = 'Abundance of exploration types',
       x = element_blank(), y = 'sequence abundance')

# sum of ET abundance by SESSION
ET_sum <- EMF %>%
  group_by(Ectomycorrhiza_exploration_type_template, session) %>%
  summarise(occurrence = n(), 'abundance_sum' = sum(abundance)) %>%
  arrange(desc(abundance_sum)) %>% 
  ungroup(Ectomycorrhiza_exploration_type_template) %>% 
  mutate('abundance_sum_prop' = abundance_sum / sum(abundance_sum)) %>% 
  ungroup()

# convert "session" from chr to number
ET_sum$session <- as.numeric(as.character(ET_sum$session))

ggplot(ET_sum, aes(session, abundance_sum, color = Ectomycorrhiza_exploration_type_template))+
  geom_point()+
  geom_line()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(labels=c('June', 'July', 'August', 'September', 'November', 'December'))+
  labs(title = 'Abundance of ETs across time',
       x = element_blank(), y = 'sequence abundance',
       color = 'Exploration type')+
  scale_colour_manual(values=unname(cols25()))

# plot ET abundance by SESSION and PLOT 
ET_sum <- EMF %>%
  group_by(Ectomycorrhiza_exploration_type_template, session, plot) %>%
  summarise(occurrence = n(), 'abundance_sum' = sum(abundance)) %>%
  arrange(desc(abundance_sum)) %>% 
  ungroup(Ectomycorrhiza_exploration_type_template) %>% 
  mutate('abundance_sum_prop' = abundance_sum / sum(abundance_sum)) %>% 
  ungroup()

# convert "session" from chr to number
ET_sum$session <- as.numeric(as.character(ET_sum$session))

ggplot(ET_sum, aes(session, abundance_sum, color = Ectomycorrhiza_exploration_type_template))+
  geom_point()+
  geom_line()+
  facet_grid(~plot)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(colour="white", fill="white"))+
  scale_x_continuous(labels=c('June', 'July', 'August', 'September', 'November', 'December'))+
  labs(title = 'Abundance of ETs across time',
       x = element_blank(), y = 'sequence abundance',
       color = 'Exploration type')+
  scale_colour_manual(values=unname(cols25()))

### Mean ET abundance ----
ET_mean <- EMF %>%
  group_by(Ectomycorrhiza_exploration_type_template, session, plot) %>%
  summarise(occurrence = n(), 'session_abun' = sum(abundance)) %>%
  arrange(desc(session_abun)) %>% 
  ungroup()

# use summarise to see how many subplots per session
subplot_ETcount <- EMF %>%
  group_by(Ectomycorrhiza_exploration_type_template, session, plot, subplot) %>%
  summarise('session_abun' = sum(abundance)) %>%
  arrange(desc(session_abun)) %>% 
  ungroup(subplot) %>% 
  summarise(n = n()) %>% 
  ungroup() 

# merge gen_mean and subplot data frames
ET_mean <- merge(ET_mean, subplot_ETcount,
                 by=c('Ectomycorrhiza_exploration_type_template', 'plot','session'))

# mutate to get average genus abundance per session and plot
ET_mean <- mutate(ET_mean, 'avg' = session_abun/n)

# convert "session" from chr to number
ET_mean$session <- as.numeric(as.character(ET_mean$session))

ggplot(ET_mean, aes(session, avg, color = Ectomycorrhiza_exploration_type_template))+
  geom_point()+
  geom_line()+
  facet_grid(~ plot)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(colour="white", fill="white"))+
  scale_x_continuous(labels=c('June', 'July', 'August', 'September', 'November', 'December'))+
  labs(title = 'Mean abundance of ETs across time',
       x = element_blank(), y = 'Mean sequence abundance',
       color = 'Exploration type')+
  scale_colour_manual(values=unname(cols25()))

## Diversity analysis ----
### Shannon's index ----