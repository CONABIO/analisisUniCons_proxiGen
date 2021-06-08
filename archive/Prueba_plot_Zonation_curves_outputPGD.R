### This script plots the performance curves and proportion of
### proxies of genetid diversity conserved under Zonation
### analysis considering the SDM, PDG, habitat and threat layers


### Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

###### Read data

## Read conservation features names
zfeatures<-read.table("../data/spatial/Zonation_output/04_MDP_vs_PDG.features_info.txt",
                      skip=1, header = TRUE)

# extract conservation features names (each SDP subdivided by its PDG)
feature_names<-substring(zfeatures$MapFileName, first= 9) # cut MDP_PDG/
feature_names<-gsub(pattern=".tif", replacement= "", feature_names) # cut .tif

# check names
feature_names[1]
feature_names[length(feature_names)]

## Read IUCN category per taxa
iucn<-read.csv("../data/spatial/Zonation_7/IUCN_threat_category.csv")


## Read Zonation output 
zcurves<-read.table("../data/spatial/Zonation_output/04_MDP_vs_PDG.curves.txt",
                    header=FALSE,
                    comment.char = "#") # firt line is comments, not real col names

###### Tidy Zonation output

## add col names manually
colnames(zcurves) <- c("Prop_landscape_lost", "cost_needed_for_top_fraction",
                       "min_prop_rem",  "ave_prop_rem", "W_prop_rem","ext_1", "ext_2", 
                       feature_names) # from column 8 onwards column names correspond to feature names

# Check few entries and cols of data       
eg<-zcurves[1:5, 1:10]
eg

### Tidy data 
# change to have one column per taxon
zcurves_tidy<-gather(zcurves, feature, prop, 8:ncol(zcurves)) %>%
  separate(feature, c("Nombre_", "PGD"), sep= "_zv_", remove=FALSE) %>%
  
  # add Area kept 
  mutate(Area_kept=(1-Prop_landscape_lost)*100) %>%
  
  # add iucn data
  left_join(iucn, by ="Nombre_") %>%
  
  # filter species that were not evaluated in iucn workshop
  filter(IUCN.threat.category!="NO EVALUADA") %>% droplevels(.) 

head(zcurves_tidy)



##### Prepare data for plot and plot

## Get mean_prop for all features

# sum of prop area for each taxa within each area_kept
all_areakept<-group_by(zcurves_tidy, Nombre_, Area_kept)  %>%
  summarise(mean_prop=mean(prop))
head(all_areakept)
tail(all_areakept)

all_areakept[100,]

# check
ggplot(all_areakept, aes(x=Area_kept, y=mean_prop, color=Nombre_)) +
  geom_point() + theme(legend.position = "none")



## Get mean_prop by IUCN category

# check iucn categories
levels(zcurves_tidy$IUCN.threat.category)

# make new categories with CR and EN in a single one
levels(zcurves_tidy$IUCN.threat.category)<-c("CR+EN","DD","CR+EN","LC","NT","VU")

## Estimate area kept for category
# CR + EN
CR_EN<-filter(zcurves_tidy, IUCN.threat.category=="CR+EN") %>%
  group_by(Area_kept) %>%
  summarise(CR_EN=mean(prop))

# DD
DD<-filter(zcurves_tidy, IUCN.threat.category=="DD") %>%
  group_by(Area_kept) %>%
  summarise(DD=mean(prop))

# LC
LC<-filter(zcurves_tidy, IUCN.threat.category=="LC") %>%
  group_by(Area_kept) %>%
  summarise(LC=mean(prop))

# NT
NT<-filter(zcurves_tidy, IUCN.threat.category=="NT") %>%
  group_by(Area_kept) %>%
  summarise(NT=mean(prop))

# VU
VU<-filter(zcurves_tidy, IUCN.threat.category=="VU") %>%
  group_by(Area_kept) %>%
  summarise(VU=mean(prop))

## Join all together for plot
zc_forplot<-list(all_areakept, CR_EN, DD, LC, NT,VU) %>%
  reduce(left_join, by="Area_kept") %>%
  gather(Category, Mean_prop, -Area_kept)
head(zc_forplot)

## Plot
ggplot(zc_forplot, aes(x=Area_kept, y=Mean_prop, color=Category)) +
  geom_point() 



