# set wd to source file location
# Load needed libraries

library(sp)
library(raster)
library(dplyr)
library(ggplot2)
library(tidyr)
library(vegan)
library(scales)

### Load Zonation outputs and proxies of gen div

# Path to rasters
myrasters <-c("../data/spatial/Zonation_output/comparacion/E_final_Todos.rank.compressed.tif", # Zonation output with species only (n=116) #1
             "../data/spatial/Zonation_output/comparacion/E_final_VegPyS.rank.compressed.tif", # Zonation output with species + ZV (n=143) #2
             "../data/spatial/Zonation_output/comparacion/E_final_HabVarios.rank.compressed.tif", # Zonation output with species and PDG (n=218) #3
                         "../data/spatial/areasProxyDivGen/PDG.tif") #Proxies div gen #6

# wt NOTE: IF NEEDED, CHANGE LAYER NAMES AND ADJUST Number (e.g. 2:8 -> 2:7; check lines 71, 260, 265)


# Stack rasters and add nicer names 
CWR<-stack(myrasters)
names(CWR)<-c("Zonation_1", "Zonation_2", "Zonation_3", "ProxiesDivGen")

#### Data description #####

ZT=0.1 #Zonation threshold (area will be calculated as (Area=1-ZT))

### plots
#colors for proxies div gen 
mycols<-c("grey", "#668325", "#9343c1", "#61b928", "#bc6cec", "#56cd5e", "#4667e6", "#bab71b", "#704dc1", "#b5cb3e", "#4c4fb8", "#83b124", "#5e7ef2", "#38a435", "#cc3ca6", "#3fcb85", "#ec49a6", "#86ca62", "#e471de", "#5d9d38", "#9b76e8", "#9ea122", "#a13ca0", "#3da35c", "#e53283", "#2b7c2c", "#b063c4", "#c7bb43", "#3b65c1", "#e2b438", "#628ceb", "#e69128", "#509be7", "#e0422f", "#34c7dd", "#c62a33", "#57d2b9", "#d1204a", "#57c597", "#ef3d76", "#56cdc9", "#ec4b5f", "#2ba198", "#e2672a", "#2baad3", "#b33f1d", "#6bbff0", "#ae6314", "#7f68c4", "#90ac42", "#774598", "#c0bd59", "#b58be8", "#4e6313", "#dd80d3", "#87be71", "#b62a7c", "#3b9065", "#c22a60", "#71b588", "#e75687", "#326a38", "#e46aad", "#628848", "#ad498c", "#b2c06d", "#535a99", "#c89d3c", "#3176ae", "#e97150", "#549bce", "#ca783f", "#7b7ec0", "#927b1f", "#a6a5e8", "#9a9c49", "#a161a4", "#a6b174", "#7d4f8c", "#e99b5c", "#206e54", "#b1394e", "#505e25", "#dc9ce0", "#786e30", "#ba7eb8", "#90652d", "#de93ba", "#caa569", "#9b3a66", "#e99d7a", "#8c4e6c", "#914c27", "#f288aa", "#b7443b", "#c6648b", "#bc7257", "#dd6e81", "#a7514b", "#e49190", "#a04f59", "#e6736d")


plot(CWR$ProxiesDivGen, main= "Areas Proxi of Genetic Diversity", col=mycols, legend=FALSE)
plot(CWR$Zonation_1>=(1-ZT), legend=FALSE, main= "All taxa")
plot(CWR$Zonation_2>=(1-ZT), legend=FALSE, main= "Taxa asociated to vegetation")
plot(CWR$Zonation_3>=(1-ZT), legend=FALSE, main= "Taxa asociated to various habitat")
#plot(CWR$Zonation_MDP_vs_PDG>=(1-ZT), legend=FALSE, main= "Zonation output SDMvsProxies")
#plot(CWR$Zonation_MDP_PDG_ADMU>=(1-ZT), legend=FALSE, main= "Zonation output SDM and Proxies as ADMU")

### Describe proxies div gen
table(getValues(CWR$ProxiesDivGen))



##### Analyses per sp #####

# Create where to store results for all spp
sol_tidy_spp<-character(0) # area of proxies of div gen for each zonation solution
sol_prop_spp<-character(0) # proportion of proxi of div gen for each zolnation solution in relation to its area in the sp distribution
sols_summary_spp<-character(0) # diversity indexes and mean of proportion 

# Remember zonation threshold used:
ZT

# List of spp
#correr hasta la l?nea 317
spp<-list.files("../data/spatial/modelos_darwin_all_final/")
spp


for (i in 1:length(spp)){ 

spname<-gsub(".tif", "", spp[i])
spname<-gsub("_", " ", spname)
print(spname)
  
# add to stacks of raster the raster with the sp distribution model

CWR_sp<-stack(CWR, paste0("../data/spatial/modelos_darwin_all/", spp[i])) #  sp model
names(CWR_sp)[5]<-"ENM_sp"

# plot
plot(CWR_sp$ENM_sp, legend=FALSE, main= paste(spname, "distribution"))

# Define areas with Proxies of DivGen for that sp
# this works because the ENM_sp is coded 1= sp presence 0= absence
sp_ProxiesDivGen<- CWR_sp$ProxiesDivGen * CWR_sp$ENM_sp

# write to raster
writeRaster(sp_ProxiesDivGen, format= "GTiff", file= paste0("../data/spatial/areasProxyDivGen/crop_to_sp/", "PDG_", spp[i]), overwrite=TRUE)

############## Quantify overlaps

### convert to rasters to points to have a "df of overlaps"
CWR_overlaps<-as.data.frame(rasterToPoints(CWR_sp))

##### calcs for each solution

##### 1) Proxies of Div Gen in base ENM of sp
## subset
x<-dplyr::filter(CWR_overlaps, ENM_sp==1) %>% arrange(., ProxiesDivGen)
x$ProxiesDivGen<-as.character(x$ProxiesDivGen)

# save results in df
sol_base<-as.data.frame(table(x$ProxiesDivGen))
names(sol_base)<-c("ProxiDivGen", "ENM_sp")


### plots

## bar plot 
#ggplot(x, aes(ProxiesDivGen)) + geom_bar() + ggtitle(paste("Areas Proxi of Genetic Diversity for", spp[i]))


##### 2) CWRsp + Zonation_SPP

### calcs
## subet
x<-dplyr::filter(CWR_overlaps, ENM_sp==1, Zonation_1>=(1-ZT)) %>% arrange(., ProxiesDivGen)
x$ProxiesDivGen<-as.character(x$ProxiesDivGen)

# Save to common data frame
sol_Z_1<-as.data.frame(table(x$ProxiesDivGen))
names(sol_Z_1)<-c("ProxiDivGen", "Zonation_1")

## bar plot 
# ggplot(x, aes(ProxiesDivGen)) + geom_bar() + ggtitle(paste(spname,"conserved under ZonationSPP at", ZT))


##### 3) CWRsp + Zonation SPPyPDG

### calcs
## subet
x<-dplyr::filter(CWR_overlaps, ENM_sp==1, Zonation_2>=(1-ZT)) %>% arrange(., ProxiesDivGen)
x$ProxiesDivGen<-as.character(x$ProxiesDivGen)


# Save to common data frame
sol_Z_2<-as.data.frame(table(x$ProxiesDivGen))
names(sol_Z_2)<-c("ProxiDivGen", "Zonation_2")


## bar plot 
# ggplot(x, aes(ProxiesDivGen)) + geom_bar() + ggtitle(paste(spname, "conserved under ZonationZV at", ZT))


##### 4) CWRsp + Zonation SPPvsPDG publicado en sintesis

### calcs
## subet
x<-dplyr::filter(CWR_overlaps, ENM_sp==1, Zonation_3>=(1-ZT)) %>% arrange(., ProxiesDivGen)
x$ProxiesDivGen<-as.character(x$ProxiesDivGen)

# Save to common data frame
sol_Z_3<-as.data.frame(table(x$ProxiesDivGen))
names(sol_Z_3)<-c("ProxiDivGen", "Zonation_3")


## bar plot 
# ggplot(x, aes(ProxiesDivGen)) + geom_bar() + ggtitle(paste(spname, "conserved under ZonationProxiDivGen at", ZT))





####### Analyse results together


### ALL
## Prepare data
# Join all solutions
sol_all<-left_join(sol_base, sol_Z_1) %>% left_join(., sol_Z_2) %>% left_join(., sol_Z_3) %>%  
  mutate_all(funs(replace(., is.na(.), 0))) # NAs are actually 0s
head(sol_all)

# gather into tidy df
sol_tidy<-gather(sol_all, Solution, Freq, -ProxiDivGen)
head(sol_tidy)

## Plot
#print(ggplot(sol_tidy, aes(ProxiDivGen, Freq)) + geom_col(aes(fill=Solution)) + ggtitle(paste(spname, "at Zonation threshold", ZT)))

print(ggplot(sol_tidy, aes(x=ProxiDivGen, y= Freq, fill=Solution)) + geom_col(position="dodge") + ggtitle(paste(spname, "at Zonation threshold", ZT)))

#print(ggplot(sol_tidy, aes(x=ProxiDivGen, y= Freq, fill=Solution)) + geom_col(position="fill") + ggtitle(paste(spname, "at Zonation threshold", ZT)))


# save results with the rest of the spp
sol_tidy$sp<-spname
sol_tidy_spp<-rbind(sol_tidy_spp, sol_tidy)


## Estimate Proportion of Zones Proxies Div Gen under conservation for each Zonation solution

# estimate proportion
sol_prop<-mutate(sol_all,
                 Zonation_1=Zonation_1/ENM_sp, 
                 Zonation_2=Zonation_2/ENM_sp,
                 Zonation_3=Zonation_3/ENM_sp)
head(sol_prop)

# tidy
sol_prop<-gather(sol_prop, Solution, Proportion_ProxiDG, -ProxiDivGen, -ENM_sp)
tail(sol_prop)

    
# save results with the rest of the spp
sol_prop$sp<-spname
sol_prop<-as.data.frame(sol_prop)
sol_prop_spp<-rbind(sol_prop_spp, sol_prop)

## make summary of proportion and area

#prop
myprop<-group_by(sol_prop, Solution) %>% 
  summarise(., mean.prop = mean(Proportion_ProxiDG), median.prop=median(Proportion_ProxiDG))

#area
tot.area.sp<-sum(sol_all$ENM_sp)

myarea<-sol_tidy %>% 
  group_by(., Solution) %>% summarise(., Area=sum(Freq)) %>%
  mutate(., Prop_to_AreaSP=Area/tot.area.sp)

sol_prop_summary<-left_join(myarea, myprop, by="Solution")       
sol_prop_summary<-as.data.frame(sol_prop_summary)

### Estimate Simpson diversity index

## we need an abundance matrix with species (ProxisDivGen in this case) in columns
# and abundance by treatment (Zonation output in this case) in rows. Therefore:

# transpose data to produce abundace matrix as expected by vegan
row.names(sol_all)<-sol_all[,1]
abundance.matrix<-t(sol_all[, 2:5])
class(abundance.matrix)

### Estimate Diversity Indices

indices.div<-data.frame("Solution"=colnames(sol_all[2:5]))

# Richness
indices.div$Richness<-rowSums(abundance.matrix>0)

# Shannon
indices.div$shannon<-diversity(abundance.matrix, index="shannon") # shannon is default

# Simpson
indices.div$simpson<-diversity(abundance.matrix, index="simpson")


#### Distance among conserving the whole distribution of the sp ("optimun scenario possible") and Zonation solutions

sol_dist<-data.frame(Solution=c("ENM_sp", "Zonation_SPP", "Zonation_SPPyPDG", "Zonation_sintesis" ,"Zonation_new" ,"Zonation_newBLP"),
                     Dist.to.Optimun=c(as.vector(dist(rbind(sol_all$ENM_sp, sol_all$ENM_sp))),
                                       as.vector(dist(rbind(sol_all$ENM_sp, sol_all$Zonation_1))),
                                       as.vector(dist(rbind(sol_all$ENM_sp,sol_all$ Zonation_2))),
                                       as.vector(dist(rbind(sol_all$ENM_sp, sol_all$Zonation_3)))))


### join summary results together 

sols_summary<-left_join(indices.div, sol_prop_summary, "Solution") %>% left_join(., sol_dist, "Solution")
print(sols_summary)

###  save with the rest of the spp
sols_summary$sp<-spname
sols_summary_spp<-rbind(sols_summary_spp, sols_summary)

}

##### Analyse results together

head(sol_tidy_spp)
head(sol_prop_spp)
head(sols_summary_spp)

# write results 
###EL OUTPUT DEBE GAURDARSE EN 3_final_scenarios
write.csv(sol_tidy_spp, "../data/comparations_output/area/10p_sol_tidy_spp.txt", row.names=FALSE, quote = FALSE)
write.csv(sol_prop_spp, "../data/comparations_output/area/10p_sol_prop_spp.txt", row.names=FALSE, quote = FALSE)
write.csv(sols_summary_spp, "../data/comparations_output/area/10p_sols_summary_spp.txt", row.names=FALSE, quote = FALSE)

### Plots

## Violin plots

# dist to optimun

filter(sols_summary_spp, Solution!="ENM_sp") %>%
       ggplot(., aes(x=Solution, y=Dist.to.Optimun, color=Solution)) + geom_violin() + geom_jitter(size = 0.3) +
       stat_summary(fun.y=mean, geom="point", color="red") +
       theme(axis.text.x = element_blank()) +
       scale_y_continuous(name="Distance to Optimun")

# Simpson

filter(sols_summary_spp, Solution!="ENM_sp") %>%
  ggplot(., aes(x=Solution, y=simpson, color=Solution)) + geom_violin() + geom_jitter(size = 0.3) +
  stat_summary(fun.y=mean, geom="point", color="red") +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(name="Simpson Index")

# Area sp conserved

filter(sols_summary_spp, Solution!="ENM_sp") %>%
  ggplot(., aes(x=Solution, y=Prop_to_AreaSP, color=Solution)) + geom_violin() + geom_jitter(size = 0.3) +
  stat_summary(fun.y=mean, geom="point", color="red") +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(name="Proportion of distribution of the sp conserved")

# Mean prop area proxies by sp ponserved

filter(sols_summary_spp, Solution!="ENM_sp") %>%
  ggplot(., aes(x=Solution, y=mean.prop, color=Solution)) + geom_violin() + geom_jitter(size = 0.3) +
  stat_summary(fun.y=mean, geom="point", color="red") +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(name="Mean proportion of the area of each proxy conserved by sp")


# Median prop area proxies by sp ponserved

filter(sols_summary_spp, Solution!="ENM_sp") %>%
  ggplot(., aes(x=Solution, y=median.prop, color=Solution)) + geom_violin() + geom_jitter(size = 0.3) +
  stat_summary(fun.y=mean, geom="point", color="red") +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(name="Median proportion of the area of each proxy conserved by sp")

# Median and mean together

filter(sols_summary_spp, Solution!="ENM_sp") %>%
  gather(., SummaryStat, value, -Solution, -c(Richness:Prop_to_AreaSP), -Dist.to.Optimun, -sp) %>%
  ggplot(., aes(x=interaction(SummaryStat, Solution), y=value, color=SummaryStat)) + geom_violin() + geom_jitter(size = 0.1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
         
 
# Proportion of the distribution of the sp conserved vs MEAN proportion of each proxy conserved by sp

filter(sols_summary_spp, Solution!="ENM_sp") %>%
  ggplot(., aes(x=Prop_to_AreaSP, y=mean.prop, color=Solution)) + geom_point() +
  scale_x_continuous(name="Proportion of distribution of the sp conserved") +
  scale_y_continuous(name="Mean proportion of the area of each proxy conserved by sp")


# Proportion of the distribution of the sp conserved vs MEADIAN proportion of each proxy conserved by sp

filter(sols_summary_spp, Solution!="ENM_sp") %>%
  ggplot(., aes(x=Prop_to_AreaSP, y=median.prop, color=Solution)) + geom_point() +
  scale_x_continuous(name="Proportion of distribution of the sp conserved") +
  scale_y_continuous(name="Median proportion of the area of each proxy conserved by sp")


## with regression lines

filter(sols_summary_spp, Solution!="ENM_sp") %>%
  ggplot(., aes(x=Prop_to_AreaSP, y=mean.prop, color=Solution)) + geom_point(size=0.7) +
  scale_x_continuous(name="Proportion of species distribution (%)") +
  scale_y_continuous(name="Mean proportion of area of proxy represented by species distribution (%)", breaks = seq(0, 1, by = 0.25), expand = c(0, 0)) + 
  geom_smooth(method=loess, aes(fill=Solution)) +
  theme_bw()+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = c(.99, 0.01), legend.justification = c("right", "bottom"))


scale_fill_discrete(name = "Scenarios", labels = c("SDM (n=116)", "SDM+LZ (n=143)", "SDM+PGD (n=218)", "SDM*PGD (n=5004)", "SDM and PDG as ADMU (n=117)"))
