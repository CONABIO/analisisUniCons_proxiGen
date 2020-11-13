# Load needed libraries

library(raster)
library(dplyr)
library(ggplot2)
library(tidyr)
library(vegan)

# Load spatial data

# Path to rasters
myrasters<-c("../data/spatial/Zonation_output/E1_02_caz.rank.compressed1.tif", # Zonation output proxis div gen
             "../data/spatial/Zonation_output/03_caz.rank.compressed.tif", # Zonation output Zonas Vida only
             "../data/spatial/Zonation_output/02_caz.rank.compressed.tif", # Zonation output with species only
             "../data/spatial/modelosDarwinZea/Zea_mays_parviglumis.tif", #  sp model
             "../data/spatial/areasProxyDivGen/tif/proxiesDivGen_merged/ProxisDivGen_merged_masked.tif") #Proxies div gen


# Stack rasters and add nicer names 
CWR_sp<-stack(myrasters)
names(CWR_sp)<-c("ZonationProxiDG", "ZonationZV", "ZonationSPP", "ENM_sp", "ProxiesDivGen")


##### Analyses #####

#### Basic plots #####

#colors for proxies div gen 
mycols<-c("#668325", "#9343c1", "#61b928", "#bc6cec", "#56cd5e", "#4667e6", "#bab71b", "#704dc1", "#b5cb3e", "#4c4fb8", "#83b124", "#5e7ef2", "#38a435", "#cc3ca6", "#3fcb85", "#ec49a6", "#86ca62", "#e471de", "#5d9d38", "#9b76e8", "#9ea122", "#a13ca0", "#3da35c", "#e53283", "#2b7c2c", "#b063c4", "#c7bb43", "#3b65c1", "#e2b438", "#628ceb", "#e69128", "#509be7", "#e0422f", "#34c7dd", "#c62a33", "#57d2b9", "#d1204a", "#57c597", "#ef3d76", "#56cdc9", "#ec4b5f", "#2ba198", "#e2672a", "#2baad3", "#b33f1d", "#6bbff0", "#ae6314", "#7f68c4", "#90ac42", "#774598", "#c0bd59", "#b58be8", "#4e6313", "#dd80d3", "#87be71", "#b62a7c", "#3b9065", "#c22a60", "#71b588", "#e75687", "#326a38", "#e46aad", "#628848", "#ad498c", "#b2c06d", "#535a99", "#c89d3c", "#3176ae", "#e97150", "#549bce", "#ca783f", "#7b7ec0", "#927b1f", "#a6a5e8", "#9a9c49", "#a161a4", "#a6b174", "#7d4f8c", "#e99b5c", "#206e54", "#b1394e", "#505e25", "#dc9ce0", "#786e30", "#ba7eb8", "#90652d", "#de93ba", "#caa569", "#9b3a66", "#e99d7a", "#8c4e6c", "#914c27", "#f288aa", "#b7443b", "#c6648b", "#bc7257", "#dd6e81", "#a7514b", "#e49190", "#a04f59", "#e6736d")


plot(CWR_sp$ProxiesDivGen, main= "Areas Proxi of Genetic Diversity", col=mycols, legend=FALSE)
plot(CWR_sp$ZonationSPP>=0.98, legend=FALSE, main= "Zonation output Spp models")
plot(CWR_sp$ZonationZV>=0.98, legend=FALSE, main= "Zonation output Life Zones")
plot(CWR_sp$ZonationProxiDG>=0.98, legend=FALSE, main= "Zonation output Proxies Genetic Diversity")
plot(CWR_sp$ENM_sp, legend=FALSE, main= "Z. m. parviglumis distribution")

## Describe proxies div gen
table(getValues(CWR_sp$ProxiesDivGen))

# Define areas with Proxies of DivGen for that sp
# this works because the ENM_sp is coded 1= sp presence 0= absence
sp_ProxiesDivGen<- CWR_sp$ProxiesDivGen * CWR_sp$ENM_sp

# write to raster
writeRaster(sp_ProxiesDivGen, format= "GTiff", file="../data/spatial/areasProxyDivGen/crop_to_sp/parvi_ProxiesDivGen")

###### Examine overlaps #####

#### Plot MAPs with overlaps 

## Zoom to sp extent

# Define sp extent
#plot(CWR_sp$ENM_sp)
#e <- drawExtent()
e<-extent(c(-11814531, -10715683 , 1760353, 2699175))

## plot Proxies div gen  and sp distribution
plot(CWR_sp$ProxiesDivGen, ext=e, col=mycols)

# set colors same than for all Mx
mycols3<-mycols[unique(getValues(sp_ProxiesDivGen))[3:length(unique(getValues(sp_ProxiesDivGen)))]] 

## plot maps for overlaping each Zonation solution

# Proxies in sp distribution
plot(sp_ProxiesDivGen, ext=e, col=c("grey94", mycols3), legend=FALSE, main= "Areas Proxi of Genetic Diversity for Z. m. parviglumis")

# Zonation Spp solution
plot(sp_ProxiesDivGen, ext=e, col=c("grey94", mycols3), legend=FALSE, main= "Z. m. parviglumis conserved under ZonationSPP 10%")
plot(CWR_sp$ZonationSPP, ext=e, add=TRUE, breaks=c(0, 0.98, 1), oldstyle=TRUE, col=c("#00000000", rgb(0,0,0,.4)), legend = FALSE)

# Zonation ZV solution
plot(sp_ProxiesDivGen, ext=e, col=c("grey94", mycols3), legend=FALSE, main= "Z. m. parviglumis conserved under ZonationZV 10%")
plot(CWR_sp$ZonationZV, ext=e, add=TRUE, breaks=c(0, 0.98, 1), oldstyle=TRUE, col=c("#00000000", rgb(0,0,0,.4)), legend = FALSE)

# Zonation Proxies solution
plot(sp_ProxiesDivGen, ext=e, col=c("grey94", mycols3), legend=FALSE, main= "Z. m. parviglumis conserved under ZonationZV 10%")
plot(CWR_sp$ZonationProxiDG, ext=e, add=TRUE, breaks=c(0, 0.98, 1), oldstyle=TRUE, col=c("#00000000", rgb(0,0,0,.4)), legend = FALSE)



################### Quantify overlaps

### convert to rasters to points to have a "df of overlaps"
CWR_overlaps<-as.data.frame(rasterToPoints(CWR_sp))

##### calcs for each solution

##### 1) Proxies of Div Gen in base ENM of sp
## subset
x<-dplyr::filter(CWR_overlaps, ENM_sp==1) %>% arrange(., ProxiesDivGen)
x$ProxiesDivGen<-as.character(x$ProxiesDivGen)

# estimates
which_sp_proxies<-unique(x$ProxiesDivGen)
which_sp_proxies

n_sp_proxies<-length(unique(x$ProxiesDivGen))
n_sp_proxies

# save results in df
sol_base<-as.data.frame(table(x$ProxiesDivGen))
names(sol_base)<-c("ProxiDivGen", "ENM_sp")
sol_base

### plots

## bar plot 
ggplot(x, aes(ProxiesDivGen)) + geom_bar() + ylim(0, 35000) + ggtitle("Areas Proxi of Genetic Diversity for Z. m. parviglumis")

# Nicer Plot
p<- ggplot(x, aes(ProxiesDivGen)) + geom_bar() + ylim(0, 35000) + ggtitle("Areas Proxi of Genetic Diversity for Z. m. parviglumis")

p + theme(axis.title.x = element_text(face="bold", size=20),
          axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
          axis.title.y = element_text(face="bold", size=20),
          axis.text.y  = element_text(size=16)) 



##### 2) CWRsp + Zonation spp models

### calcs
## subet
x<-dplyr::filter(CWR_overlaps, ENM_sp==1, ZonationSPP>=0.90) %>% arrange(., ProxiesDivGen)
x$ProxiesDivGen<-as.character(x$ProxiesDivGen)

# estimates
which_sp_proxies<-unique(x$ProxiesDivGen)
which_sp_proxies

n_sp_proxies<-length(unique(x$ProxiesDivGen))
n_sp_proxies


# Save to common data frame
sol_ZSPP<-as.data.frame(table(x$ProxiesDivGen))
names(sol_ZSPP)<-c("ProxiDivGen", "ZonationSPP10")
sol_ZSPP

# set colors same than for all Mx
mycols3<-mycols[unique(getValues(sp_ProxiesDivGen))[3:length(unique(getValues(sp_ProxiesDivGen)))]] 

## bar plot 
ggplot(x, aes(ProxiesDivGen)) + geom_bar() + ylim(0, 35000) + ggtitle("Z. m. parviglumis conserved under ZonationSPP 10%")




##### 3) CWRsp + Zonation Zonas de vida

### calcs
## subet
x<-dplyr::filter(CWR_overlaps, ENM_sp==1, ZonationZV>=0.90) %>% arrange(., ProxiesDivGen)
x$ProxiesDivGen<-as.character(x$ProxiesDivGen)

# estimates
which_sp_proxies<-unique(x$ProxiesDivGen)
which_sp_proxies

n_sp_proxies<-length(unique(x$ProxiesDivGen))
n_sp_proxies


# Save to common data frame
sol_ZZV<-as.data.frame(table(x$ProxiesDivGen))
names(sol_ZZV)<-c("ProxiDivGen", "ZonationZV10")
sol_ZZV


# set colors same than for all Mx
mycols3<-mycols[unique(getValues(sp_ProxiesDivGen))[3:length(unique(getValues(sp_ProxiesDivGen)))]] 

## bar plot 
ggplot(x, aes(ProxiesDivGen)) + geom_bar() + ylim(0, 35000) + ggtitle("Z. m. parviglumis conserved under ZonationZV 10%")



##### 4) CWRsp + Zonation Zonas Proxies Div Gen

### calcs
## subet
x<-dplyr::filter(CWR_overlaps, ENM_sp==1, ZonationProxiDG>=0.90) %>% arrange(., ProxiesDivGen)
x$ProxiesDivGen<-as.character(x$ProxiesDivGen)

# estimates
which_sp_proxies<-unique(x$ProxiesDivGen)
which_sp_proxies

n_sp_proxies<-length(unique(x$ProxiesDivGen))
n_sp_proxies


# Save to common data frame
sol_ZProxiDG<-as.data.frame(table(x$ProxiesDivGen))
names(sol_ZProxiDG)<-c("ProxiDivGen", "ZonationProxiDG10")
sol_ZProxiDG


# set colors same than for all Mx
mycols3<-mycols[unique(getValues(sp_ProxiesDivGen))[3:length(unique(getValues(sp_ProxiesDivGen)))]] 

## bar plot 
ggplot(x, aes(ProxiesDivGen)) + geom_bar() + ylim(0, 35000) + ggtitle("Z. m. parviglumis conserved under ZonationProxiDivGen 10%")


####### Analyse results together


### ZonationSPP and ProxiDivGen
## Prepare data
# Join solutions
sol_2<-left_join(sol_ZSPP, sol_ZProxiDG) %>% 
  mutate_all(funs(replace(., is.na(.), 0))) # NAs are actually 0s
head(sol_2)

# gather into tidy df
sol_tidy<-gather(sol_2, Solution, Freq, -ProxiDivGen)
head(sol_tidy)


## Plot
p<- ggplot(sol_tidy, aes(ProxiDivGen, Freq)) + geom_col(aes(fill=Solution))
p

p + theme(axis.title.x = element_text(face="bold", size=20),
          axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
          axis.title.y = element_text(face="bold", size=20),
          axis.text.y  = element_text(size=16)) + 
  
  scale_fill_discrete(name="Zonation output",
                      breaks=c("ZonationProxiDG10", "ZonationSPP10"),
                      labels=c("Using proxies of genetic diversity", "Using SDM models only")) +
  
  theme(legend.title = element_text(face="bold", size=20),
        legend.text = element_text(size=18))


### ALL
## Prepare data
# Join all solutions
sol_all<-left_join(sol_base, sol_ZSPP) %>% left_join(., sol_ZZV) %>% left_join(., sol_ZProxiDG) %>% 
  mutate_all(funs(replace(., is.na(.), 0))) # NAs are actually 0s
head(sol_all)

# gather into tidy df
sol_tidy<-gather(sol_all, Solution, Freq, -ProxiDivGen)
head(sol_tidy)

## Plot
ggplot(sol_tidy, aes(ProxiDivGen, Freq)) + geom_col(aes(fill=Solution))

## Estimate Proportion of Zones Proxies Div Gen under conservation for each Zonation solution

# estimate proportion
sol_prop<-mutate(sol_all, ZonationSPP10=ZonationSPP10/ENM_sp, 
                 ZonationZV10=ZonationZV10/ENM_sp,
                 ZonationProxiDG10=ZonationProxiDG10/ENM_sp) %>% select(-ENM_sp)
head(sol_prop)

# tidy
sol_prop<-gather(sol_prop, Solution, Proportion_ProxiDG, -ProxiDivGen)
head(sol_prop)

# view summary
sol_prop<-group_by(sol_prop, Solution)
summarise(sol_prop, mean.prop = mean(Proportion_ProxiDG))
summarise(sol_prop, mean.prop = mean(Proportion_ProxiDG), sum.prop=sum(Proportion_ProxiDG))

### Estimate Simpson diversity index

## we need an abundance matrix with species (ProxisDivGen in this case) in columns
# and abundance by treatment (Zonation output in this case) in rows. Therefore:

# transpose data to produce abundace matrix as expected by vegan
row.names(sol_all)<-sol_all[,1]
abundance.matrix<-t(sol_all[, 2:5])
class(abundance.matrix)

### Estimate Indices

indices.div<-data.frame("zonification"=colnames(sol_all[2:5]))

# Richness
indices.div$Richness<-rowSums(abundance.matrix>0)

# Shannon
indices.div$shannon<-diversity(abundance.matrix, index="shannon") # shannon is default

# Simpson
indices.div$simpson<-diversity(abundance.matrix, index="simpson")

## see results
indices.div
