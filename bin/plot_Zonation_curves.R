### This script plots the performance curves and proportion of
### proxies of genetid diversity conserved under Zonation
### analysis considering the SDM, PDG, habitat and threat layers


### Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

####               
#### For all Taxa (Figure main text)
####

###### Read data

## Read conservation features names
zfeatures<-read.table("../data/spatial/Zonation_final_solutions/E_final_Todos.features_info.txt",
                      skip=1, header = TRUE)

# extract conservation features names (each SDP subdivided by its PDG)
feature_names<-substring(zfeatures$MapFileName, first= 9) # cut MDP_PDG/
feature_names<-gsub(pattern=".tif", replacement= "", feature_names) # cut .tif

# check names
feature_names[1]
feature_names[length(feature_names)]

## Read IUCN category per taxa
iucn<-read.csv("../data/spatial/Zonation_final_solutions/IUCN_threat_category.csv")

## Read Zonation curves output 
zcurves<-read.table("../data/spatial/Zonation_final_solutions/E_final_Todos.curves.txt",
                  header=FALSE,
                  comment.char = "#") # firt line is comments, not real col names

###### Tidy Zonation curves output
             
## add col names manually
colnames(zcurves) <- c("Prop_landscape_lost", "cost_needed_for_top_fraction",
         "min_prop_rem",  "ave_prop_rem", "W_prop_rem","ext_1", "ext_2", 
         feature_names) # from column 8 onwards column names correspond to feature names
    
# Check few entries and cols of data       
zcurves[1:5, 1:10]


### Tidy data 
   # change to have one column per taxon
zcurves_tidy<-gather(zcurves, feature, prop, 8:ncol(zcurves)) %>%
              separate(feature, c("Nombre_", "PGD"), sep= "_zv_") %>%
   
   # add Area kept 
   mutate(Area_kept=(1-Prop_landscape_lost)*100) %>%
   
   # add iucn data
   left_join(iucn, by ="Nombre_") %>%

   # filter species that were not evaluated in iucn workshop
   filter(IUCN.threat.category!="NO EVALUADA") %>% droplevels(.)

head(zcurves_tidy)


##### Prepare data for plot and plot

## Get mean_prop for all features
all_areakept<-group_by(zcurves_tidy, Nombre_, Area_kept) %>%
             summarise(all=mean(prop))
head(all_areakept)


## Mean prop for all taxa per Area kept

# Plot for each taxa
ggplot(all_areakept, aes(x=Area_kept, y=all, color=Nombre_)) +
  geom_line() + theme(legend.position = "none")

# Plot average and sd
ggplot(all_areakept, aes(x=Area_kept, y=all)) +
  stat_summary(fun.y = mean,
               geom = "ribbon",
               alpha = .3,
               fun.ymax = function(x) mean(x) + sd(x),
               fun.ymin = function(x) mean(x) - sd(x)) +
  stat_summary(fun.y="mean", geom="line")

# Plot mean, median and mix max
ggplot(all_areakept, aes(x=Area_kept, y=all)) +
  stat_summary(fun.y = mean,
               geom = "ribbon",
               alpha = .3,
               fun.ymax = max,
               fun.ymin = min) +
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.y="median", geom = "line", linetype=2)

# Plot by genus and iucn category
# get genus
by_genus<-left_join(all_areakept, iucn, by="Nombre_") %>%
          separate(Nombre_, 
                   into=c("Genus", "Species"), 
                   sep="_", extra="merge") %>%
          droplevels()
                   
head(by_genus)

# colors for Iucn categories
iucn.cols<-c("#D81E05", "#FC7F3F", "#F9E814", "#CCE226", "#60C659", "#D1D1C6") # respectively for "CR", "EN", "VU", "NT", "LC", "DD"

# plot
ggplot(by_genus, aes(x=Area_kept, y=all, color=IUCN.threat.category)) +
       geom_point(size=.7) + facet_wrap( ~ Genus, ncol=2) +
       scale_color_manual(values= iucn.cols,
                          breaks= c("CR", "EN", "VU", "NT", "LC", "DD"),
                          name= "IUCN threat category",
                          guide = guide_legend(override.aes = list(size = 3))) +
  
     # nicer labels
      labs(x="Percentage of Mexican territory", y= "Mean proportion of the area of each PDG within each taxon") +
      theme_bw() + 
      
     # larger text
      theme(strip.text.x = element_text(size=13),
            axis.title = element_text(size=14),
            axis.text = element_text(size=13),
            legend.text = element_text(size=13),
            legend.title = element_text(size=13, face="bold")) 
      
  

## Estimate area kept for category
# CR
CR<-filter(zcurves_tidy, IUCN.threat.category=="CR") %>%
      group_by(Nombre_, Area_kept) %>%
      summarise(prop=mean(prop)) %>%
      mutate(Category="CR")

# EN
EN<-filter(zcurves_tidy, IUCN.threat.category=="EN") %>%
  group_by(Nombre_, Area_kept) %>%
  summarise(prop=mean(prop)) %>%
  mutate(Category="EN")

# DD
DD<-filter(zcurves_tidy, IUCN.threat.category=="DD") %>%
  group_by(Nombre_, Area_kept) %>%
  summarise(prop=mean(prop)) %>%
  mutate(Category="DD")

# LC
LC<-filter(zcurves_tidy, IUCN.threat.category=="LC") %>%
  group_by(Nombre_, Area_kept) %>%
  summarise(prop=mean(prop)) %>%
  mutate(Category="LC")

# NT
NT<-filter(zcurves_tidy, IUCN.threat.category=="NT") %>%
  group_by(Nombre_, Area_kept) %>%
  summarise(prop=mean(prop)) %>%
  mutate(Category="NT")

# VU
VU<-filter(zcurves_tidy, IUCN.threat.category=="VU") %>%
  group_by(Nombre_, Area_kept) %>%
  summarise(prop=mean(prop)) %>%
  mutate(Category="VU")

## Join all together for plot
zc_forplot<-rbind(CR, EN, DD, LC, NT, VU)
head(zc_forplot)


### Plot

# plot
ggplot(zc_forplot, aes(x=Area_kept, y=prop, color=Category)) +
  stat_summary(fun = mean, geom = "line", size=1.5) +
  #stat_summary(fun = median, geom = "line", linetype=2) +
  scale_color_manual(values= iucn.cols,
                     breaks= c("CR", "EN", "VU", "NT", "LC", "DD"),
                     name= "IUCN threat category",
                     guide = guide_legend(override.aes = list(size = 3))) +
  # nicer labels
  labs(x="Percentage of Mexican territory", y= "Mean proportion of the area of each PDG within each taxon") +
 
   # nicer background
  theme_bw() +
  
  # larger text
  theme(
        axis.title = element_text(size=14),
        axis.text = element_text(size=13),
        legend.text = element_text(size=13),
        legend.title = element_text(size=13, face="bold")) 


####               
#### Scenarios with all taxa, taxa exclusively distributing in natural vegetation, and taxa associated to different habitats, i.e. natural vegetation, agricultural and urban areas.
#### (Figure Sup Mats)

## Function to transform the Zonation curves to the right format for plotting
generate_zc_forplot<-function(input_features, input_zcurves){
  ### This function as input the csv of Zonation features and curves,
  ### as well as the list of IUCN category for each taxa to
  ### transform them to the right format to plot the curves grouping taxa by IUCN category 

## Read IUCN category per taxa
  iucn<-read.csv("../data/spatial/Zonation_final_solutions/IUCN_threat_category.csv")
  
## Read conservation features names
zfeatures<-read.table(input_features,
                      skip=1, header = TRUE)

# extract conservation features names (each SDP subdivided by its PDG)
feature_names<-substring(zfeatures$MapFileName, first= 9) # cut MDP_PDG/
feature_names<-gsub(pattern=".tif", replacement= "", feature_names) # cut .tif

## Read Zonation curves output 
zcurves<-read.table(input_zcurves,
                    header=FALSE,
                    comment.char = "#") # firt line is comments, not real col names


###### Tidy Zonation curves output

## add col names manually
colnames(zcurves) <- c("Prop_landscape_lost", "cost_needed_for_top_fraction",
                       "min_prop_rem",  "ave_prop_rem", "W_prop_rem","ext_1", "ext_2", 
                       feature_names) # from column 8 onwards column names correspond to feature names

### Tidy data 
# change to have one column per taxon
zcurves_tidy<-gather(zcurves, feature, prop, 8:ncol(zcurves)) %>%
  separate(feature, c("Nombre_", "PGD"), sep= "_zv_") %>%
  
  # add Area kept 
  mutate(Area_kept=(1-Prop_landscape_lost)*100) %>%
  
  # add iucn data
  left_join(iucn, by ="Nombre_") %>%
  
  # filter species that were not evaluated in iucn workshop
  filter(IUCN.threat.category!="NO EVALUADA") %>% droplevels(.)

## Estimate area kept for category
# CR
CR<-filter(zcurves_tidy, IUCN.threat.category=="CR") %>%
  group_by(Nombre_, Area_kept) %>%
  summarise(prop=mean(prop)) %>%
  mutate(Category="CR")

# EN
EN<-filter(zcurves_tidy, IUCN.threat.category=="EN") %>%
  group_by(Nombre_, Area_kept) %>%
  summarise(prop=mean(prop)) %>%
  mutate(Category="EN")

# DD
DD<-filter(zcurves_tidy, IUCN.threat.category=="DD") %>%
  group_by(Nombre_, Area_kept) %>%
  summarise(prop=mean(prop)) %>%
  mutate(Category="DD")

# LC
LC<-filter(zcurves_tidy, IUCN.threat.category=="LC") %>%
  group_by(Nombre_, Area_kept) %>%
  summarise(prop=mean(prop)) %>%
  mutate(Category="LC")

# NT
NT<-filter(zcurves_tidy, IUCN.threat.category=="NT") %>%
  group_by(Nombre_, Area_kept) %>%
  summarise(prop=mean(prop)) %>%
  mutate(Category="NT")

# VU
VU<-filter(zcurves_tidy, IUCN.threat.category=="VU") %>%
  group_by(Nombre_, Area_kept) %>%
  summarise(prop=mean(prop)) %>%
  mutate(Category="VU")

## Join all together for plot
zc_forplot<-rbind(CR, EN, DD, LC, NT, VU)
return(zc_forplot)
}


### Generate zonation curves for plot for each dataset

# all
zc_forplot_all<-generate_zc_forplot(input_features<-"../data/spatial/Zonation_final_solutions/E_final_Todos.features_info.txt",
                                    input_zcurves<-"../data/spatial/Zonation_final_solutions/E_final_Todos.curves.txt")

# taxa exclusively distributing in natural vegetation
zc_forplot_Veg<-generate_zc_forplot(input_features<-"../data/spatial/Zonation_final_solutions/E_final_VegPyS.features_info.txt",
                                    input_zcurves<-"../data/spatial/Zonation_final_solutions/E_final_VegPyS.curves.txt")

# taxa associated to different habitats, i.e. natural vegetation, agricultural and urban areas.
zc_forplot_HabVarios<-generate_zc_forplot(input_features<-"../data/spatial/Zonation_final_solutions/E_final_HabVarios.features_info.txt",
                                    input_zcurves<-"../data/spatial/Zonation_final_solutions/E_final_HabVarios.curves.txt")


### Plot 
# colors for Iucn categories
iucn.cols<-c("#D81E05", "#FC7F3F", "#F9E814", "#CCE226", "#60C659", "#D1D1C6") # respectively for "CR", "EN", "VU", "NT", "LC", "DD"

# plot a) all taxa
plot.a<-ggplot(zc_forplot_all, aes(x=Area_kept, y=prop, color=Category)) +
  stat_summary(fun = mean, geom = "line", size=1.5) +
  scale_color_manual(values= iucn.cols,
                     breaks= c("CR", "EN", "VU", "NT", "LC", "DD"),
                     name= "IUCN threat category",
                     guide = guide_legend(override.aes = list(size = 3))) +
  # nicer labels
  labs(x="", y= "") +
  ggtitle("a)") +
  
  # nicer background
  theme_bw() +
  
  # larger text
  theme(
    plot.title = element_text(size=16, face="bold"),
    axis.title = element_text(size=14),
    axis.text = element_text(size=13),
    legend.text = element_text(size=13),
    legend.title = element_text(size=13, face="bold"))

# plot b) taxa exclusively distributing in natural vegetation
plot.b<-ggplot(zc_forplot_Veg, aes(x=Area_kept, y=prop, color=Category)) +
  stat_summary(fun = mean, geom = "line", size=1.5) +
  scale_color_manual(values= iucn.cols,
                     breaks= c("CR", "EN", "VU", "NT", "LC", "DD"),
                     name= "IUCN threat category",
                     guide = guide_legend(override.aes = list(size = 3))) +
  # nicer labels
  labs(x="", y= "Mean proportion of the area of each PDG within each taxon") +
  ggtitle("b)") +
  
  # nicer background
  theme_bw() +
  
  # larger text
  theme(
    plot.title = element_text(size=16, face="bold"),
    axis.title = element_text(size=14),
    axis.text = element_text(size=13),
    legend.text = element_text(size=13),
    legend.title = element_text(size=13, face="bold")) 

# plot for taxa associated to different habitats, i.e. natural vegetation, agricultural and urban areas.
plot.c<-ggplot(zc_forplot_HabVarios, aes(x=Area_kept, y=prop, color=Category)) +
  stat_summary(fun = mean, geom = "line", size=1.5) +
  scale_color_manual(values= iucn.cols,
                     breaks= c("CR", "EN", "VU", "NT", "LC", "DD"),
                     name= "IUCN threat category",
                     guide = guide_legend(override.aes = list(size = 3))) +
  # nicer labels
  labs(x="Percentage of Mexican territory", y= "") +
  ggtitle("c)") +
  
  # nicer background
  theme_bw() +
  
  # larger text
  theme(
    plot.title = element_text(size=16, face="bold"),
    axis.title = element_text(size=14),
    axis.text = element_text(size=13),
    legend.text = element_text(size=13),
    legend.title = element_text(size=13, face="bold"))

## Multiplot

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Multiple plot function from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
  #
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  #
    library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot(plot.a, plot.b, plot.c, cols=1)

#### Get session info
sessionInfo()