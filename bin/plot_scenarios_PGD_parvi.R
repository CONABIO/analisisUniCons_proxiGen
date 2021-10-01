# Libraries
library(dplyr)
library(tidyverse)
library(forcats)


# Path to files
spatial_output <- read.table("../data/comparations_output/20_sol_prop_spp.txt", 
                    header=T, sep = ",")
head(spatial_output)

genetic_output <- read.table("../data/genetic/output/PGD_points_parvi.txt", 
            header=T, sep = ",")


# organize genetic data
PGD <- unique(genetic_output[c("PGD_points")]) #obtain unique PGD
colnames(PGD) <- "PGD" #change col names
as.data.frame(PGD)
str(PGD)


# organize spatial data
data <- as.data.frame(spatial_output) %>% 
  filter(Solution %in% c("Scenario_SDM", "Scenario_SDM_vs_PGD"))

colnames(data)[1] <- "PGD" #change col names


# include genetic and spatial data 
data2 <- inner_join(data,PGD,by="PGD") # join genetic and spatial df


# calculate mean rep values of PGD for each scenario
data3 <- aggregate(x=data2$Proportion_ProxiDG,
          by=list(data2$PGD,data2$Solution),
          FUN=mean) 
colnames(data3)[1] <- "PGD" 
colnames(data3)[2] <- "Solution" 
colnames(data3)[3] <- "Mean_Prop" 
head(data3)

# A grouped bar chart
legend_title <- "Scenario" 
ggplot(data3, aes(fill=Solution, y=Mean_Prop, x=factor(PGD))) + 
  geom_bar(position="dodge", stat="identity") +
  labs(x = "ID Proxy of genetic differentiation", 
       y = "Mean proporton of the area of \nZ. m. parviglumis within each PGD (%)")+
  theme(text = element_text(size=9)) + 
  scale_fill_hue(legend_title, labels = c("SDM", "SDM*PGD"))+
  theme_bw() +
  theme(legend.position = "top", 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
        ) 

### change order of PGD by Mean_Prop and plot again
# get order
x<-with(data3, data3[order(Solution, Mean_Prop),]) 
x
PGDs_ordered<-as.character(x[x$Solution=="Scenario_SDM", 1])

# order levels
data4<-data3
data4$PGD<-factor(data4$PGD, levels=PGDs_ordered)
data4$PGD

# plot
legend_title <- "Scenario" 

p<-ggplot(data4, aes(fill=Solution, y=Mean_Prop, x=factor(PGD))) + 
  geom_bar(position="dodge", stat="identity") +
  labs(x = "PGD", 
       y = "Mean proportion (%)")+
  scale_fill_hue(legend_title, labels = c("SDM", "SDM*PGD"))+
  theme_bw() +
  theme(legend.position = "top", 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=15)) 
p

ggsave("../figures/PGDmap_barplotPGD_parvi.png",
       plot=p, dpi=300,
       width=10, height = 10, units="cm")
