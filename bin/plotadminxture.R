library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(scales)


## Parviglumis

### Read meta data
meta<-read.delim("../data/genetic/output/ADN_pasap_3604.txt")
parviPlinksamples<-read.table("../data/genetic/output/bytaxa_parviglumis.fam")
parvi_CV<-read.delim("../data/genetic/output/admixture/parvi/Kerror_parviglumis.txt", 
                    sep=" ", header=FALSE, stringsAsFactors = FALSE) %>%
                    mutate(K = parse_number(V3)) %>% 
                    select(K, V4) %>% rename(CV=V4)

## Plot CV
p <- ggplot(parvi_CV, aes(x=K, y=CV)) + geom_point(size=1) + geom_line()  

p + theme(axis.title.x = element_text(face="bold", size=20),
          axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
          axis.title.y = element_text(face="bold", size=20),
          axis.text.y  = element_text(size=16)) 

### Get Q data and plot

for(K in c(1:42)){
## get data
  
Qval<-read.table(paste0("../data/genetic/output/admixture/parvi/bytaxa_parviglumis.", K, ".Q"))
names(Qval)<-paste0("K", 1:ncol(Qval))

## Tidy data for ggplot
# get metadata of plink samples and add to Qtable
parvimeta<-meta[meta$POBL %in% parviPlinksamples$V2, ]
Qval<-cbind(parvimeta$POB, parvimeta$INDIV, parvimeta$LONGITUDE, parvimeta$LATITUDE, Qval)
names(Qval)[1:4]<-c("POB", "INDIV", "LONGITUDE", "LATITUDE") 

# write file
write.table(Qval, file=paste0("../data/genetic/output/admixture/parvi/Qparvi_K", K, ".tsv"), sep="\t", quote = FALSE, row.names = FALSE)

# transfor to long format  
Qval_long<-gather(Qval, key=Kgroup, value=Qadmixture, 5:ncol(Qval))

## Plot
# standar
plt<-ggplot(Qval_long, aes(x=INDIV, y=Qadmixture, fill=Kgroup)) + geom_col() +
     theme(axis.text.x= element_blank())

# get color names in case it is needed match with other plots outside r
print("colors used")
print(scales::hue_pal()(K))

print(plt)     
  }


