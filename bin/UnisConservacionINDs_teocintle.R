###


## Needed libraries
library(dplyr)
library(tidyr)
#library(maptools)
#library(raster)


#### Load data

## Teocintle conservation units
cuadro27<-read.delim("../data/CUteocintles/input/cuadro27.tsv") # Cuadro 27 from informe v2 Jesús Sánchez, contains conservations units
cuadro24<-read.delim("../data/CUteocintles/input/cuadro24.tsv") # cuadro 24 Idem, conservation units Z m parviglumis
cuadro25<-read.delim("../data/CUteocintles/input/cuadro25.tsv") # cuadro 25 Idem, Idem Z m mexicana Chalco
cuadro26<-read.delim("../data/CUteocintles/input/cuadro26.tsv") # cuadro 26 Idem, Idem Z m mexicana Mesa Central, Nogogame y Dgo


## Samples used in genetic analyses
individuos_33929<-read.delim("../data/CUteocintles/input/individuos_33929.tsv") # samples in the T4465_33929.h5 original file
individuos_33929plink<-read.delim("../data/CUteocintles/input/plink.fam", header=FALSE)


## Geographic metada
geometa<-read.delim("../data/CUteocintles/input/teocintles_geografia.tsv") # geographic metadata


######## Tidy data

##### Clean geometa

# Filter maize and missing data
levels(geometa$Especie)
geometa<-filter(geometa, Especie!="Z. m. mays", Especie!="Zea mays ssp. mays", Especie!="")
geometa<-droplevels(geometa)
complete.cases(geometa$LATITUDE)

# Fix levels
levels(geometa$Especie)
levels(geometa$Especie)[1]<-"Zea mays ssp. mexicana"
levels(geometa$Especie)[2]<-"Zea mays ssp. parviglumis"
levels(geometa$Especie)


## add "Taxa" column inc Sp + Raza data
# check levels
levels(geometa$Raza)
levels(geometa$Especie)

# add taxa
geometa<-mutate(geometa, Taxa=paste(Especie, Raza))
geometa$Taxa<-as.factor(geometa$Taxa)

#clean and check levels
levels(geometa$Taxa)
levels(geometa$Taxa)[2]<-"Zea diploperennis-Na"
levels(geometa$Taxa)[4]<-"Zea luxurians SFU"
levels(geometa$Taxa)[13]<-"Zea nicaraguensis"
levels(geometa$Taxa)[14]<-"Zea perennis"
levels(geometa$Taxa)[15]<-"Zea perennis-Ur"
levels(geometa$Taxa)

## Parc data should be character
geometa$Parc<-as.character(geometa$Parc)
cuadro27$Parc<-as.character(cuadro27$Parc)
cuadro24$Poblaciones<-as.character(cuadro24$Poblaciones)
individuos_33929$Parc<-as.character(individuos_33929$Parc)



##### Add geographic location to conservation units


## To cuadro 27 conservation units for all Zea
conUnisAll<-left_join(cuadro27, geometa, by=c("Parc" = "Parc")) %>%
           dplyr::select(., Parc:Altitud, Especie, LATITUDE, LONGITUDE, FECHA, Taxa) %>%# keep unly useful  columns
           dplyr::select(., Parc, Unidad, Taxa, Especie, Raza.x:Altitud, LATITUDE:FECHA) %>%
           rename(., Raza=Raza.x, LUGAR=LUGAR.x, MUNICIPIO=MUNICIPIO.x, ESTADO=ESTADO.x)

write.table(conUnisAll, "../data/CUteocintles/output/conUnisAllZea.tsv", quote = FALSE, row.names = FALSE, sep="\t")

## To cuadro 24 conservation units from Zea m parviglumis

# change 100* to 100 (don't know why there is a *)
cuadro24$Poblaciones[8]<-"110"

# tidy and add geometa
conUnisParvi<-dplyr::select(cuadro24, -Estado) %>% # remove Estado column, not useful
           separate_rows(., Poblaciones, sep=", ") %>% # have Pobs as a tidy variable
           rename(., Parc=Poblaciones) %>% # To match with Parc from geometa
           left_join(., geometa, by=c("Parc" = "Parc")) # add geometa data

write.table(conUnisParvi, "../data/CUteocintles/output/conUnisParvi.tsv", quote = FALSE, row.names = FALSE, sep="\t")

## To cuadros 25 y 26 conservation units Z m mexicana

conUnisMexi<-rbind(cuadro25, cuadro26) %>%
             dplyr::select(., -Estado, -No.) %>% # remove Estado and No. column, not useful
             separate_rows(., Poblaciones, sep=", ") %>% # have Pobs as a tidy variable
             rename(., Parc=Poblaciones) %>% # To match with Parc from geometa
             left_join(., geometa, by=c("Parc" = "Parc")) %>% # add geometa data
             mutate(., CUnitNo.=Unidad.de.conservaci.n)

levels(conUnisMexi$CUnitNo.)<-c(1:46)
             
write.table(conUnisMexi, "../data/CUteocintles/output/conUnisMexi.tsv", quote = FALSE, row.names = FALSE, sep="\t")
             

#### Filter h5 sample individuals and add geographic data

## Teocintles full metadata
teocintes33929<-filter(individuos_33929, ANALISIS != "ELIMINADA") %>% #keep only teocintes used for analyses
                dplyr::select(., -ANALISIS) %>%
                left_join(., geometa, by=c("Parc" = "Parc"))
 
write.table(teocintes33929, "../data/CUteocintles/output/teocintes33929_meta.tsv", quote = FALSE, row.names = FALSE, sep="\t")

## Parviglumis full metadata
parvi33929meta<-filter(teocintes33929, Especie=="Zea mays ssp. parviglumis")

write.table(parvi33929meta, "../data/CUteocintles/output/parvi33929_meta.tsv", quote = FALSE, row.names = FALSE, sep="\t")

## Parviglumis only samples names (to filter plink with only parvi)

parvi33929inds<-individuos_33929plink[individuos_33929plink$V2 %in% parvi33929meta$INDIV, ] %>%
                dplyr::select(., V1, V2) #V2 = ind ID, V1 = Family ID

write.table(parvi33929inds, "../data/CUteocintles/output/parvi33929indIDs.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE)


## Mexicana full metadata

mexicana33929meta<-filter(teocintes33929, Especie=="Zea mays ssp. mexicana")

write.table(mexicana33929meta, "../data/CUteocintles/output/mexicana33929_meta.tsv", quote = FALSE, row.names = FALSE, sep="\t")

## Mexicana only samples names (to filter plink with only mexicana)

mexicana33929inds<-individuos_33929plink[individuos_33929plink$V2 %in% mexicana33929meta$INDIV, ] %>%
  dplyr::select(., V1, V2) #V2 = ind ID, V1 = Family ID

write.table(mexicana33929inds, "../data/CUteocintles/output/mexicana33929indIDs.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE)


