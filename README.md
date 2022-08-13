# README

Code and settings for the analyses of the paper **_Incorporating evolutionary and threat processes into crop wild relatives conservation_** by Wolke Tobón-Niedfeldt, Alicia Mastretta-Yanes, Tania Urquiza-Haas, Bárbara Goettsch, Angela P. Cuervo-Robayo, Esmeralda Urquiza-Haas, M. Andrea Orjuela-R, Francisca Acevedo Gasman, Oswaldo Oliveros-Galindo, Caroline Burgeff, Diana Rivera-Rodríguez, José de Jesús Sánchez González, Jesús Alarcón-Guerrero, Araceli Aguilar-Meléndez, Flavio Aragón Cuevas, Valeria Alavez, Gabriel Alejandre-Iturbide, Carlos-H. Avendaño-Arrazate, César Azurdia Pérez, Alfonso Delgado-Salinas, Pablo Galán, Manuel González-Ledesma, Jesús Hernández-Ruíz, Francisco G. Lorea-Hernández, Rafael Lira Saade, Aarón Rodríguez, Dagoberto Rodríguez Delcid, José Ariel Ruiz-Corra, Juan José Santos Pérez, Ofelia Vargas-Ponce, Melania Vega, Ana Wegier, Martín Quintana-Camargo, José Sarukhán and
Patricia Koleff.

The data used by the scripts below is provided in a separate Dryad repository ().

### Dependencies

Analyses were carried out in Zonation version 4 and in R version 3.5.1. The following R packages were used:
`purrr_0.3.4`, `tidyr_1.0.2`   `dplyr_1.0.2`   `ggplot2_3.3.3`, `readr_1.4.0`, `gridExtra_2.3`,  `ggnewscale_0.4.5`, `scatterpie_0.1.5`, `pophelper_2.3.1`, `rgdal_1.4-8`, `raster_3.4-5`, `sp_1.4-4`, `rgl_0.107.10` and `pcadapt_4.3.3`. Additional dependencies can be consulted in the html notebook of the r scritps (session info section).


## Analyses and Figures made in R

Code for analyses and plots made in R is available in the `/bin` directory of this repository in the .R and .Rmd scripts detailed below. For each .Rmd there is a .html version of it with an R notebook showing the code and output.


### `PCAdapt_PGD_teocintles.Rmd`

Performs the PCA with the genetic data of *Z. m. parviglumis*. 

Used to produce Figure 3c (3D PCA of genetic data).

**Data used by this script:**

* Genetic data for PCA: `../data/genetic/output/bytaxa_parviglumis.bed`
* Samples metadata including PGD where they fell: `../data/genetic/output/PGD_points_parvi.txt"`

**Data available at the Dryad repository as:** `genetic.tar.gz`.


### `plot_admixture_PGD_teocintles.Rmd`

Explores admixture results to perform: 

* Figure 3a (admixture plot) and Figure 3b (map with points).

* Supplementary Figure 5 (map with admixture pies).

**Data used by this script:**

* *Zea mays parviglumis* SDM: `/data/spatial/modelosDarwinZea/Zea_mays_parviglumis.tif`

* proxies of genetic differentiation croped for *Z. m. parviglumis* SDM:`"../data/spatial/areasProxyDivGen/crop_to_sp/PGD_Zea_mays_parviglumis.tif"`

* *Z. m. parviglumis* CV error (admixture output): `"../data/genetic/output/admixture/parvi/Kerror_parviglumis.txt"`

* *Z. m. parviglumis* Q files for K= 25 and K= 13 (admixture output):
`"../data/genetic/output/admixture/parvi/bytaxa_parviglumis.25.Q"` and `"../data/genetic/output/admixture/parvi/bytaxa_parviglumis.13.Q"`.

* Passport data for each sample: `"../data/genetic/output/ADN_pasap_3604.txt"`

* Fam file with sample list: `"../data/genetic/output/bytaxa_parviglumis.fam"`

**Data available at the Dryad repository as:** `genetic.tar.gz`.

### `plot_scenarios_PGD_parvi.R`

Plots the barplot of Figure 3d.

**Data used by this script:**

* Spatial output of zonation solution for 20% Mexico: `/data/comparations_output/sol_prop_spp.txt`

* Samples metadata including PGD where they fell: `/data/genetic/output/PGD_points_parvi.txt`

**Data available at the Dryad repository as:** `Comparing_output_txt.zip` and `genetic.tar.gz`.

The admixture plot, map, PCA and bar plot were joined using Inkscape to make:

* Figure 3. Genetic diversity of Zea mays subsp. parviglumis, a maize CWR, represented in the proxies of genetic differentiation (PGD).


### `plot_scenarios_PDG.R`

Performs the plot an analyses:

* Figure 4. Performance of five systematic conservation planning scenarios to represent conservation features of Mesoamerican crop wild relatives, considering 20% of Mexico’s terrestrial area.

* Supplementary analyses suggested during review (only stats are mentioned in the ms text)

**Data used by this script:** `/data/comparations_output/sols_summary_spp.txt` produced by the script `spatial_analyses_zonationVSproxiesdivgen_all_ms.R`

**Data available at the Dryad repository as:** `Comparing_output_txt.zip`

### `plot_Zonation_curves.R`

Performs the analyses and plots of the following figures:

* Figure 6. Performance curves quantifying the representation of proxies of genetic differentiation within the distribution of crop wild relatives in Mexico, based on the hierarchical landscape priority rank map.

* Supplementary Figure 9. Performance curves of how the proportion of taxa distribution ranges increased with the amount of land area. 

**Data used by this script:**

* Read IUCN category per taxa: `/data/spatial/Zonation_final_solutions/IUCN_threat_category.csv"`. 

* Zonation curves output and features for all taxa: `/data/spatial/Zonation_final_solutions/E_final_Todos.curves.txt` and `/data/spatial/Zonation_final_solutions/E_final_Todos.features_info.txt`.  

* Zonation curves output and features but for taxa exclusively distributing in natural vegetation: `/data/spatial/Zonation_final_solutions/E_final_VegPyS.features_info.txt` and `/data/spatial/Zonation_final_solutions/E_final_VegPyS.curves.txt`.

* Zonation curves output and features but for taxa associated to different habitats, i.e. natural vegetation, agricultural and urban areas: `"../data/spatial/Zonation_final_solutions/E_final_HabVarios.features_info.txt"` and `"../data/spatial/Zonation_final_solutions/E_final_HabVarios.curves.txt"`.

**Data available at the Dryad repository as:** `IUCN_threat_category.csv` and `Zonation_output_final.tar.gz`

### `spatial_analyses_zonationVSproxiesdivgen_all_ms.R`

For each taxon crops the proxies of gen div to the species distribution models. Then performs analyses to compare different ways to incorporate the PDG into the Zonation analyses. 

**Data used by this script:**

* Zonation outputs and proxies of gen div:

```
/data/spatial/Zonation_output/01_MDP.rank.compressed.tif, # Zonation output with species only (n=116) #1
/data/spatial/Zonation_output/02_MDP_ZV.rank.compressed.tif, # Zonation output with species + LZ (n=143) #2
/data/spatial/Zonation_output/03_MDP_PDG.rank.compressed.tif, # Zonation output with species and PDG (n=218) #3
/data/spatial/Zonation_output/04_MDP_vs_PDG.rank.compressed.tif, # Zonation output with species vs PDG (n=5004) #4
/data/spatial/Zonation_output/05_MDP_PDG_ADMU.rank.compressed.tif, # Zonation output with species and PDG as ADMU (n=116+1) #5
/data/spatial/areasProxyDivGen/PDG.tif") #Proxies div gen 
```

* Species distribution models of each species: 
`/data/spatial/modelos_darwin_all_final/*.tif`

**Input data available at the Dryad repository as:** `Zonation_output_final.tar.gz`, `PGD_rasters.tar.gz` and `SDM_rasters.tar.gz`.

Output of this script:

* The rasters of proxies of genetic differentiation cropped for each taxa: `/data/spatial/areasProxyDivGen/crop_to_sp/`.

* Area of proxies of div gen for each zonation solution: `/data/comparations_output/sol_tidy_spp.txt`

* Proportion of proxi of div gen for each zolnation solution in relation to its area in the sp distribution: `/data/comparations_output/sol_prop_spp.txt`

* Diversity indexes and mean of proportion: `/data/comparations_output/sols_summary_spp.txt`

**Output data available at the Dryad repository as:** `PGD_croptoSDM_rasters.tar.gz` and `Comparing_output_txt.zip`


## Zonation analyses
Zonation configuration files are in:`/bin/Zonation_config_files`.

File contents as follows:

* Batch files, contain the command line to call Zonation in the cluster: `.sh` files

* Run settings file, defines the analysis settings: `E_final.dat`

* Biodiversity feature files, listing the biodiversity layers: `.spp` files

* Condition file; listing the habitat layers: `habitat_features.txt`

* Groups file; links the biodiversity features to habitat features using habitat groups: `habitat_group_features.txt`

* Species of Special Interest file; lists observation records for species without potential distribution model: `SSI_list.txt`

* Zonation configuration file: `zonation_config.txt`.

Data needed for the Zonation analyses are available at the Dryad repository as: `Zonation_input_conservationfeatures.tar.gz`, `Zonation_input_habitat.tar.gz` and `Zonation_input_spp_points.zip`.

Note: the taxa *Physalis ixocarpa* and *Solanum endiense* are mentioned in the spp list, but where excluded (0 in the first column) from the final analyses because the first turned to be a synonym and the second was not evaluated by the ICUN.

## Non-coding analyses and Figures

Figures 2, 5 and Supplementary Figures 1,2,4,8,10, 12, 13 and 14 were made in ArcMap with the spatial data and Zonation raster outputs available in the Dryad repository.

## Species distribution modelling

Species distribution modelling was done by a previous study using the custom function `MNE.R`, made by Juan Barrios and available at [https://github.com/AngelaCrow/MNE_ParientesSilvestresCultivos](https://github.com/AngelaCrow/MNE_ParientesSilvestresCultivos). See that repository for further details.

Data used:

* CWR ocurrence points in .csv. These were dowloaded for each taxon from the Sistema Nacional de Información sobre Biodiversidad de México [(SNIB)](http://www.snib.mx/) and [Global Biodiversity Information Facility](https://www.gbif.org/) and curated as explained in Methods.
* 19 bioclimatic variables form [Worlclim Version 1.4](http://worldclim.org/) from **1**.
* Terrestial ecoregions proposed by [Olson et al. en el 2001](https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world)
* Evapotranspiration, aridity index, annual radiation data from **2-3**.
* Bare soil and cultivated areas from **4**.

1.	Hijmans, R. J., Cameron, S. E., Parra, J. L., Jones, P. G. & Jarvis, A. Very high resolution interpolated climate surfaces for global land areas. Int. J. Climatol. 25, 1965–1978 (2005).
2.	Zomer, R. J., Trabucco, A., Van Straaten, O. & Bossio, D. A. Carbon, Land and Water: A Global Analysis of the Hydrologic Dimensions of Climate Change Mitigation through Afforestation/Reforestation. International Water Management Institute 101, (2006).
3.	Trabucco, A. & Zomer, R. Global Aridity Index (Global-Aridity) and Global Potential Evapo-Transpiration (Global-PET) Geospatial Database. (2009).
4.	Tuanmu, M. N. & Jetz, W. A global 1-km consensus land-cover product for biodiversity and ecosystem modelling. Glob. Ecol. Biogeogr. 23, 1031–1045 (2014).

Species distribution models can be downloaded from: http://www.conabio.gob.mx/informacion/gis/) under the category: Biodiversidad > Agrobiodiversidad y agroecosistemas > Parientes Silvestres de Cultivos

