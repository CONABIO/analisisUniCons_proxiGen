# README

Code for the analyses of the paper "Incorporating evolutionary and threat processes into crop wild relatives conservation" by Wolke Tobón-Niedfeldt, 
Alicia Mastretta-Yanes, Tania Urquiza-Haas, Bárbara Goettsch, Angela P. Cuervo Robayo, Esmeralda Urquiza-Haas, M. Andrea Orjuela-R, Francisca Acevedo Gasman, Oswaldo Oliveros-Galindo, Caroline Burgeff, Diana Rivera-Rodríguez, José de Jesús Sánchez González, Jesús Alarcón-Guerrero, José Sarukhán and Patricia Koleff.


## Data

Data is available at the Dryad repository XXXXX (available upon aceptance). There you could find the data described below. The scripts in `/bin` expect the data organized in the following way:

### Spatial raw rasters

* Rasters (.tif) of species distribution models for each of the cwr taxa uses for analyses should be in the directory:
 `/data/spatial/modelos_darwin_all_final`

Raster names have the taxon name.
 
* Rasters (.tif) of Holdrige live zones should be in the directory:
 `/data/spatial/zv27`

Rasters are named `zv_[1:27].tif`

 * Rasters (.tif) with proxies of genetic diversity (PDG) should be in the directory:
 `/data/spatial/areasProxyDivGen`
 
 Rasters file names are as follow: each of the Holdrige life zones (`zv_[1:27]`) are subdivided with a `_[letter].tif` for each of the PGD for that life zone.


### Zonation 

**Zonation input data** should be in the directory:
`data/spatial/Zonation_input/`

File contect is the following:

* Biodiversity rasters (.tif); result of the combination of species distribution models (SDM, MDP in Spanish) and proxies of genetic diversity (PDG, PGD in Spanish): MDP_PDG. Raster names have taxon name and PGD, according to each of the Holdrige life zones (zv_[1:27]) subdivided with a _[letter].tif for each of the PGD for that life zone.

* Habitat rasters (.asc): condition_habitat2_final. Raster are named Hab_[taxon_name].asc

* Species of Special Interest files (.txt); observation records for taxa without potential distribution model: puntos. Raster names have taxon name.  


**Zonation output rasters** should be in the directory:
`/data/spatial/Zonation_output/`

File content is the following:

* Zonation output with species only (n=116): `01_MDP.rank.compressed.tif"` 
* Zonation output with species + LZ (n=143): `02_MDP_ZV.rank.compressed.tif"`
* Zonation output with species and PDG (n=218): `03_MDP_PDG.rank.compressed.tif`
* Zonation output with species vs PDG (n=5004): `04_MDP_vs_PDG.rank.compressed.tif`
* Zonation output with species and PDG as ADMU (n=116+1): `05_MDP_PDG_ADMU.rank.compressed.tif`
 

**Zonation output curves** should be in the directory::
`/data/spatial/Zonation_final_solutions/`

File content is the following:

* Zonation curves: `E_final_Todos.curves.txt`
* Conservation features: `E_final_Todos.features_info.txt`
* Zonation curves output and features but for taxa exclusively distributing in natural vegetation: 
`E_final_VegPyS.features_info.txt"` and `E_final_VegPyS.curves.txt`.
* Zonation curves output and features but for taxa associated to different habitats, i.e. natural vegetation, agricultural and urban areas: `E_final_HabVarios.features_info.txt` and `E_final_HabVarios.curves.txt`.

### Genetic

Admixture output data for Z. m. parviglumis is in the directory:
`/data/genetic/output/admixture/parvi/`

File content is the following:

* Admixture output K error for CV analysis: `Kerror_parviglumis.txt`
* Admixture Q files for K= 25 and K= 13: `bytaxa_parviglumis.25.Q` and `bytaxa_parviglumis.13.Q"`


### Other
* Table with IUCN category per taxa: 
`/data/spatial/Zonation_final_solutions/IUCN_threat_category.csv`
 


## Non-coding analyses and Figures

Figures 2, 5 and Supplementary Figures 1-3, 5-7, 10 and 11 were made in ArcMap with the spatial data and zonation raster outputs detailed above (georeferenced data, species distribution models, zonation outputs, among others).

## Zonation analyses
Zonation configuration files are in:`/bin/Zonation_files`.

File contents as follows:

* Batch files, contain the command line to call Zonation in the cluster: `.sh` files
* Run settings file, defines the analysis settings: `E_final.dat`
* Biodiversity feature files, list the biodiversity layers: `.spp` files
* Condition file; lists the habitat layers: `habitat_features.txt`
* Groups file; links the biodiversity features to habitat features using habitat groups: `habitat_group_features.txt`
* Species of Special Interest file; lists observation records for species without potential distribution model: `SSI_list.txt`

## Analyses in R and Figures

Code for analyses and plots made in R is available in the `/bin` directory of this repository in the scripts detailed below. For each script there is a .html version of it with an R notebook showing the code and output.


#### 1) `plot_admixture_PGD_teocintles.R`

Explores admixture results to perform the plots and map of: 

* Figure 3. Genetic diversity of Zea mays ssp. parviglumis represented in the proxies of genetic diversity (PGD).

* Supplementary Fig. 4. Potential species distribution model of Z. mays ssp. parviglumis as given by PGD (background colours).

The admixture plot, text and map of Fig. 3 were joined in a single figure using Inkscape.


Data used by this script:


* *Zea mays parviglumis* SDM: `"../data/spatial/modelosDarwinZea/Zea_mays_parviglumis.tif"`
* Proxies of Genetic Diversity croped for *Z. m. parviglumis* SDM:`"../data/spatial/areasProxyDivGen/crop_to_sp_final/PDG_Zea_mays_parviglumis.tif"`
* *Z. m. parviglumis* CV error (admixture output): `"../data/genetic/output/admixture/parvi/Kerror_parviglumis.txt"`
* *Z. m. parviglumis* Q files for K= 25 and K= 13 (admixture output):
`"../data/genetic/output/admixture/parvi/bytaxa_parviglumis.25.Q"` and `"../data/genetic/output/admixture/parvi/bytaxa_parviglumis.13.Q"`.

#### 2) `spatial_analyses_zonationVSproxiesdivgen_all_ms.R`

For each taxon crops the proxies of gen div to the species distribution models. Then performs analyses to compare different ways to incorporate the PDG into the Zonation analyses. 

Data used by this script:

* Zonation outputs and proxies of gen div:
```
"../data/spatial/Zonation_output/comparacion/01_MDP.rank.compressed.tif", # Zonation output with species only (n=116) #1
"../data/spatial/Zonation_output/comparacion/02_MDP_ZV.rank.compressed.tif", # Zonation output with species + LZ (n=143) #2
"../data/spatial/Zonation_output/comparacion/03_MDP_PDG.rank.compressed.tif", # Zonation output with species and PDG (n=218) #3
"../data/spatial/Zonation_output/comparacion/04_MDP_vs_PDG.rank.compressed.tif", # Zonation output with species vs PDG (n=5004) #4
"../data/spatial/Zonation_output/comparacion/05_MDP_PDG_ADMU.rank.compressed.tif", # Zonation output with species and PDG as ADMU (n=116+1) #5
"../data/spatial/areasProxyDivGen/PDG.tif") #Proxies div gen 
```

* Species distribution models of each species: `"../data/spatial/modelos_darwin_all_final/*.tif"`

Output of this script:

* The rasters of proxies of genetic diversity cropped for each taxa: `"../data/spatial/areasProxyDivGen/crop_to_sp/"`.
* Area of proxies of div gen for each zonation solution: `"../data/comparations_output/sol_tidy_spp.txt"`
* Proportion of proxi of div gen for each zolnation solution in relation to its area in the sp distribution: `"../data/comparations_output/sol_prop_spp.txt"`
* Diversity indexes and mean of proportion: `"../data/comparations_output/sols_summary_spp.txt"`

#### 3) `plot_scenarios_PDG.R`

Uses the data `sols_summary_spp.txt` produced by the script `spatial_analyses_zonationVSproxiesdivgen_all_ms.R` to perform the plots of the following figure:

* Fig 4. Performance of five scenarios to represent conservation features of Mesoamerican crop wild relatives, considering 20% of Mexico’s terrestrial area.


#### 4) `plot_Zonation_curves.R`
Performs the analyses and plots of the following figures:

* Fig. 6. Performance curves quantifying the proportion of crop wild relatives within the scenario considering potential species distribution models subdivided by proxies of genetic diversity.

* Supplementary Fig 8. Performance curves quantifying the proportion of taxa
distribution ranges considered for each scenario, grouped by IUCN Red List
Category.

* Supplementary Fig. 9. Performance curves quantifying the proportion of taxa
distribution ranges considering all priority taxa

Data used for this script:

* Read IUCN category per taxa: `"../data/spatial/Zonation_final_solutions/IUCN_threat_category.csv"`
* Zonation curves output and features for all taxa: `"../data/spatial/Zonation_final_solutions/E_final_Todos.curves.txt"` and `"../data/spatial/Zonation_final_solutions/E_final_Todos.features_info.txt"`
* Zonation curves output and features but for taxa exclusively distributing in natural vegetation: `"../data/spatial/Zonation_final_solutions/E_final_VegPyS.features_info.txt"` and `"../data/spatial/Zonation_final_solutions/E_final_VegPyS.curves.txt"`.
* Zonation curves output and features but for taxa associated to different habitats, i.e. natural vegetation, agricultural and urban areas: `"../data/spatial/Zonation_final_solutions/E_final_HabVarios.features_info.txt"` and `"../data/spatial/Zonation_final_solutions/E_final_HabVarios.curves.txt"`.
