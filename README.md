# Python code for 'Inconsistent Regulatory Mapping Quietly Threatens Rivers and Streams'
Traduction de ce README en français: [![fr](https://img.shields.io/badge/lang-fr-blue.svg)](https://github.com/messamat/cartographie_cours_deau/blob/master/README.FR.md)  

This repository contains python code associated with _Messager, M. L., Pella, H., & Datry, T. (2024). 
Inconsistent regulatory mapping quietly threatens rivers and streams. Environmental Science & Technology. 
https://doi.org/10.1021/acs.est.4c01859_

A copy of the accepted version (after peer-review) of this article is available under a CC-BY-NC license
at: https://hal.inrae.fr/hal-04700243  
A French translation of the article is also available at: https://hal.inrae.fr/hal-04699448

## Abstract
Even the most stringent environmental law cannot  protect a river if its tributaries remain exposed to pollution and other
threats upstream. Excluding a subset of watercourses from legal protection therefore threatens to alter freshwater ecosystems across 
entire river networks and the services they provide, such as drinking water and flood regulation. 
Considerable attention has been devoted to defining the scope of environmental laws protecting watercourses. 
Yet how these definitions are implemented through regulatory mapping, the cartography of waterbodies that legally
qualify as watercourses and are thus protected, has not been examined outside of the United States.
Here, we demonstrate the consequences of regulatory mapping on the extent of river networks that are protected, 
using France as a case study. By assembling the first map of France’s watercourses protected under the Water Law, 
we estimate that a quarter of previously mapped hydrographic segments were excluded from protection and found stark 
geographical variations in the extent of protected ecosystems. Headwater and nonperennial segments are disproportionately 
excluded by 28% compared to their prevalence (67%) in the overall hydrographic network, with potentially far-reaching 
implications for biodiversity and people. We expect regulatory frameworks in most countries to be equally susceptible to 
local interpretation of legal definitions.

## Introduction

This repository includes the portions of the analysis conducted in Python, which encompass all spatial formatting of the
data prior to data analysis. This analysis workflow needs to be conducted prior to conducting data analysis in R with code in the 
following repository: https://github.com/messamat/cartographie_cours_deau_R. 

These scripts are annotated but could be challenging to follow. If you encounter any trouble, please don't hesitate
to contact Mathis L. Messager for comments and clarifications by email or to log an issue in github.

Files needed to run this analysis are either downloaded directly from the code or were provided by Directions
Departementales des Territoires. Please reach out to Mathis L. Messager for assistance with obtaining these data.

### Prerequisites
All GIS analyses in this study require an ESRI ArcGIS Pro license including the Spatial Analyst extension, 
which itself requires a Windows OS. We used the Python Arcpy module associated with ArcGIS Pro 3.1 in Python 3.9.

## Workflow
In subsequent codes, common acronyms include:
- ce: "cours d'eau", means watercourse
- nce: "non cours d'eau", means non-watercourse
- net: network
- deps: departments
- coms: "communes", means municipality
- ddt: Direction Départementale des Territoires (administrative entity in charge of mapping watercourses at the departmental level)
- cats: catchments ("bassin versant topographique BD Topage")
  
### Utility codes
- [setup_classement.py](https://github.com/messamat/cartographie_cours_deau/blob/master/setup_classement.py) : 
import libraries, define folder structure, and basic utility functions. Called from within other scripts.

### Download and format cartographic data from departments
This part of the workflow is entirely conducted with open-source packages.
- [1_googlesearch_DDTwebsite_bydept.py](https://github.com/messamat/cartographie_cours_deau/blob/master/1_googlesearch_DDTwebsite_bydept.py): 
automatically search google for "cartographie des cours d'eau [department name]" and returns the top 3 results for inspection.
- [2_get_DDThtml_bydept.py](https://github.com/messamat/cartographie_cours_deau/blob/master/2_get_DDThtml_bydept.py): 
download an HTML snapshot of the DDT website for archival

An analysis of all departmental websites and of the inter-ministry data catalogue yielded an extensive metadata table 
(available upon request to the authors). This metadata table constitutes the basis for the subsequent workflow 
downloading and formatting the departmental maps of watercourses. 

- [3_download_data_geoide.py](https://github.com/messamat/cartographie_cours_deau/blob/master/3_download_data_geoide.py): 
download and unzip all departmental maps based on URLs from a metadata table
- [4_get_metadata_lyrs.py](https://github.com/messamat/cartographie_cours_deau/blob/master/4_get_metadata_lyrs.py): 
create a csv table listing all departmental datasets that are available in local data repository, count the number of features classes in each,
and list all attribute names.
- [05_QC_metadata.py](https://github.com/messamat/cartographie_cours_deau/blob/master/5_QC_metadata.py): 
remove invalid geometries, merge Yonne layers (watercourses and non-watercourses),
check that network attribute names and categories match metadata for subsequent manual inspection.
- [6_format_merge_DDTnetworks.py](https://github.com/messamat/cartographie_cours_deau/blob/master/6_format_merge_DDTnetworks.py):
  - Format each departmental map (create editable layer in .gdb format with standard name,
    remove features with incompatible formats or invalid geometry, standardize encoding, standardize attribute names,
    delete duplicate lines keeping the ones with the most information, recategorize watercourses status attribute into five categories:
    cours d'eau [watercourse], non cours d'eau [non-watercourse], indéterminé [uncategorized], inexistant [inexistent], and hors département [outside departmental
    boundaries]), compute number of NAs per column, reproject to a standard coordinate system, convert polygons to centerlines).
  - Merge all harmonized departmental maps
  - Create unique ID (UID).

### Download and process reference networks and socio-environmental variables
This part of the workflow requires an ArcGIS Pro license.

- [7_download_ancillary_data.py](https://github.com/messamat/cartographie_cours_deau/blob/master/7_download_ancillary_data.py):
download administrative and hydrographic boundaries, other (reference) hydrographic networks for comparison (BD TOPO, BD Carthage, BCAE),
download other point-based monitoring data (river drying, macroinvertebrates, and fish data), download socio-environmental drivers.
- [8_preprocess_reference_networks.py](https://github.com/messamat/cartographie_cours_deau/blob/master/8_preprocess_reference_networks.py): 
prepare reference hydrographic data: intersect sub-catchments with department boundaries (this intersection is the finer unit of analysis for this study),
create single feature class for each network, remove duplicate geometries, spatially match BD TOPO and BD Carthage with departmental maps,
intersect all networks (including departmental maps) with catchment-department units of analysis.

- [9_preprocess_ancillary_data.py](https://github.com/messamat/cartographie_cours_deau/blob/master/9_preprocess_ancillary_data.py): 
preprocess all ancillary data (including dasymetric interpolation of irrigated area and population)
and compute summary statistics by catchment-department unit of analysis, and by department
