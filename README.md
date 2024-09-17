# Python code for 'Inconsistent Regulatory Mapping Quietly Threatens Rivers and Streams'

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
The /data folder in the figshare repository contains raw data and the directory structure enables users to reproduce our 
study using the scripts herein.

### Prerequisites
All GIS analyses in this study require an ESRI ArcGIS Pro license including the Spatial Analyst extension, 
which itself requires a Windows OS. We used the Python Arcpy module associated with ArcGIS Pro 3.1 in Python 3.9.

## Workflow
### Utility codes
- [setup_classement.py](https://github.com/messamat/cartographie_cours_deau/blob/master/setup_classement.py) : 
import libraries, define folder structure, and basic utility functions. Called from within other scripts.

### Download and format cartographic data from departments
- [1_googlesearch_DDTwebsite_bydept.py](https://github.com/messamat/cartographie_cours_deau/blob/master/1_googlesearch_DDTwebsite_bydept.py): 
automatically search google for "cartographie des cours d'eau [department name]" and returns the top 3 results for inspection.
- [2_get_DDThtml_bydept.py](https://github.com/messamat/cartographie_cours_deau/blob/master/2_get_DDThtml_bydept.py): 
download 
- [3_download_data_geoide.py](https://github.com/messamat/cartographie_cours_deau/blob/master/3_download_data_geoide.py): 
awra
- [4_get_metadata_lyrs.py](https://github.com/messamat/cartographie_cours_deau/blob/master/4_get_metadata_lyrs.py): 
compute 
- [05_QC_metadata.py](https://github.com/messamat/cartographie_cours_deau/blob/master/5_QC_metadata.py): 
downscale 
- [6_format_merge_DDTnetworks.py](https://github.com/messamat/cartographie_cours_deau/blob/master/6_format_merge_DDTnetworks.py): 
downscale 

### Download and process reference networks and socio-environmental variables
- [7_download_ancillary_data.py](https://github.com/messamat/cartographie_cours_deau/blob/master/7_download_ancillary_data.py): 
compute 
- [8_preprocess_reference_networks.py](https://github.com/messamat/cartographie_cours_deau/blob/master/8_preprocess_reference_networks.py): 
downscale 
- [9_preprocess_ancillary_data.py](https://github.com/messamat/cartographie_cours_deau/blob/master/9_preprocess_ancillary_data.py): 
downscale 
