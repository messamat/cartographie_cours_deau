# Code Python pour 'Une cartographie réglementaire incohérente menace les rivières et les ruisseaux'
[![en](https://img.shields.io/badge/lang-en-red.svg)](https://github.com/messamat/cartographie_cours_deau/blob/master/README.md)

Ce dépot contient le code python associé avec l'article scientifique _Messager, M. L., Pella, H., & Datry, T. (2024). 
Inconsistent regulatory mapping quietly threatens rivers and streams. Environmental Science & Technology. 
https://doi.org/10.1021/acs.est.4c01859_

Une copie en anglais de la version accepté de l'article (après revue par des pairs) sous licence CC-BY-NC
à l'adresse suivante: https://hal.inrae.fr/hal-04700243  
Une traduction en français est également disponible: https://hal.inrae.fr/hal-04699448

## Résumé
Même la législation environnementale la plus stricte ne peut protéger un cours d'eau si ses affluents restent 
exposés à la pollution et à d'autres menaces en amont. Exclure un sous-ensemble de cours d'eau de la protection
juridique menace donc d'altérer les écosystèmes d'eau douce de réseaux fluviaux entiers et les services qu'ils 
fournissent, tels que l'eau potable et la régulation des crues. Une attention considérable a été accordée à la 
définition du champ d'application des lois environnementales protégeant les cours d'eau. Cependant, la manière
dont ces définitions sont mises en œuvre par le biais de la cartographie réglementaire, c'est-à-dire la cartographie 
des masses d'eau qui sont légalement considérées comme des cours d'eau et donc protégées, n'a pas été étudiée 
en dehors des États-Unis. Nous démontrons ici les conséquences de la cartographie réglementaire sur l'étendue 
des réseaux fluviaux protégés, en utilisant la France comme étude de cas. En assemblant la première carte des 
cours d'eau français protégés au titre de la Loi sur l'eau, nous estimons qu'un quart des segments hydrographiques 
précédemment cartographiés ont été exclus, et constatons de fortes variations géographiques dans l'étendue des 
écosystèmes protégés. Les segments de tête de bassin et les segments non pérennes sont disproportionnellement 
exclus de 28% par rapport à leur prévalence (67 %) dans l'ensemble du réseau hydrographique, avec des implications
potentiellement considérables pour la biodiversité et les populations humaines. Nous nous attendons à ce que les 
cadres réglementaires de la plupart des pays soient également sensibles à l'interprétation locale des définitions 
juridiques. 

## Introduction

Ce dépôt comprend les parties de l'analyse effectuée en Python, qui englobent toute la mise en forme spatiale 
des données avant l'analyse des données. Ce flux de travail doit être effectué avant de procéder à l'analyse 
des données dans R avec le code du dépôt suivant : https://github.com/messamat/cartographie_cours_deau_R.  

Ces scripts sont annotés mais peuvent être difficiles à suivre. Si vous rencontrez des difficultés, 
n'hésitez pas à contacter Mathis L. Messager pour des commentaires et des clarifications par e-mail
ou à signaler un problème sur GitHub.

Les fichiers nécessaires pour exécuter cette analyse sont soit téléchargés directement depuis le code, 
soit fournis par les Directions Départementales des Territoires. Veuillez contacter Mathis L. Messager 
pour obtenir de l'aide concernant ces données. 

### Prérequis
Certains analyses SIG de cette étude nécessitent une licence ESRI ArcGIS Pro incluant l'extension 
Spatial Analyst, qui elle-même nécessite un système d'exploitation Windows. Nous avons utilisé Python 3.9 et
le module Python Arcpy associé à ArcGIS Pro 3.1.

## Flux de travail
Dans les codes suivants, les acronymes courants incluent :
- ce : cours d'eau
- nce : non cours d'eau
- net : réseau
- deps : départements
- coms : communes
- ddt : Direction Départementale des Territoires (entité administrative chargée de la cartographie des cours d'eau au niveau départemental)
- cats : bassin versant topographique BD Topage

### Codes utilitaires
- [setup_classement.py](https://github.com/messamat/cartographie_cours_deau/blob/master/setup_classement.py) : 
importe les bibliothèques, définit la structure des dossiers et les fonctions utilitaires de base. Appelé depuis d'autres scripts.

### Télécharge et formate les données cartographiques des départements
Cette partie du flux de travail est entièrement réalisée avec des packages open-source.
- [1_googlesearch_DDTwebsite_bydept.py](https://github.com/messamat/cartographie_cours_deau/blob/master/1_googlesearch_DDTwebsite_bydept.py) : 
recherche automatiquement sur Google "cartographie des cours d'eau [nom du département]" et renvoie les 3 meilleurs résultats pour inspection.
- [2_get_DDThtml_bydept.py](https://github.com/messamat/cartographie_cours_deau/blob/master/2_get_DDThtml_bydept.py) : 
télécharge une capture HTML du site Web de la DDT pour archivage.

Une analyse de tous les sites Web départementaux et du catalogue de données interministériel a donné lieu 
à une table de métadonnées exhaustive (disponible sur demande aux auteurs). Cette table de métadonnées constitue 
la base du flux de travail suivant pour télécharger et formater les cartes départementales des cours d'eau.

- [3_download_data_geoide.py](https://github.com/messamat/cartographie_cours_deau/blob/master/3_download_data_geoide.py) : 
télécharge et décompresse toutes les cartes départementales basées sur les URL d'une table de métadonnées.
- [4_get_metadata_lyrs.py](https://github.com/messamat/cartographie_cours_deau/blob/master/4_get_metadata_lyrs.py) : 
crée un tableau CSV répertoriant toutes les base de données départementales ayant été téléchargées (automatiquement ou manuellement par envoie d'une DDT),
compte le nombre de géométries dans chaque base de données et répertorie tous les noms d'attributs.
- [05_QC_metadata.py](https://github.com/messamat/cartographie_cours_deau/blob/master/5_QC_metadata.py) : 
supprime les géométries non valides, fusionne les couches de l'Yonne (cours d'eau et non cours d'eau),
vérifie que les noms des attributs de réseau et les catégories correspondent aux métadonnées pour une inspection manuelle ultérieure.
- [6_format_merge_DDTnetworks.py](https://github.com/messamat/cartographie_cours_deau/blob/master/6_format_merge_DDTnetworks.py) :
  - Formatte chaque carte départementale (créer une couche modifiable au format .gdb avec un nom standard,
    supprime les entitées ayant des formats incompatibles ou une géométrie non valide, standardise l'encodage,
    standardise les noms d'attributs, supprime les lignes en double en gardant celles avec le plus d'informations,
    recatégorise l'attribut de statut des cours d'eau en cinq catégories : cours d'eau, non cours d'eau, indéterminé,
    inexistant, et hors département), calcule le nombre de valeurs nulles par colonne, reprojete dans un système de
    coordonnées standard, convertit les polygones en lignes centrales.
  - Fusionne toutes les cartes départementales harmonisées.
  - Crée un identifiant unique (UID) pour chaque ligne.

### Télécharge et traite les réseaux hydrographiques de référence et les variables socio-environnementales
Cette partie du flux de travail nécessite une licence ArcGIS Pro.

- [7_download_ancillary_data.py](https://github.com/messamat/cartographie_cours_deau/blob/master/7_download_ancillary_data.py) : 
télécharge les limites administratives et hydrographiques, d'autres réseaux hydrographiques de référence pour comparaison (BD TOPO, BD Carthage, BCAE),
télécharge d'autres données de surveillance (assèchement des rivières, macroinvertébrés, et données sur les poissons),
 télécharge les facteurs socio-environnementaux.
- [8_preprocess_reference_networks.py](https://github.com/messamat/cartographie_cours_deau/blob/master/8_preprocess_reference_networks.py) : 
prépare les données hydrographiques de référence : intersecte les sous-bassins (bassin versant topographique BD Topage) avec les limites des départements
(cette intersection est l'unité d'analyse la plus fine de cette étude), crée une couche unique pour chaque réseau hydrographique,
supprime les géométries en double, met les géométries des réseaux BD TOPO et BD Carthage en correspondance spatiale avec les cartes départementales,
 intersecte tous les réseaux hydrographiques (y compris les cartes départementales) avec les unités d'analyse des bassins versants-départements.

- [9_preprocess_ancillary_data.py](https://github.com/messamat/cartographie_cours_deau/blob/master/9_preprocess_ancillary_data.py) : 
pré-traite toutes les données auxiliaires (y compris l'interpolation dasymétrique de la surface irriguée et de la population)
et calcule les statistiques récapitulatives par unité d'analyse bassin versant-département, et par département.
