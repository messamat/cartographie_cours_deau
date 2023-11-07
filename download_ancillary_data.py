import os
import codecs
from setup_classement import  *
out_dir = os.path.join(datdir, 'données_auxiliaires')
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

#Download administrative boundaries
admin_express_url = "https://wxs.ign.fr/x02uy2aiwjo9bm8ce5plwqmr/telechargement/prepackage/ADMINEXPRESS_SHP_WGS84G_PACK_2023-05-04$ADMIN-EXPRESS_3-2__SHP_WGS84G_FRA_2023-05-03/file/ADMIN-EXPRESS_3-2__SHP_WGS84G_FRA_2023-05-03.7z"
admin_express_zip_path = os.path.join(out_dir, 'ADMIN-EXPRESS_3-2__SHP_WGS84G_FRA_2023-05-03.7z')
admin_express_unzipped_path = os.path.join(out_dir, 'ADMIN-EXPRESS_3-2__SHP_WGS84G_FRA_2023-05-03')
if not (os.path.exists(admin_express_unzipped_path) or os.path.exists(admin_express_zip_path)):
    with open(admin_express_zip_path, "wb") as file:
        # get request
        print(f"Downloading {Path(admin_express_url).name}")
        response = requests.get(admin_express_url, verify=False)
        file.write(response.content)

################  Download hydrographic basis for classification  #########################################################
# According to Note 2015 "Appui Onema et IGN à l’inventaire des cours d’eau police de l’eau –Volet information géographique"

#Download Basins BD Topage 2023 ----------------------------------------------------------------------------------------
#From https://www.sandre.eaufrance.fr/atlas/srv/fre/catalog.search#/metadata/5ffa1e5b-4e33-4365-91f6-3d1011e466d6
topage_dir = os.path.join(out_dir, 'topage')
if not os.path.exists(topage_dir):
    os.mkdir(topage_dir)

topage_bv_url = ("https://services.sandre.eaufrance.fr/telechargement/geo/ETH/BDTopage/2023/BassinVersantTopographique/"
                 "BassinVersantTopographique_FXX-gpkg.zip")
out_topage_bv_zip_path = os.path.join(topage_dir, os.path.split(topage_bv_url)[1])
out_topage_bv_unzipped_path = os.path.splitext(out_topage_bv_zip_path)[0]
if not (os.path.exists(out_topage_bv_zip_path) or os.path.exists(out_topage_bv_unzipped_path)):
    with open(out_topage_bv_zip_path, "wb") as file:
        # get request
        print(f"Downloading {Path(topage_bv_url).name}")
        response = requests.get(topage_bv_url, verify=False)
        file.write(response.content)

#Download Carthage 2014 ------------------------------------------------------------------------------------------------
"""Download from http://services.sandre.eaufrance.fr/telechargement/geo/ETH/BDCarthage/FXX/2014/"""
carthage_dir = os.path.join(out_dir, 'carthage')
if not os.path.exists(carthage_dir):
    os.mkdir(carthage_dir)

carthage_url = ("http://services.sandre.eaufrance.fr/telechargement/geo/ETH/BDCarthage/FXX/2014/arcgis/FranceEntiere/"
                   "COURS_D_EAU_FXX-shp.zip")
out_carthage_zip_path = os.path.join(carthage_dir, os.path.split(carthage_url)[1])
out_carthage_unzipped_path = os.path.splitext(out_carthage_zip_path)[0]
if not (os.path.exists(out_carthage_zip_path) or os.path.exists(out_carthage_unzipped_path)):
    with open(out_carthage_zip_path, "wb") as file:
        # get request
        print(f"Downloading {Path(carthage_url).name}")
        response = requests.get(carthage_url, verify=False)
        file.write(response.content)


#Download BD topo 2015 (151 edition) -----------------------------------------------------------------------------------
"""Download from https://geoservices.ign.fr/bdtopo#telechargement2015"""
bdtopo_dir = os.path.join(out_dir, 'BDTOPO')
if not os.path.exists(bdtopo_dir):
    os.mkdir(bdtopo_dir)

for i in range(1,96):
    if i != 20:
        bdtopo_url = ("https://wxs.ign.fr/859x8t863h6a09o9o6fy4v60/telechargement/inspire/BDTOPO_TOUSTHEMES_FR151_2015"
             "-03-26$BDTOPO_2-1_TOUSTHEMES_SHP_LAMB93_D0{0}_2015-03-26/file/BDTOPO_2-1_TOUSTHEMES_SHP_LAMB93_"
             "D0{0}_2015-03-26.7z").format(str(i).zfill(2))
        out_bdtopo_zip_path = os.path.join(bdtopo_dir, os.path.split(bdtopo_url)[1])
        out_bdtopo_unzipped_path = os.path.splitext(out_bdtopo_zip_path)[0]
        if not (os.path.exists(out_bdtopo_zip_path) or os.path.exists(out_bdtopo_unzipped_path)):
            with open(out_bdtopo_zip_path, "wb") as file:
                # get request
                print(f"Downloading {Path(bdtopo_url).name}")
                response = requests.get(bdtopo_url, verify=False)
                file.write(response.content)

#Download BCAE 2023 ----------------------------------------------------------------------------------------------------
#https://geoservices.ign.fr/bcae
bcae_dir = os.path.join(out_dir, 'BCAE_20231106')
if not os.path.exists(bcae_dir):
    os.mkdir(bcae_dir)

bcae_url = ('https://wxs.ign.fr/zhy6cayp63shfo3dzoxz5lk7/telechargement/prepackage/BCAE_PACK_2023$BCAE4_SHP_RGF93LAMB93_'
            'FXX_2023-01-01/file/BCAE4_SHP_RGF93LAMB93_FXX_2023-01-01.7z')
out_bcae_zip_path = os.path.join(bcae_dir, os.path.split(bcae_url)[1])
out_bcae_unzipped_path = os.path.splitext(out_bcae_zip_path)[0]
if not (os.path.exists(out_bcae_zip_path) or os.path.exists(out_bcae_unzipped_path)):
    with open(out_bcae_zip_path, "wb") as file:
        # get request
        print(f"Downloading {Path(bcae_url).name}")
        response = requests.get(bcae_url, verify=False)
        file.write(response.content)

#Download RHT - from Hervé Pella ---------------------------------------------------------------------------------------

################  Download point-based data  ###########################################################################
#Download ONDE ----------------------------------------------------------------------------------------------------
"https://onde.eaufrance.fr/content/t%C3%A9l%C3%A9charger-les-donn%C3%A9es-des-campagnes-par-ann%C3%A9e"
onde_dir = os.path.join(out_dir, 'onde')
if not os.path.exists(onde_dir):
    os.mkdir(onde_dir)

for i in range(2012, 2024):
    onde_url = ("https://onde.eaufrance.fr/sites/default/files/fichiers-telechargeables/onde_france_{}.zip").format(i)
    out_onde_zip_path = os.path.join(onde_dir, os.path.split(onde_url)[1])
    out_onde_unzipped_path = os.path.splitext(out_onde_zip_path)[0]
    if not (os.path.exists(out_onde_zip_path) or os.path.exists(out_onde_unzipped_path)):
        with open(out_onde_zip_path, "wb") as file:
            # get request
            print(f"Downloading {Path(onde_url).name}")
            response = requests.get(onde_url, verify=False)
            file.write(response.content)

#Download hydrobiology data on Naiades  --------------------------------------------------------------------------------
#https://naiades.eaufrance.fr/france-entiere#/
#naiades
naiade_dir = os.path.join(out_dir, 'naiade')
if not os.path.exists(naiade_dir):
    os.mkdir(naiade_dir)

naiade_hydrobio_url = "https://naiades.eaufrance.fr/reports/reportsperyear/HB/Naiades_Export_France_Entiere_HB.zip"
out_naiade_hydrobio_zip_path = os.path.join(naiade_dir, os.path.split(naiade_hydrobio_url)[1])
out_naiade_hydrobio_unzipped_path = os.path.splitext(out_naiade_hydrobio_zip_path)[0]
if not (os.path.exists(out_naiade_hydrobio_zip_path) or os.path.exists(out_naiade_hydrobio_unzipped_path)):
    with open(out_naiade_hydrobio_zip_path, "wb") as file:
        # get request
        print(f"Downloading {Path(naiade_hydrobio_url).name}")
        response = requests.get(naiade_hydrobio_url, verify=False)
        file.write(response.content)

#Download Aspe fish data on hubeau API  --------------------------------------------------------------------------------

aspe_dir = os.path.join(out_dir, 'aspe')
if not os.path.exists(aspe_dir):
    os.mkdir(aspe_dir)

out_aspe_unzipped_path = os.path.join(aspe_dir, 'aspe_etat_piscicole_observations.csv')

api_poisson_code = 206
page_poisson = 1
while api_poisson_code==206:
    api_poissons = requests.get("https://hubeau.eaufrance.fr/api/v1/etat_piscicole/observations.csv?"
                                "bbox=-10.0&bbox=35.0&bbox=15.0&bbox=50&page={}&size=20000".format(page_poisson))
    out_aspe_page = "{0}_20231106_{1}.csv".format(os.path.splitext(out_aspe_unzipped_path)[0],
                                              page_poisson)
    print("Saving {}...".format(os.path.split(out_aspe_page)[1]))
    with codecs.open(out_aspe_page, "w", encoding="utf-8") as f:
        f.write(api_poissons.text)
        f.close()
    api_poisson_code=api_poissons.status_code
    page_poisson += 1

tab_pages_list = Path(aspe_dir).glob("{}_20231106_*_.csv".format(os.path.splitext(out_aspe_unzipped_path)[0]))

################  Download potential drivers to relate to  #############################################################
#Crop and Agricultural parcels (to get boundaries and crops for since 2015)
#https://geoservices.ign.fr/rpg
"https://wxs.ign.fr/0zf5kvnyfgyss0dk5dvvq9n7/telechargement/prepackage/RPG_PACK_DIFF_FXX_2022_01$RPG_2-0__GPKG_LAMB93_FXX_2022-01-01/file/RPG_2-0__GPKG_LAMB93_FXX_2022-01-01.7z""

#---- BDHaie ------
"""https://geoservices.ign.fr/bdhaie"""

#---- Land cover ------------
#BD foret
"""https://geoservices.ign.fr/bdforet"""
#CES Occupation des sols (Inglada et al. 2017)
"https://www.theia-land.fr/product/carte-doccupation-des-sols-de-la-france-metropolitaine/"


#---- Topography (DEM) -------------
#Average slope. Topographic wetness index. Curvature. Upstream area
"""https://geoservices.ign.fr/rgealti#telechargement5m"""

#---- Densité de population -------------
#Revenus, pauvreté et niveau de vie en 2019 - Données carroyées (middle year of 2015-2023)
"""https://www.insee.fr/fr/statistiques/7655503?sommaire=7655515"""

#---- Irrigation par commune -------------
#Recensement agricole par commune 2020: https://stats.agriculture.gouv.fr/cartostat/#bbox=252869,6953760,876913,518649&c=indicator&i=cult_2020_2.partirrig20_&selcodgeo=17336&t=A02&view=map11


#---- Prélèvement d'eau ---------
"""https://hubeau.eaufrance.fr/page/api-prelevements-eau"""


#---- Upstream glaciers -------------

#---- Sols #----
#Propriétés de granulométrie des sols
"""https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.57745/N4E4NE"""

#Réservoir en eau utile des sols
"""https://www.theia-land.fr/product/carte-du-reservoir-en-eau-utile-des-sols-de-france-metropolitaine/"""

#Dominant soil type
#1/250,000 national map is not available: https://artificialisation.developpement-durable.gouv.fr/bases-donnees/referentiel-regional-pedologique
#1/10001000
"""https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/BPN57S"""

#---- Lithology ------------
#BD Charm-50
"""https://www.geocatalogue.fr/Detail.do?fileIdentifier=94636790-8615-11dc-9e02-0050568151b7"""
#S_FGEOL Formations géologique polygons

#Intermittency - Snelder

#Intermittency - ONDE

#Climate - AURHELY?
# Write to JP Vidal or Eric Sauquet

#Barriers
#Amber Atlas#
"""https://amber.international/european-barrier-atlas/"""

#Mean bed substrate size - RHT

############### Will not get ##########################
#Soil hydraulic conductivity
#Dissection density is generally inversely related to the hydraulic conductivity of the underlying soil [Montgomery and Dietrich, 1989].

#CarHab is not available nationally yet
#1/250000 soil map either
#High-resolution AI land cover either
