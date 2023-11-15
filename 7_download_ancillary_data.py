import os
import codecs
import bs4 as BeautifulSoup
#import py7zr
import zipfile

from setup_classement import  *

out_dir = os.path.join(datdir, 'données_auxiliaires')
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

def standard_download_zip(in_url, out_rootdir, out_name):
    download_dir = os.path.join(out_rootdir, out_name)
    if not os.path.exists(download_dir):
        os.mkdir(download_dir)

    zip_path = os.path.join(download_dir, os.path.split(in_url)[1])
    unzipped_path = os.path.splitext(zip_path)[0]
    if not (os.path.exists(zip_path) or os.path.exists(unzipped_path)):
        print(f"Downloading {Path(in_url).name}")
        response = requests.get(in_url, verify=False)
        with open(zip_path, "wb") as file:
            # get request
            file.write(response.content)
    else:
        print("{} already exists. Skipping...".format(unzipped_path))

#Download administrative boundaries
standard_download_zip(in_url=("https://wxs.ign.fr/x02uy2aiwjo9bm8ce5plwqmr/telechargement/prepackage/"
                              "ADMINEXPRESS_SHP_TERRITOIRES_PACK_2023-10-16$ADMIN-EXPRESS_3-2__SHP_LAMB93_FXX_2023-10-16/f"
                              "ile/ADMIN-EXPRESS_3-2__SHP_LAMB93_FXX_2023-10-16.7z"),
                      out_rootdir=out_dir,
                      out_name="admin_express")

################  Download hydrographic basis for classification  #########################################################
# According to Note 2015 "Appui Onema et IGN à l’inventaire des cours d’eau police de l’eau –Volet information géographique"

#Download Basins BD Topage 2023 ----------------------------------------------------------------------------------------
#From https://www.sandre.eaufrance.fr/atlas/srv/fre/catalog.search#/metadata/5ffa1e5b-4e33-4365-91f6-3d1011e466d6
standard_download_zip(in_url=("https://services.sandre.eaufrance.fr/telechargement/geo/ETH/BDTopage/2023/"
                              "BassinVersantTopographique/BassinVersantTopographique_FXX-gpkg.zip"),
                      out_rootdir=out_dir,
                      out_name="topage")

#Download HydroBASINS
standard_download_zip(in_url=("https://data.hydrosheds.org/file/hydrobasins/customized_with_lakes/"
                              "hybas_lake_eu_lev01-12_v1c.zip"),
                      out_rootdir=out_dir,
                      out_name="hydrosheds")

#Download Carthage 2014 ------------------------------------------------------------------------------------------------
"""Download from http://services.sandre.eaufrance.fr/telechargement/geo/ETH/BDCarthage/FXX/2014/"""
standard_download_zip(
    in_url=("http://services.sandre.eaufrance.fr/telechargement/geo/ETH/BDCarthage/FXX/2014/arcgis/FranceEntiere/"
                   "COURS_D_EAU_FXX-shp.zip"),
    out_rootdir=out_dir,
    out_name="carthage")

standard_download_zip(
    in_url=("http://services.sandre.eaufrance.fr/telechargement/geo/ETH/BDCarthage/FXX/2014/arcgis/FranceEntiere/"
            "TRONCON_HYDROGRAPHIQUE_FXX-shp.zip"),
    out_rootdir=out_dir,
    out_name="carthage")

#Download BD topo 2015 (151 edition) -----------------------------------------------------------------------------------
"""Download from https://geoservices.ign.fr/bdtopo#telechargement2015"""

for i in range(1,96):
    if i != 20:
        standard_download_zip(
            in_url=("https://wxs.ign.fr/859x8t863h6a09o9o6fy4v60/telechargement/inspire/BDTOPO_TOUSTHEMES_FR151_2015"
             "-03-26$BDTOPO_2-1_TOUSTHEMES_SHP_LAMB93_D0{0}_2015-03-26/file/BDTOPO_2-1_TOUSTHEMES_SHP_LAMB93_"
             "D0{0}_2015-03-26.7z").format(str(i).zfill(2)),
            out_rootdir=out_dir,
            out_name="bdtopo151")

#Download BCAE 2023 ----------------------------------------------------------------------------------------------------
#https://geoservices.ign.fr/bcae
standard_download_zip(
    in_url=('https://wxs.ign.fr/zhy6cayp63shfo3dzoxz5lk7/telechargement/prepackage/BCAE_PACK_2023$BCAE4_SHP_RGF93LAMB93_'
            'FXX_2023-01-01/file/BCAE4_SHP_RGF93LAMB93_FXX_2023-01-01.7z'),
    out_rootdir=out_dir,
    out_name="BCAE_20231106")

#Download RHT - from Hervé Pella ---------------------------------------------------------------------------------------

################  Download point-based data  ###########################################################################
#Download ONDE ----------------------------------------------------------------------------------------------------
"https://onde.eaufrance.fr/content/t%C3%A9l%C3%A9charger-les-donn%C3%A9es-des-campagnes-par-ann%C3%A9e"

for i in range(2012, 2024):
    standard_download_zip(
        in_url=("https://onde.eaufrance.fr/sites/default/files/fichiers-telechargeables/onde_france_{}.zip").format(i),
        out_rootdir=out_dir,
        out_name="onde")

#Download hydrobiology data on Naiades  --------------------------------------------------------------------------------
#https://naiades.eaufrance.fr/france-entiere#/
#naiades
standard_download_zip(
    in_url=("https://naiades.eaufrance.fr/reports/reportsperyear/HB/Naiades_Export_France_Entiere_HB.zip"),
    out_rootdir=out_dir,
    out_name="naiade")

#Download Aspe fish data on hubeau API  --------------------------------------------------------------------------------
"""Irz P, Vigneron T, Poulet N, Cosson E, Point T, Baglinière E, Porcher J-P. 2022. A long-term monitoring database on 
fish and crayfish species in French rivers. Knowl. Manag. Aquat. Ecosyst., 423, 25. (DOI: 10.1051/kmae/2022021)
 Repository: https://zenodo.org/records/7099129"""
standard_download_zip(
    in_url=("https://zenodo.org/records/7099129/files/aspe_data_paper_accepted.zip"),
    out_rootdir=out_dir,
    out_name="aspe")

#Download Global ETOPO  --------------------------------------------------------------------------------
#https://fred.igb-berlin.de/data/package/837
#https://fred.igb-berlin.de/data/file_download/1692

################  Download potential drivers to relate to  #############################################################
#Crop and Agricultural parcels (to get boundaries and crops for since 2015)----------------------------------------------------
#https://geoservices.ign.fr/rpg
standard_download_zip(
    in_url=("https://wxs.ign.fr/0zf5kvnyfgyss0dk5dvvq9n7/telechargement/prepackage/RPG_PACK_DIFF_FXX_2022_01$RPG_2-0__"
            "GPKG_LAMB93_FXX_2022-01-01/file/RPG_2-0__GPKG_LAMB93_FXX_2022-01-01.7z"),
    out_rootdir=out_dir,
    out_name="ign_rpg")

#---- BDHaie (from BDTopo) -------------------------------------------------------------------------------------------------
"""https://geoservices.ign.fr/bdhaie"""
for region in [11, 24, 27, 28, 32, 44, 52, 53, 75, 76, 84, 93]:
    standard_download_zip(
        in_url=("https://wxs.ign.fr/859x8t863h6a09o9o6fy4v60/telechargement/prepackage/"
                "BDTOPOV3-TOUSTHEMES-REGION_GPKG_PACK_233$BDTOPO_3-3_TOUSTHEMES_GPKG_LAMB93_R{0}_2023-09-15/file/"
                "BDTOPO_3-3_TOUSTHEMES_GPKG_LAMB93_R{0}_2023-09-15.7z").format(region),
        out_rootdir=out_dir,
        out_name="bdtopo233")

#Download BD topo 2019 (pack 191) for batiments (to spread population) ------------------------------------------------------------
"""Download from https://geoservices.ign.fr/bdtopo#telechargement2015"""
for i in range(1,96):
    if i != 20:
        standard_download_zip(
            in_url=("https://wxs.ign.fr/859x8t863h6a09o9o6fy4v60/telechargement/prepackage/"
                    "BDTOPOV3-TOUSTHEMES-DEPARTEMENT-PACK_191$BDTOPO_3-0_TOUSTHEMES_SHP_LAMB93_D0{0}_2019-03-15/file/"
                    "BDTOPO_3-0_TOUSTHEMES_SHP_LAMB93_D0{0}_2019-03-15.7z").format(str(i).zfill(2)),
            out_rootdir=out_dir,
            out_name="bdtopo191")

#---- Land cover -------------------------------------------------------------------------------------------------------
#BD foret---------------------------------------------------------------------------------------------------------------------
"""https://geoservices.ign.fr/bdforet"""

in_url = "https://geoservices.ign.fr/bdforet"
session = requests.Session()
# ... whatever other requests config you need here
response = session.get(in_url, verify=False)
soup = BeautifulSoup.BeautifulSoup(response.text, "html.parser")
corsica_urls = ["https://wxs.ign.fr/24pq1j404fo3y3nhf6rc9i1z/telechargement/inspire/BDFORETV2-PACK_14-09-2022$BDFORET_2-0__SHP_LAMB93_{}".format(suffix) for suffix in
                 ["D02A_2017-05-10/file/BDFORET_2-0__SHP_LAMB93_D02A_2017-05-10.7z",
                  "D02B_2016-02-16/file/BDFORET_2-0__SHP_LAMB93_D02B_2016-02-16.7z"]
                ]
for res in soup.findAll('a'):  # images, css, etc..
    if res.has_attr('href'):  # check inner tag (file object) MUST exists
        if re.search("BDFORETV2[-]PACK_14[-]09[-]2022.*[.]7z$", res['href']):
            print(res['href'])
            if not res['href'] in corsica_urls:
                standard_download_zip(
                    in_url=res['href'],
                    out_rootdir=out_dir,
                    out_name="bdforet_v2")
del session, response, soup

#CES Occupation des sols (Inglada et al. 2017)------------------------------------------------------------------------------
"https://www.theia-land.fr/product/carte-doccupation-des-sols-de-la-france-metropolitaine/"
oso_dir = os.path.join(out_dir, 'oso')
if not os.path.exists(oso_dir):
    os.mkdir(oso_dir)
for year,repo in {#2018:"https://zenodo.org/records/3613415/files/Classif_Seed_0.tif", 2018 layer is buggy between 47N and 49N. Cannot use it
                  2019:"https://zenodo.org/records/6538321/files/Classif_Seed_0.tif",
                  2020:"https://zenodo.org/records/6538861/files/Classif_Seed_0_2020.tif",
                  2021:"https://zenodo.org/records/6538910/files/Classif_Seed_0_2021.tif"}.items():
    standard_download_zip(
        in_url=repo,
        out_rootdir=oso_dir,
        out_name='oso_{}'.format(year))

#---- Climate ---------------------------------------------------------------------------------------------------------
#Global aridity index
'''Download Global Aridity Index and Potential Evapotranspiration ET0 Climate Database v2'''
standard_download_zip(
    in_url="https://figshare.com/ndownloader/files/34377269",
    out_rootdir=out_dir,
    out_name='gaiv3')

#---- Topography (DEM) --------------------------------------------------------------------------------------------------------
#Average slope. Topographic wetness index. Curvature. Upstream area
"""https://geoservices.ign.fr/bdalti"""

session = requests.Session()
# ... whatever other requests config you need here
response = session.get("https://geoservices.ign.fr/bdalti", verify=False)
soup = BeautifulSoup.BeautifulSoup(response.text, "html.parser")
for res in soup.findAll('a'):  # images, css, etc..
    if res.has_attr('href'):  # check inner tag (file object) MUST exists
        if re.search("BDALTI-25M_PACK_FXX_.*[.]7z$", res['href']):
            standard_download_zip(
                in_url=re.findall("https.*[7z]", res['href'])[0],
                out_rootdir=out_dir,
                out_name="bdalti")
del session, response, soup


#---- Densité de population -------------------------------------------------------------------------------------------
#Revenus, pauvreté et niveau de vie en 2019 - Données carroyées (middle year of 2015-2023)
"""https://www.insee.fr/fr/statistiques/7655503?sommaire=7655515"""
standard_download_zip(
    in_url="https://www.insee.fr/fr/statistiques/fichier/7655503/Filosofi2019_carreaux_nivNaturel_shp.zip",
    out_rootdir=out_dir,
    out_name="insee_pop")

standard_download_zip(
    in_url="https://www.insee.fr/fr/statistiques/fichier/7655475/Filosofi2019_carreaux_200m_shp.zip",
    out_rootdir=out_dir,
    out_name="insee_pop")


#---- Irrigation par commune -------------------------------------------------------------------------------------------
#Recensement agricole par commune 2020: https://stats.agriculture.gouv.fr/cartostat/#bbox=252869,6953760,876913,518649&c=indicator&i=cult_2020_2.partirrig20_&selcodgeo=17336&t=A02&view=map11
#Manually downloaded in /agreste

#---- Prélèvement d'eau ----------------------------------------------------------------------------------------------------
"""https://hubeau.eaufrance.fr/page/api-prelevements-eau"""
bnpe_dir = os.path.join(out_dir, 'bnpe')
if not os.path.exists(bnpe_dir):
    os.mkdir(bnpe_dir)

#Download ouvrages
out_bnpe_unzipped_path = os.path.join(bnpe_dir, 'bnpe_ouvrages.csv')

if not os.path.exists(out_bnpe_unzipped_path):
    api_ouvrage_code = 206
    page_ouvrage = 1

    api_ouvrage_count = requests.get("https://hubeau.eaufrance.fr/api/v1/prelevements/referentiel/ouvrages?&size=2").json()["count"]

    while api_ouvrage_code==206:
        out_bnpe_page = "{0}_20231107_{1}.csv".format(os.path.splitext(out_bnpe_unzipped_path)[0], page_ouvrage)

        if not os.path.exists(out_bnpe_page):
            api_ouvrage = requests.get("https://hubeau.eaufrance.fr/api/v1/prelevements/referentiel/ouvrages.csv?"
                                        "bbox=-10.0&bbox=35.0&bbox=15.0&bbox=50&page={}".format(page_ouvrage))
            print("Requesting page {}".format(page_ouvrage))
            print("Saving {}...".format(os.path.split(out_bnpe_page)[1]))
            with open(out_bnpe_page, "w", encoding="utf-8") as f:
                    f.write(api_ouvrage.text)
                    f.close()
            api_ouvrage_code = api_ouvrage.status_code
        else:
            print("{} already exists. Skipping...".format(os.path.split(out_bnpe_page)[1]))
        page_ouvrage += 1

    tab_pages_list = list(Path(bnpe_dir).glob("{}_20231107_*.csv".format(
        os.path.splitext(os.path.split(out_bnpe_unzipped_path)[1])[0])))
    bnpe_pdlist = []
    for filename in tab_pages_list:
        print('Reading {}...'.format(filename))
        df_bnpe = pd.read_csv(filename, encoding='utf-8', sep=";")
        bnpe_pdlist.append(df_bnpe)
    bnpe_pdall = pd.concat(bnpe_pdlist, axis=0, ignore_index=True)
    bnpe_pdall.to_csv(out_bnpe_unzipped_path)
    for filename in tab_pages_list:
        os.remove(filename)

#Download chroniques
out_bnpe_chroniques_unzipped_path = os.path.join(bnpe_dir, 'bnpe_chroniques.csv')

if not os.path.exists(out_bnpe_chroniques_unzipped_path):
    api_chroniques_code = 206
    page_chroniques = 1

    api_chroniques_count = requests.get("https://hubeau.eaufrance.fr/api/v1/prelevements/chroniques?&size=2").json()["count"]

    while api_chroniques_code==206:
        out_bnpe_page = "{0}_20231107_{1}.csv".format(os.path.splitext(out_bnpe_chroniques_unzipped_path)[0],
                                                      page_chroniques)

        if not os.path.exists(out_bnpe_page):
            api_chroniques = requests.get("https://hubeau.eaufrance.fr/api/v1/prelevements/chroniques.csv?"
                                        "bbox=-10.0&bbox=35.0&bbox=15.0&bbox=50&page={}".format(page_chroniques))
            print("Requesting page {}".format(page_chroniques))
            print("Saving {}...".format(os.path.split(out_bnpe_page)[1]))
            with open(out_bnpe_page, "w", encoding="utf-8") as f:
                    f.write(api_chroniques.text)
                    f.close()
            api_chroniques_code = api_chroniques.status_code
        else:
            print("{} already exists. Skipping...".format(os.path.split(out_bnpe_page)[1]))
        page_chroniques += 1

    tab_pages_list = list(Path(bnpe_dir).glob("{}_20231107_*.csv".format(
        os.path.splitext(os.path.split(out_bnpe_chroniques_unzipped_path)[1])[0])))
    bnpe_pdlist = []
    for filename in tab_pages_list:
        print('Reading {}...'.format(filename))
        df_bnpe = pd.read_csv(filename, encoding='utf-8', sep=";")
        bnpe_pdlist.append(df_bnpe)
    bnpe_pdall = pd.concat(bnpe_pdlist, axis=0, ignore_index=True)
    bnpe_pdall.to_csv(out_bnpe_chroniques_unzipped_path)
    for filename in tab_pages_list:
        os.remove(filename)

#Can also download from Sandre
"https://www.sandre.eaufrance.fr/atlas/srv/eng/catalog.search#/metadata/2a6fdd18-39ef-4d8a-8e84-84ef86e4a159"
"https://services.sandre.eaufrance.fr/telechargement/geo/PRL/2022/SHP/OuvragePrel_FRA-shp.zip"

#---- Sols #----
#Propriétés de granulométrie des sols - GlobalSoilMap-------------------------------------------------------------------
#Download manually from
"""https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.57745/N4E4NE"""

#Réservoir en eau utile des sols----------------------------------------------------------------------------------------
"""https://www.theia-land.fr/product/carte-du-reservoir-en-eau-utile-des-sols-de-france-metropolitaine/"""
#Download manually from
"https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/9IRARJ"

#Dominant soil type-----------------------------------------------------------------------------------------------------
#1/250,000 national map is not available: https://artificialisation.developpement-durable.gouv.fr/bases-donnees/referentiel-regional-pedologique
#1/10001000
#Download manually from
"""https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/BPN57S"""

#---- Lithology --------------------------------------------------------------------------------------------------------
#BD Charm-50
"""https://www.geocatalogue.fr/Detail.do?fileIdentifier=94636790-8615-11dc-9e02-0050568151b7"""
#S_FGEOL Formations géologique polygons

for i in range(1,96):
    if i != 20:
        standard_download_zip(
            in_url=("http://infoterre.brgm.fr/telechargements/BDCharm50/GEO050K_HARM_0{}.zip").format(str(i).zfill(2)),
            out_rootdir=out_dir,
            out_name="bdcharm50")
standard_download_zip(
    in_url="http://infoterre.brgm.fr/telechargements/BDCharm50/GEO050K_HARM_075_077_078_091_092_093_094_095.zip",
    out_rootdir=out_dir,
    out_name="bdcharm50")

#Intermittency - Snelder -----------------------------------------------------------------------------------------------
#Got data from Hervé Pella

#Intermittency - ONDE  -------------------------------------------------------------------------------------------------
#Will post-process

#Climate - AURHELY?  ----------------------------------------------------------------------------------------------------
# Write to JP Vidal or Eric Sauquet

#Barriers - Amber Atlas# ----------------------------------------------------------------------------------------------------
"""https://amber.international/european-barrier-atlas/"""
#Download from Figshare:
"""https://figshare.com/articles/dataset/AMBER_Atlas_of_Instream_Barriers_in_Europe/12629051"""
amber_url = "https://figshare.com/ndownloader/articles/12629051/versions/5"
standard_download_zip(
    in_url=amber_url,
    out_rootdir=out_dir,
    out_name="amber")

#Mean bed substrate size - RHT  ----------------------------------------------------------------------------------------
#Got data from Hervé Pella

#################################### UNZIP ANCILLARY DATA ##############################################################
#Not working for now -- will work on it
# files_tounzip = []
# for ext in ['zip', '7z']:
#     files_tounzip.extend(Path(out_dir).rglob('*.{}'.format(ext)))
# for f in files_tounzip:
#     print(f)
#     if f.suffix == '.7z':
#         with py7zr.SevenZipFile(f, mode='r') as z:
#             if not os.path.exists(os.path.split(str(f))[0]):
#                 z.extractall(os.path.split(str(f))[0])
#     elif f.suffix == '.zip':
#         with zipfile.ZipFile(f) as z:
#             if not os.path.exists(os.path.split(str(f))[0]):
#                 z.extractall(path=os.path.split(str(f))[0])
#     else:
#         print("{} is neither a 7z file nor a zip file".format(f))

############ POTENTIAL
#Theia Thisme data: https://thisme.cines.teledetection.fr/home
#SurfWater - 10 m - Water Surface occurence by month and by year (https://gitlab.irstea.fr/loic.lozach/thismescripts/-/blob/master/THISMEDownload.py)
#Plot moisture

############### Will not get ##########################
#Soil hydraulic conductivity
#Dissection density is generally inversely related to the hydraulic conductivity of the underlying soil [Montgomery and Dietrich, 1989].

#CarHab is not available nationally yet
#1/250000 soil map either
#High-resolution AI land cover either

#Theia irrigation is not available nationwide yet
#Theia building footprints are not available nationwide yet



#hub'eau poissons API dysfunctions at 34 with 20000 queries:##############################################################################
# aspe_dir = os.path.join(out_dir, 'aspe')
# if not os.path.exists(aspe_dir):
#     os.mkdir(aspe_dir)
#
# out_aspe_unzipped_path = os.path.join(aspe_dir, 'aspe_etat_piscicole_observations.csv')
#
# if not os.path.exists(out_aspe_unzipped_path):
#     api_poisson_code = 206
#     page_poisson = 1
#     while api_poisson_code==206:
#         out_aspe_page = "{0}_20231106_{1}.csv".format(os.path.splitext(out_aspe_unzipped_path)[0], page_poisson)
#
#         if not os.path.exists(out_aspe_page):
#             api_poisson = requests.get("https://hubeau.eaufrance.fr/api/v1/etat_piscicole/observations.csv?"
#                                         "bbox=-10.0&bbox=35.0&bbox=15.0&bbox=50&page={}&size=19790".format(page_poisson))
#             print("Requesting page {}".format(page_poisson))
#             print("Saving {}...".format(os.path.split(out_aspe_page)[1]))
#             try:
#                 with codecs.open(out_aspe_page, "w", encoding="cp1252") as f:
#                     f.write(api_poisson.text)
#                     f.close()
#             except UnicodeEncodeError:
#                 with codecs.open(out_aspe_page, "w", encoding="utf-8") as f:
#                     f.write(api_poisson.text)
#                     f.close()
#             api_poisson_code = api_poisson.status_code
#         else:
#             print("{} already exists. Skipping...".format(os.path.split(out_aspe_page)[1]))
#         page_poisson += 1
#
#     tab_pages_list = list(Path(aspe_dir).glob("{}_20231106_*.csv".format(
#         os.path.splitext(os.path.split(out_aspe_unzipped_path)[1])[0])))
#     aspe_pdlist = []
#     for filename in tab_pages_list:
#         print('Reading {}...'.format(filename))
#         try:
#             df_aspe = pd.read_csv(filename, encoding='cp1252', sep=";")
#         except UnicodeDecodeError:
#             df_aspe = pd.read_csv(filename, encoding='utf-8', sep=";")
#         aspe_pdlist.append(df_aspe)
#     aspe_pdall = pd.concat(aspe_pdlist, axis=0, ignore_index=True)
#     aspe_pdall.to_csv(out_aspe_unzipped_path)
#
# api_poisson = requests.get("https://hubeau.eaufrance.fr/api/v1/etat_piscicole/observations.csv?"
#                            "bbox=-10.0&bbox=35.0&bbox=15.0&bbox=50&page=3364&size=200".format(page_poisson))
