import arcpy

from setup_classement import  *

anci_dir = os.path.join(datdir, "données_auxiliaires")#Ancillary data directory


#------------------------------------------- UNITS OF ANALYSIS ---------------------------------------------------------
#Departments and communes - admin express
admin_dir = os.path.join(anci_dir, "admin_express", "ADMIN-EXPRESS_3-2__SHP_WGS84G_FRA_2023-05-03",
                    "1_DONNEES_LIVRAISON_2023-05-03", "ADE_3-2_SHP_WGS84G_FRA")
deps = os.path.join(admin_dir, "DEPARTEMENT.shp")
coms = os.path.join(admin_dir, "COMMUNE.shp")

#---------------------------------------------- DATA SOURCES -----------------------------------------------------------
#% of irrigated land by commune - agreste
percent_irrigated_crops_communes = os.path.join(anci_dir, "agreste", "data.csv")

#Barrier points - amber
barriers = os.path.join(anci_dir, "amber", "AMBER_BARRIER_ATLAS_V1.csv")

#Fish sampling stations - aspe
fish_stations = os.path.join(anci_dir, "aspe", "raw_data", "raw_data", "csv", "station.csv")

#DEM - BDAlti
bdalti_dir = os.path.join(anci_dir, "bdalti")
bdalti_filelist = getfilelist(bdalti_dir, "^BDALTIV2_25M_FXX_.*[.]asc$")

#Lithology
bdcharm_dir = os.path.join(anci_dir, "bdcharm50")
bdcharm_filelist = getfilelist(bdcharm_dir,"^GEO050K_HARM_[0-9]{3}_S_FGEOL_2154[.]shp$")

#Forest type
bdforet_dir = os.path.join(anci_dir, "bdforet_v2")
bdforet_filelist = getfilelist(bdforet_dir, "FORMATION_VEGETALE.shp")

#Dominant soil type 1:1,000,000
bdgsf = os.path.join(anci_dir, "bdgsf", "30169_L93", "30169_L93.shp")

#BD topo 2023
bdtopo2023_dir = os.path.join(anci_dir, "bdtopo233") #233 is the version/package number (vs e.g. 151 for 2015)
bdtopo2023_gpkglist = getfilelist(bdtopo2023_dir, "BDT_3-3_GPKG_LAMB93_R[0-9]{2}-ED2023-09-15.gpkg$")
bdtopo2023_haie_filelist = [os.path.join("main.haie") for gpkg in bdtopo2023_gpkglist]
"plan_d_eau"
"reservoir"

#Withdrawal points - bnpe
withdrawal_pts=os.path.join(anci_dir, "bnpe", "bnpe_ouvrages.csv")

#GlobalSoilModel - available water capacity ############### TO FINISH
gsm_awc_dir = os.path.join(anci_dir, "gsm")

#Farm parcels - rpg
farm_parcels = os.path.join(anci_dir, "ign_rpg", "RPG_2-0__GPKG_LAMB93_FXX_2022-01-01", "RPG",
                            "1_DONNEES_LIVRAISON_2023-08-01", "RPG_2-0_GPKG_LAMB93_FXX-2022",
                            "PARCELLES_GRAPHIQUES.gpkg", "main.parcelles_graphiques")
#print(arcpy.GetCount_management(farm_parcels)[0])

#Population count - insee
pop_count = os.path.join(anci_dir, "insee_pop", "carreaux_nivNaturel_met.shp")

#Hydrobiology sampling stations - naiade
hydrobio_stations = os.path.join(anci_dir,"naiade", "stations.csv")

#FLow intermittence observation network - ONDE
onde_filelist = getfilelist(os.path.join(anci_dir, "onde"), "onde_france_[0-9]{4}[.]csv$")

#Land cover - theia oso
lc_dir = os.path.join(anci_dir, 'oso')
lc_filedict = {yr: getfilelist(os.path.join(lc_dir, "oso_{}".format(yr)),"Classif_Seed_.*[.]tif$")
               for yr in [2018, 2019, 2020, 2021]}

#Snelder intermittence
snelder_ires = os.path.join(anci_dir, "snelder", "rhtvs2_all_phi_qclass.shp")

#Preprocessing geodatabase
pregdb = os.path.join(resdir, "preprocessing_ancillary_data.gdb")
if not arcpy.Exists(pregdb):
    arcpy.CreateFileGDB_management(out_folder_path=os.path.split(pregdb)[0],
                                   out_name=os.path.split(pregdb)[1])


#--------------------------------- DEM - BDALTI ------------------------------------------------------------------------
bdalti_mosaic = os.path.join(pregdb, "bdalti_25m_mosaic")
#Based on metadata from https://geoservices.ign.fr/documentation/donnees/alti/bdalti
#https://geoservices.ign.fr/sites/default/files/2021-07/DC_BDALTI_2-0.pdf
bdalti_sr = arcpy.SpatialReference(5698)

if not arcpy.Exists(bdalti_mosaic):
    print("Processing {}...".format(bdalti_mosaic))
    arcpy.MosaicToNewRaster_management(input_rasters=bdalti_filelist,
                                       output_location=os.path.split(bdalti_mosaic)[0],
                                       raster_dataset_name_with_extension=os.path.split(bdalti_mosaic)[1],
                                       coordinate_system_for_the_raster=bdalti_sr,
                                       pixel_type="16_BIT_SIGNED",
                                       number_of_bands=1
                                       )

#Compute slope
bdalti_slope = os.path.join(pregdb, "bdalti_25m_slope")
if not arcpy.Exists(bdalti_slope):
    print("Processing {}...".format(bdalti_slope))
    Slope(in_raster=bdalti_mosaic, output_measurement="DEGREE").save(bdalti_slope)

#Compute curvature
bdalti_tcurvature = os.path.join(pregdb, "bdalti_25m_curvature_tangential")
if not arcpy.Exists(bdalti_tcurvature):
    print("Processing {}...".format(bdalti_tcurvature))
    SurfaceParameters(in_raster=bdalti_mosaic,
                      parameter_type="TANGENTIAL_CURVATURE",
                      local_surface_type = "QUADRATIC",
                      z_unit='METER'
                      ).save(bdalti_tcurvature)

bdalti_pcurvature = os.path.join(pregdb, "bdalti_25m_curvature_profile")
if not arcpy.Exists(bdalti_pcurvature):
    print("Processing {}...".format(bdalti_pcurvature))
    SurfaceParameters(in_raster=bdalti_mosaic,
                      parameter_type="PROFILE_CURVATURE",
                      local_surface_type = "QUADRATIC",
                      z_unit='METER'
                      ).save(bdalti_pcurvature)

#Compute topographic wetness index in the future - https://mapscaping.com/topographic-wetness-index-in-arcgis-pro/

#--------------------------------- Lithology - BDCHARM ------------------------------------------------------------------------
{}

bdcharm_mosaic = os.path.join(pregdb, "lithology_bdcharm_mosaic")
bdalti_sr = arcpy.SpatialReference(5698)

if not arcpy.Exists(bcae_fr):
    arcpy.management.Merge(inputs=bcae_filelist, output=bcae_fr)

if not arcpy.Exists(bdalti_mosaic):
    print("Processing {}...".format(bdalti_mosaic))
    arcpy.MosaicToNewRaster_management(input_rasters=bdcharm_filelist,
                                       output_location=os.path.split(bdalti_mosaic)[0],
                                       raster_dataset_name_with_extension=os.path.split(bdalti_mosaic)[1],
                                       coordinate_system_for_the_raster=bdalti_sr,
                                       pixel_type="16_BIT_SIGNED",
                                       number_of_bands=1
                                       )





