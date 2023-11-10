import arcpy.analysis

from setup_classement import  *

anci_dir = os.path.join(datdir, "données_auxiliaires") #Ancillary data directory
pregdb = os.path.join(resdir, "preprocessing_ancillary_data.gdb") #Preprocessing geodatabase

#------------------------------------------- UNITS OF ANALYSIS ---------------------------------------------------------
#Departments and communes - admin express
admin_dir = os.path.join(anci_dir, "admin_express", "ADMIN-EXPRESS_3-2__SHP_LAMB93_FXX_2023-10-16",
                         "ADMIN-EXPRESS", "1_DONNEES_LIVRAISON_2023-10-16", "ADE_3-2_SHP_LAMB93_FXX")
deps = os.path.join(admin_dir, "DEPARTEMENT.shp")
coms = os.path.join(admin_dir, "COMMUNE.shp")

cats_hybasdeps = os.path.join(pregdb, 'BV_hybas0809_depsinters') #BV joined to HydroBASINS and intersected with Departements

#---------------------------------------------- DATA SOURCES -----------------------------------------------------------
#% of irrigated land by commune - agreste
irrigation_communes = os.path.join(anci_dir, "agreste", "data_irrigation_SAU_communes.csv")
irrigation_departements = os.path.join(anci_dir, "agreste", "data_irrigation_SAU_departements.csv")

#Barrier points - amber
barriers = os.path.join(anci_dir, "amber", "AMBER_BARRIER_ATLAS_V1.csv")

#Fish sampling stations - aspe
fish_stations = os.path.join(anci_dir, "aspe", "raw_data", "raw_data", "csv", "station.csv")

#DEM - BDAlti
bdalti_dir = os.path.join(anci_dir, "bdalti")

#Lithology
bdcharm_dir = os.path.join(anci_dir, "bdcharm50")

#Forest type
bdforet_dir = os.path.join(anci_dir, "bdforet_v2")

#Dominant soil type 1:1,000,000
bdgsf = os.path.join(anci_dir, "bdgsf", "30169_L93", "30169_L93.shp")

#BD topo 2023
bdtopo2023_dir = os.path.join(anci_dir, "bdtopo233") #233 is the version/package number (vs e.g. 151 for 2015)
"plan_d_eau"
"reservoir"

#Withdrawal points - bnpe
withdrawal_stations_bnpe=os.path.join(anci_dir, "bnpe", "bnpe_ouvrages.csv")

#GlobalSoilModel - available water capacity ############### TO FINISH
gsm_awc_dir = os.path.join(anci_dir, "gsm")

#Farm parcels - rpg
farm_parcels = os.path.join(anci_dir, "ign_rpg", "RPG_2-0__GPKG_LAMB93_FXX_2022-01-01", "RPG",
                            "1_DONNEES_LIVRAISON_2023-08-01", "RPG_2-0_GPKG_LAMB93_FXX-2022",
                            "PARCELLES_GRAPHIQUES.gpkg", "main.parcelles_graphiques")
#print(arcpy.GetCount_management(farm_parcels)[0])

#Buildings
bdtopo2019_dir = os.path.join(anci_dir, 'bdtopo191')

#Population count - insee
pop_count_variable_mesh = os.path.join(anci_dir, "insee_pop", "carreaux_nivNaturel_met.shp")
pop_count_200m_mesh = os.path.join(anci_dir, "insee_pop", "arreaux_200m_shp")


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

#-------------------------------- OUTPUTS ------------------------------------------------------------------------------
#Scratch gdb
tempgdb = os.path.join(resdir, "scratch.gdb")
if not arcpy.Exists(tempgdb):
    arcpy.CreateFileGDB_management(out_folder_path=os.path.split(tempgdb)[0],
                                   out_name=os.path.split(tempgdb)[1])

#Template sr
sr_template = arcpy.Describe(coms).SpatialReference

#Output layers
bdalti_mosaic = os.path.join(pregdb, "bdalti_25m_mosaic")
bdalti_slope = os.path.join(pregdb, "bdalti_25m_slope")
bdalti_tcurvature = os.path.join(pregdb, "bdalti_25m_curvature_tangential")
bdalti_pcurvature = os.path.join(pregdb, "bdalti_25m_curvature_profile")

bdcharm_fr = os.path.join(pregdb, "lithology_bdcharm_fr")
bdcharm_bvinters = os.path.join(pregdb, "lithology_bdcharm_fr_bvinters")
bdcharm_bvinters_tab = os.path.join(resdir, "lithology_bdcharm_fr_bvinters.csv")

bdforet_fr = os.path.join(pregdb, "bdforet_fr")
bdforet_bvinters = os.path.join(pregdb, "bdforet_fr_bvinters")
bdforet_bvinters_tab = os.path.join(resdir, "bdforet_fr_bvinters.csv")

bdhaies_fr = os.path.join(pregdb, "bdhaies_fr")
bdhaies_bvinters = os.path.join(pregdb, "bdhaies_bvinters")
bdhaies_bvinters_tab = os.path.join(pregdb, "bdhaies_bvinters.csv")

irrig_coms_formatted = os.path.join(resdir, "agreste_irrig_deps_formatted.csv")
irrig_coms_poly = os.path.join(pregdb, "communes_irrig")
irrig_coms_bvinters = os.path.join(pregdb, "communes_irrig_bvinters")
irrig_coms_bvinters_tab = os.path.join(resdir, "communes_irrig_bvinters.csv")

barriers_pts_proj = os.path.join(pregdb, 'barriers_amber_proj')
barriers_pts_bvinters = os.path.join(pregdb, 'barriers_amber_bvinters')
barriers_pts_bvinters_tab = os.path.join(resdir, 'barriers_amber_bvinters.csv')

withdrawal_pts_proj = os.path.join(pregdb, "withdrawals_bnpe_proj")
withdrawal_pts_bvinters = os.path.join(pregdb, "withdrawals_bnpe_proj_bvinters")
withdrawal_pts_bvinters_tab = os.path.join(resdir, "withdrawals_bnpe_proj_bvinters.csv")

onde_stations_pts = os.path.join(pregdb, "onde_stations")
onde_stations_pts_bvinters = os.path.join(pregdb, "onde_stations_bvinters")
onde_stations_pts_bvinters_tab = os.path.join(resdir, "onde_stations_bvinters.csv")

snelder_ires_bvinters = os.path.join(tempgdb, "snelder_ires_bvinters")
snelder_ires_bvinters_tab = os.path.join(resdir, "snelder_ires_bvinters.csv")

fish_stations_pts = os.path.join(pregdb, 'fish_stations_aspe')
hydrobio_stations_pts = os.path.join(pregdb, "hydrobio_stations_naiade")

buildings_filtered_fr = os.path.join(pregdb, "buildings_filtered_fr")

#--------------------------------- DEM - BDALTI ------------------------------------------------------------------------
#Based on metadata from https://geoservices.ign.fr/documentation/donnees/alti/bdalti
#https://geoservices.ign.fr/sites/default/files/2021-07/DC_BDALTI_2-0.pdf
bdalti_sr = arcpy.SpatialReference(5698)

if not arcpy.Exists(bdalti_mosaic):
    print("Processing {}...".format(bdalti_mosaic))
    bdalti_filelist = getfilelist(bdalti_dir, "^BDALTIV2_25M_FXX_.*[.]asc$")
    arcpy.MosaicToNewRaster_management(input_rasters=bdalti_filelist,
                                       output_location=os.path.split(bdalti_mosaic)[0],
                                       raster_dataset_name_with_extension=os.path.split(bdalti_mosaic)[1],
                                       coordinate_system_for_the_raster=bdalti_sr,
                                       pixel_type="16_BIT_SIGNED",
                                       number_of_bands=1
                                       )

#Compute slope
if not arcpy.Exists(bdalti_slope):
    print("Processing {}...".format(bdalti_slope))
    Slope(in_raster=bdalti_mosaic, output_measurement="DEGREE").save(bdalti_slope)

#Compute curvature
"""Tangential curvature: Positive values indicate areas of diverging surface flow. Negative tangential curvatures 
indicates areas of converging surface flow. A positive tangential (normal contour) curvature indicates that the 
surface is convex at that cell perpendicular to the direction of the slope. A negative curvature indicates that the 
surface is concave at that cell in the direction perpendicular to the slope. A value of 0 indicates that the surface 
is flat."""
if not arcpy.Exists(bdalti_tcurvature):
    print("Processing {}...".format(bdalti_tcurvature))
    SurfaceParameters(in_raster=bdalti_mosaic,
                      parameter_type="TANGENTIAL_CURVATURE",
                      local_surface_type = "QUADRATIC",
                      z_unit='METER'
                      ).save(bdalti_tcurvature)
"""Profile curvature: Positive values indicate areas of acceleration of surface flow and erosion. Negative profile
 curvature indicates areas of slowing surface flow and deposition. A positive profile (normal slope line) curvature 
 indicates that the surface is convex at that cell in the direction of the slope. A negative curvature indicates that 
 the surface is concave at that cell in that same direction. A value of 0 indicates that the surface is flat."""
if not arcpy.Exists(bdalti_pcurvature):
    print("Processing {}...".format(bdalti_pcurvature))
    SurfaceParameters(in_raster=bdalti_mosaic,
                      parameter_type="PROFILE_CURVATURE",
                      local_surface_type = "QUADRATIC",
                      z_unit='METER'
                      ).save(bdalti_pcurvature)

#Compute topographic wetness index in the future - https://mapscaping.com/topographic-wetness-index-in-arcgis-pro/

#--------------------------------- Lithology - BDCHARM -----------------------------------------------------------------
if not arcpy.Exists(bdcharm_fr):
    print("Processing {}...".format(bdcharm_fr))
    bdcharm_filelist = getfilelist(bdcharm_dir, "^GEO050K_HARM_[0-9]{3}_S_FGEOL_2154[.]shp$")
    arcpy.management.Merge(inputs=bdcharm_filelist, output=bdcharm_fr)

if not arcpy.Exists(bdcharm_bvinters_tab):
    arcpy.analysis.Intersect([bdcharm_fr, cats_hybasdeps],
                             out_feature_class=bdcharm_bvinters,
                             join_attributes="ALL")
    arcpy.CopyRows_management(in_rows=bdcharm_bvinters,
                              out_table=bdcharm_bvinters_tab)
    arcpy.Delete_management(bdcharm_bvinters)

#--------------------------------- Forest type - BD foret  -------------------------------------------------------------
if not arcpy.Exists(bdforet_fr):
    print("Processing {}...".format(bdforet_fr))
    bdforet_filelist = getfilelist(bdforet_dir, "FORMATION_VEGETALE.shp$")
    arcpy.management.Merge(inputs=bdforet_filelist, output=bdforet_fr)

if not arcpy.Exists(bdforet_bvinters_tab):
    arcpy.analysis.Intersect([bdforet_fr, cats_hybasdeps],
                             out_feature_class=bdforet_bvinters,
                             join_attributes="ALL")
    arcpy.CopyRows_management(in_rows=bdforet_bvinters,
                              out_table=bdforet_bvinters_tab)
    arcpy.Delete_management(bdforet_bvinters)


#--------------------------------- Hedges - BDHaies  -------------------------------------------------------------------
if not arcpy.Exists(bdhaies_fr):
    print("Processing {}...".format(bdhaies_fr))
    bdtopo2023_gpkglist = getfilelist(bdtopo2023_dir, "BDT_3-3_GPKG_LAMB93_R[0-9]{2}-ED2023-09-15.gpkg$")
    bdhaies_filelist = [os.path.join(gpkg, "main.haie") for gpkg in bdtopo2023_gpkglist]
    arcpy.management.Merge(inputs=bdhaies_filelist, output=bdhaies_fr)

if not arcpy.Exists(bdhaies_bvinters_tab):
    arcpy.analysis.Intersect([bdhaies_fr, cats_hybasdeps],
                             out_feature_class=bdhaies_bvinters,
                             join_attributes="ALL")
    arcpy.CopyRows_management(in_rows=bdhaies_bvinters,
                              out_table=bdhaies_bvinters_tab)
    arcpy.Delete_management(bdhaies_bvinters)

#--------------------------------- Get % irrigated area by commune -----------------------------------------------------
if not arcpy.Exists(irrig_coms_poly):
    #Read tables
    irrig_deps_pd = pd.read_csv(filepath_or_buffer=irrigation_departements, sep=";", skiprows=2)
    irrig_coms_pd = pd.read_csv(filepath_or_buffer=irrigation_communes, sep=";", skiprows=2)

    # Format column names
    irrig_deps_pd.columns = ["{}_dep".format(
        re.sub("[*=,()]", "",
               re.sub("(\\s[-]\\s)|\\s", "_", col)))
        for col in irrig_deps_pd.columns]
    irrig_coms_pd.columns = ["{}_com".format(
        re.sub("[*=,()]", "",
               re.sub("(\\s[-]\\s)|\\s", "_", col)))
        for col in irrig_coms_pd.columns]

    #Merge tables
    irrig_coms_pd['Code_dep'] = irrig_coms_pd.Code_com.str[0:2]
    irrig_comsdepsmerge_pd = irrig_coms_pd.merge(irrig_deps_pd, on="Code_dep", how='left')

    #Remove French guyana
    irrig_comsdepsmerge_pd = irrig_comsdepsmerge_pd.loc[irrig_comsdepsmerge_pd['Code_dep'] != '97']

    irrig_comsdepsmerge_pd['irrig_area_com'] =  (irrig_comsdepsmerge_pd.SAU_en_2020_com
                                                 * irrig_comsdepsmerge_pd.Part_des_superficies_irriguées_dans_la_SAU_2020_com/100.0)
    irrig_comsdepsmerge_pd['irrig_area_com_depsum'] = irrig_comsdepsmerge_pd.groupby('Code_dep')['irrig_area_com'].transform('sum')
    irrig_comsdepsmerge_pd['irrig_area_dep'] = (irrig_comsdepsmerge_pd.SAU_en_2020_dep
                                                 * irrig_comsdepsmerge_pd.Part_des_superficies_irriguées_dans_la_SAU_2020_dep/100.0)
    #In Doubs, 0% of SAU is irrigated so won't interpolate
    irrig_comsdepsmerge_pd.loc[(irrig_comsdepsmerge_pd.irrig_area_com_depsum > irrig_comsdepsmerge_pd.irrig_area_dep)]

    irrig_comsdepsmerge_pd['irrig_area_comsumdep_diff'] = (irrig_comsdepsmerge_pd['irrig_area_dep']
                                                           - irrig_comsdepsmerge_pd['irrig_area_com_depsum'])
    irrig_comsdepsmerge_pd.loc[irrig_comsdepsmerge_pd.irrig_area_comsumdep_diff <0, 'irrig_area_comsumdep_diff'] = 0 #Make sure 0 for Doubs

    #Interpolate area of irrigated agriculture with communes with NAs due to anonymity by allocating remaining irrigated
    #area in the department proportionaly to the SAU (utilized agricultural surface area)
    #If 1000 ha of irrigated area are unassigned in a department, and a commune has 10% of the SAU with NA in the department
    #then it is assigned 10ha
    irrig_comsdepsmerge_pd['SAU_with_NAirrig_com'] = (
            irrig_comsdepsmerge_pd['Part_des_superficies_irriguées_dans_la_SAU_2020_com'].isna()
            *irrig_comsdepsmerge_pd.SAU_en_2020_com)
    irrig_comsdepsmerge_pd['SAU_with_NAirrig_depsum'] = (irrig_comsdepsmerge_pd.groupby('Code_dep')['SAU_with_NAirrig_com'].
                                                         transform('sum'))

    (irrig_comsdepsmerge_pd['irrig_area_comsumdep_diff']>irrig_comsdepsmerge_pd['SAU_with_NAirrig_depsum']).sum()

    irrig_comsdepsmerge_pd['irrig_area_com_interp'] = (irrig_comsdepsmerge_pd['irrig_area_comsumdep_diff']
                                                       *irrig_comsdepsmerge_pd['SAU_with_NAirrig_com']
                                                       /irrig_comsdepsmerge_pd['SAU_with_NAirrig_depsum'])

    irrig_comsdepsmerge_pd['irrig_area_com_all'] = (irrig_comsdepsmerge_pd.irrig_area_com.fillna(0) +
                                                    irrig_comsdepsmerge_pd.irrig_area_com_interp.fillna(0))

        #Check that numbers add up
        # check = pd.concat(
        #     [irrig_comsdepsmerge_pd.groupby('Code_dep')['irrig_area_com_all'].sum(),
        #      irrig_comsdepsmerge_pd.groupby('Code_dep')['irrig_area_dep'].unique()],
        #     axis=1
        # )
        # check['check'] = (check.irrig_area_com_all != check.irrig_area_dep)
    #Export table
    if not arcpy.Exists(irrig_coms_formatted):
        irrig_comsdepsmerge_pd.to_csv(irrig_coms_formatted)

    #Join it to the commune shapefile
    arcpy.MakeFeatureLayer_management(coms, 'coms_lyr', where_clause="NOT (INSEE_DEP IN ('2A', '2B'))")
    arcpy.management.AddJoin(in_layer_or_view='coms_lyr', in_field='INSEE_COM',
                             join_table=irrig_coms_formatted, join_field='Code_com')
    arcpy.CopyFeatures_management('coms_lyr', irrig_coms_poly)

if not arcpy.Exists(irrig_coms_bvinters_tab):
    #Intersect with bv and export
    arcpy.analysis.Intersect([irrig_coms_poly, cats_hybasdeps],
                             out_feature_class=irrig_coms_bvinters,
                             join_attributes="ALL")
    arcpy.CopyRows_management(in_rows=irrig_coms_bvinters,
                              out_table=irrig_coms_bvinters_tab)
    arcpy.Delete_management(irrig_coms_bvinters)

#--------------------------------- Create barrier points --------------------------------------------------------------
if not arcpy.Exists(barriers_pts_proj):
    barriers_pts = os.path.join(pregdb, 'barriers_amber_wgs84')
    arcpy.management.XYTableToPoint(in_table=barriers, out_feature_class=barriers_pts,
                                    x_field='Longitude_WGS84', y_field='Latitude_WGS84',
                                    coordinate_system=arcpy.SpatialReference(4326))
    arcpy.management.Project(barriers_pts, barriers_pts_proj, out_coor_system=sr_template)
    arcpy.Delete_management(barriers_pts)

if not arcpy.Exists(barriers_pts_bvinters_tab):
    #Intersect with bv and export
    arcpy.analysis.Intersect([barriers_pts_proj, cats_hybasdeps],
                             out_feature_class=barriers_pts_bvinters,
                             join_attributes="ALL")
    arcpy.CopyRows_management(in_rows=barriers_pts_bvinters,
                              out_table=barriers_pts_bvinters_tab)
    arcpy.Delete_management(barriers_pts_bvinters)

#--------------------------------- Create water withdrawal points  -----------------------------------------------------
if not arcpy.Exists(withdrawal_pts_proj):
    withdrawal_pts = os.path.join(pregdb, 'withdrawal_bnpe_wgs84')
    arcpy.management.XYTableToPoint(in_table=withdrawal_stations_bnpe, out_feature_class=withdrawal_pts,
                                    x_field='longitude', y_field='latitude',
                                    coordinate_system=arcpy.SpatialReference(4326))
    arcpy.management.Project(withdrawal_pts, withdrawal_pts_proj, out_coor_system=sr_template)
    arcpy.Delete_management(withdrawal_pts)

if not arcpy.Exists(withdrawal_pts_bvinters_tab):
    #Intersect with bv and export
    arcpy.analysis.Intersect([withdrawal_pts_proj, cats_hybasdeps],
                             out_feature_class=withdrawal_pts_bvinters,
                             join_attributes="ALL")
    arcpy.CopyRows_management(in_rows=withdrawal_pts_bvinters,
                              out_table=withdrawal_pts_bvinters_tab)
    arcpy.Delete_management(withdrawal_pts_bvinters)

#--------------------------------- Create ONDE points  ---------------------------------------------------------
if not arcpy.Exists(onde_stations_pts):
    onde_catpd = pd.concat(
        [pd.read_csv(onde_yr) for onde_yr in onde_filelist]
    )
    onde_catpd.columns = [re.sub('[<>]', '', col) for col in onde_catpd.columns]

    onde_cat = os.path.join(resdir, "onde_obs_all.csv")
    onde_catpd.to_csv(onde_cat)

    arcpy.management.XYTableToPoint(in_table=onde_cat, out_feature_class=onde_stations_pts,
                                    x_field='CoordXSiteHydro', y_field='CoordYSiteHydro',
                                    coordinate_system=sr_template)
    arcpy.management.DeleteIdentical(in_dataset=onde_stations_pts, fields='Shape', xy_tolerance='1 meter')

if not arcpy.Exists(onde_stations_pts_bvinters_tab):
    #Intersect with bv and export
    arcpy.analysis.Intersect([onde_stations_pts, cats_hybasdeps],
                             out_feature_class=onde_stations_pts_bvinters,
                             join_attributes="ALL")
    arcpy.CopyRows_management(in_rows=onde_stations_pts_bvinters,
                              out_table=onde_stations_pts_bvinters_tab)
    arcpy.Delete_management(onde_stations_pts_bvinters)

#--------------------------------- Snelder intersection  --------------------------------------------------------------
if not arcpy.Exists(snelder_ires_bvinters_tab):
    #Intersect with bv and export
    arcpy.analysis.Intersect([snelder_ires, cats_hybasdeps],
                             out_feature_class=snelder_ires_bvinters,
                             join_attributes="ALL")
    arcpy.CopyRows_management(in_rows=snelder_ires_bvinters,
                              out_table=snelder_ires_bvinters_tab)
    arcpy.Delete_management(snelder_ires_bvinters)

#--------------------------------- Format INSEE population data  -------------------------------------------------------
#Filter and merge buildings
if not arcpy.Exists(buildings_filtered_fr):
    buildings_filelist = getfilelist(bdtopo2019_dir, "BATIMENT.shp$")
    #buildings_shp = buildings_filelist[0]

    temp_buildings_list= []
    for buildings_shp in buildings_filelist:
        temp_lyr=os.path.join(tempgdb,
                              "{}_batiments".format(
                                  re.sub('[-]', '_',
                                         Path(buildings_shp).parts[-3])
                              ))
        if not arcpy.Exists(temp_lyr):
            print("Processing {}".format(os.path.split(temp_lyr)[1]))

            arcpy.MakeFeatureLayer_management(in_features=buildings_shp,
                                   out_layer='buildings_lyr',
                                   where_clause=("((USAGE1 IN ('Résidentiel', 'Indifférencié')) "
                                                 "OR (USAGE2 IN ('Résidentiel', 'Indifférencié'))) "
                                                 "AND (NATURE IN ('Indifférenciée', 'Château')) "
                                                 "AND (LEGER = 'Non') "
                                                 "AND (ETAT='En service') "
                                                 "AND ((HAUTEUR > 2) OR (HAUTEUR = 0))"))
            arcpy.CopyFeatures_management(in_features='buildings_lyr',
                                          out_feature_class=temp_lyr
                                          )
            temp_buildings_list.append(temp_lyr)
            arcpy.AddGeometryAttributes_management(temp_lyr, Geometry_Properties='AREA', Area_Unit='SQUARE_METERS')
            with arcpy.da.UpdateCursor(temp_lyr, ['POLY_AREA']) as cursor:
                for row in cursor:
                    if row[0] < 20:
                        cursor.deleteRow()
            arcpy.Delete_management('buildings_lyr')

    print("Merging all building layers...")
    arcpy.Merge_management(inputs=temp_buildings_list, output=buildings_filtered_fr)
    print("Deleting temporary building layers...")
    for temp_lyr in temp_buildings_list:
        arcpy.Delete_management(temp_lyr)

#Intersect buildings with variable mesh size

#For residential buildings without # of households, assign based on height

#For indifférencié buildings, assign 1 household?

#Apportion population in variable mesh census based on number of households in each building

#Check for outliers

#Spatial join buildings to 200-m mesh

#Compute total population in each 200-m mesh

#COnvert 200-m mesh to raster


#--------------------------------- Create fish station points  ---------------------------------------------------------
if not arcpy.Exists(fish_stations_pts):
    arcpy.management.XYTableToPoint(in_table=fish_stations, out_feature_class=fish_stations_pts,
                                    x_field='sta_coordonnees_x', y_field='sta_coordonnees_y',
                                    coordinate_system=sr_template)

#--------------------------------- Create hydrobiology station points  ---------------------------------------------------------
if not arcpy.Exists(hydrobio_stations_pts):
    hydrobio_stations_copy = os.path.join(resdir, "hydrobio_stations_nariade_copy.csv")
    pd.read_csv(hydrobio_stations, sep=";", encoding='cp1252').to_csv(hydrobio_stations_copy, encoding='utf-8', sep=";")
    arcpy.management.XYTableToPoint(in_table=hydrobio_stations_copy, out_feature_class=hydrobio_stations_pts,
                                    x_field='CoordXStationMesureEauxSurface', y_field='CoordYStationMesureEauxSurface',
                                    coordinate_system=sr_template)
    arcpy.Delete_management(hydrobio_stations_copy)