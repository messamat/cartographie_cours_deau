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

#Land cover
lc_dir = os.path.join(anci_dir, 'oso')

#Global aridity index
gai_dir = os.path.join(anci_dir, 'gaiv3', 'Global-AI_v3_monthly')

#Dominant soil type 1:1,000,000
bdgsf = os.path.join(anci_dir, "bdgsf", "30169_L93", "30169_L93.shp")

#BD topo 2023
bdtopo2023_dir = os.path.join(anci_dir, "bdtopo233") #233 is the version/package number (vs e.g. 151 for 2015)
"plan_d_eau"
"reservoir"

#Withdrawal points - bnpe
withdrawal_stations_bnpe=os.path.join(anci_dir, "bnpe", "bnpe_ouvrages.csv")

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

#Snelder intermittence
snelder_ires = os.path.join(anci_dir, "snelder", "rhtvs2_all_phi_qclass.shp")

#-------------------------------- OUTPUTS ------------------------------------------------------------------------------
#Scratch gdb
tempgdb = os.path.join(resdir, "scratch.gdb")
if not arcpy.Exists(tempgdb):
    arcpy.CreateFileGDB_management(out_folder_path=os.path.split(tempgdb)[0],
                                   out_name=os.path.split(tempgdb)[1])

#gdb to hold land cover statistics
lcav_dir = os.path.join(resdir, 'oso_lc_stats')
if not arcpy.Exists(lcav_dir):
    os.mkdir(lcav_dir)

#Template sr
sr_template = arcpy.SpatialReference(2154)

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
bdhaies_bvinters_tab = os.path.join(resdir, "bdhaies_bvinters.csv")

oso_veg = os.path.join(lcav_dir, 'oso_veg')
oso_imp = os.path.join(lcav_dir, 'oso_imp')
oso_agr = os.path.join(lcav_dir, 'oso_agr')
oso_scr = os.path.join(lcav_dir, 'oso_scr')

gai_yr = os.path.join(pregdb, "ai_v3_yrav")
gai_summer = os.path.join(pregdb, "ai_v3_summerav")

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
buildings_popvariable =os.path.join(pregdb, "buildings_pop_nivNaturel_inters")
pop_count_variable_mesh = os.path.join(anci_dir, "insee_pop", "carreaux_nivNaturel_met.shp")
pop_count_variable_mesh_proj = os.path.join(anci_dir, tempgdb, "carreaux_nivNaturel_met_proj")
pop_count_200m_fishnet = os.path.join(anci_dir, tempgdb, "carreaux_200m")
pop_count_200m_fishnet_lambert = os.path.join(anci_dir, tempgdb, "carreaux_200m_lambert")
pop_count_200m_fishnet_buildingsjoin = os.path.join(anci_dir, tempgdb, "carreaux_200m_buildingsjoin")
pop_ras_200m = os.path.join(pregdb, 'insee_pop_interp200m')

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

#--------------------------------- Land cover --------------------------------------------------------------------------
#Land cover - theia oso
lc_filedict = {yr: getfilelist(os.path.join(lc_dir, "oso_{}".format(yr)),"Classif_Seed_.*[.]tif$")[0]
               for yr in [2019, 2020, 2021]}
#Create dictionary of classes https://www.theia-land.fr/en/product/land-cover-map/
lc_class_dict = {1:'urba1', 2:'urba2', 3:'indus', 4:'roads', 5:'wioil', 6:'straw', 7:'spoil', 8:"soy", 9:"sunfl",
              10:"corn", 11:"rice", 12:"roots", 13:"pastu", 14:"orchd", 15:"vinyd", 16:"forbr", 17:"forco",
              18:"grass", 19:"heath", 20:"rocks", 21:"beach", 22:'glasn', 23:'water'}

#Tile raster (SplitRaster function does not work with arcpy)
with arcpy.EnvManager(snapRaster=lc_filedict[2019]):
    for cl in lc_class_dict:
        out_cl = os.path.join(lcav_dir, 'oso_cl{}.tif'.format(str(cl).zfill(2)))
        #start = time.time()
        if not arcpy.Exists(out_cl):
            print("Processing {}...".format(out_cl))
            CellStatistics(in_rasters_or_constants=[(Raster(lc_filedict[2019])==cl),
                                                    (Raster(lc_filedict[2020])==cl),
                                                    (Raster(lc_filedict[2021])==cl)],
                           statistics_type='SUM').save(out_cl)
        #print(time.time() - start)

if not arcpy.Exists(oso_veg):
    CellStatistics(
        in_rasters_or_constants=[os.path.join(lcav_dir, 'oso_cl{}.tif'.format(str(cl).zfill(2)))
                                 for cl in [16,17,18,19]],
        statistics_type='SUM').save(oso_veg)

if not arcpy.Exists(oso_imp):
    CellStatistics(
        in_rasters_or_constants=[os.path.join(lcav_dir, 'oso_cl{}.tif'.format(str(cl).zfill(2))) for cl in [1,2,3,4]],
        statistics_type='SUM').save(oso_imp)

if not arcpy.Exists(oso_agr):
    CellStatistics(
        in_rasters_or_constants=[os.path.join(lcav_dir, 'oso_cl{}.tif'.format(str(cl).zfill(2))) for cl in range(5,16)],
        statistics_type='SUM').save(oso_agr)

if not arcpy.Exists(oso_scr):
    CellStatistics(
        in_rasters_or_constants=[os.path.join(lcav_dir, 'oso_cl{}.tif'.format(str(cl).zfill(2))) for cl in range(8,13)],
        statistics_type='SUM').save(oso_scr)

#--------------------------------- Global aridity index  ---------------------------------------------------------------
gai_filedict = {mo: os.path.join(gai_dir, "ai_v3_{}.tif".format(str(mo).zfill(2))) for mo in range(1,13)}

with arcpy.EnvManager(extent="-5, 42, 9, 52", snapRaster=gai_filedict[1]):
    if not arcpy.Exists(gai_yr):
        CellStatistics(in_rasters_or_constants=list(gai_filedict.values()),
                       statistics_type='MEAN'
                       ).save(gai_yr)
    if not arcpy.Exists(gai_summer):
        CellStatistics(in_rasters_or_constants=[gai_filedict[mo] for mo in [7,8,9]],
                       statistics_type='MEAN'
                       ).save(gai_summer)

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
    arcpy.AlterField_management(buildings_filtered_fr, 'POLY_AREA',
                                'POLY_AREA_WHOLE', 'POLY_AREA_WHOLE')
    print("Deleting temporary building layers...")
    for temp_lyr in temp_buildings_list:
        arcpy.Delete_management(temp_lyr)

#Intersect buildings with variable mesh size
if not arcpy.Exists(buildings_popvariable):
    arcpy.Intersect_analysis([buildings_filtered_fr, pop_count_variable_mesh],
                             out_feature_class=buildings_popvariable,
                             join_attributes='ALL')

#Pre-process attributes
arcpy.AddGeometryAttributes_management(buildings_popvariable, Geometry_Properties='AREA', Area_Unit='SQUARE_METERS')
arcpy.AlterField_management(buildings_popvariable, 'POLY_AREA',
                            'POLY_AREA_MESHED', 'POLY_AREA_MESHED')

if not 'log_total' in [f.name for f in arcpy.ListFields(buildings_popvariable)]: #Compute total number of households in each quadrat
    arcpy.AddField_management(buildings_popvariable, 'log_total', 'FLOAT')
if not 'VOLUME' in [f.name for f in arcpy.ListFields(buildings_popvariable)]: #Volume of each building
    arcpy.AddField_management(buildings_popvariable, 'VOLUME', 'FLOAT')
if not 'avind_per_log' in [f.name for f in arcpy.ListFields(buildings_popvariable)]: #Average number of individuals per housing unit in quadrat
    arcpy.AddField_management(buildings_popvariable, 'avind_per_log', 'FLOAT')

buildings_lgts_meshdict = defaultdict(int)
buildings_nolgt_totalvol_meshdict = defaultdict(float)

with arcpy.da.UpdateCursor(buildings_popvariable, ['POLY_AREA_MESHED', 'POLY_AREA_WHOLE',
                                                   'log_av45', 'log_45_70', 'log_70_90',
                                                   'log_ap90', 'log_inc', 'log_total',
                                                   'HAUTEUR', 'VOLUME',
                                                   'avind_per_log', 'ind',
                                                   'idcar_nat', 'NB_LOGTS'
                                                   ]) as cursor:
    for row in cursor:
        if (row[0]/row[1])<0.5: #Remove parts under half of a building and compute number of housing units based on census
            cursor.deleteRow()
        else:
            log_sum = row[2] + row[3] + row[4] + row[5] + row[6]

            # Compute total number of households in each mesh
            row[7] = log_sum

            # Compute volume of each building
            if row[8] > 0:
                row[9] = row[8]*row[1] #Volume = height*surface area
            else:
                row[9] = 4*row[1] #If no height data is available, assign standard "hauteur au faitage" of a single-floor house

            #Compute average number of individuals per housing unit in quadrat according to census
            if log_sum > 0:
                row[10] = row[11]/log_sum
            cursor.updateRow(row)
            del log_sum

            #Count the total volume of buildings without a registered number of housing units in each census quadrat
            buildings_nolgt_totalvol_meshdict[row[12]] += row[9]
            #Count the registered number of housing units based on building attributes in each census quadrat
            buildings_lgts_meshdict[row[12]] += row[13]

#For residential buildings without # of households, assign number of housing units based on volume, the compute pop/building
if not 'NB_LOGTS_EST' in [f.name for f in arcpy.ListFields(buildings_popvariable)]:
    arcpy.AddField_management(buildings_popvariable, 'NB_LOGTS_EST', 'FLOAT')
if not 'ind_est' in [f.name for f in arcpy.ListFields(buildings_popvariable)]:
    arcpy.AddField_management(buildings_popvariable, 'ind_est', 'FLOAT')

with arcpy.da.UpdateCursor(buildings_popvariable, ['idcar_nat', 'log_total', 'NB_LOGTS_EST', 'NB_LOGTS',
                                                    'VOLUME', 'ind_est', 'avind_per_log']) as cursor:
    for row in cursor:
        #Adjust number of household units per building
        if row[1] > 0: #If there are housing units in the quadrat according to the census
            #If there are some buildings with a registered number of households in the quadrat
            if buildings_lgts_meshdict[row[0]] > 0:
                lgts_ratio_meshtobuildings = row[1]/buildings_lgts_meshdict[row[0]] #Ratio of census- vs. buildings-based number of housing units

                #(If total number of building-based housing units exceeds that in the census)
                #OR
                #(If total number of building-based housing units is inferior to that in the census,
                # And there are no buildings without housing units.)
                # ----> Adjust number of housing units in buildings with registered number of housing
                if ((lgts_ratio_meshtobuildings<1) or
                    ((lgts_ratio_meshtobuildings>1) and (buildings_nolgt_totalvol_meshdict[row[0]] == 0))
                ):
                    row[2] = row[3]*lgts_ratio_meshtobuildings

                #If total number of building-based household units is inferior to that in the census
                # AND there are buildings without household units, assign units to those based on their volume
                elif ((lgts_ratio_meshtobuildings>1) and (buildings_nolgt_totalvol_meshdict[row[0]] > 0)):
                    if row[3] == 0:
                        row[2] = (row[1]-buildings_lgts_meshdict[row[0]]
                                  )*(row[4]/buildings_nolgt_totalvol_meshdict[row[0]]) #Diff in number of housing units*proportion of unassigned building volume in quadrat in this building
                    else: #Keep the registered number of housing units for the others
                        row[2] = row[3]
                #Otherwise, keep the same number of household units in building
                else:
                    row[2] = row[3]
            #If there are no buildings with a registered number of households in the quadrat
            else:
                #Assign units based on volume
                row[2] = row[1] * (row[4] / buildings_nolgt_totalvol_meshdict[row[0]])  # Diff in number of housing units*proportion of unassigned building volume in quadrat in this building
            cursor.updateRow(row)

            #Adjust number of household units per building if there are people in the quadrat
            if row[6] > 0:
                row[5] = row[2]*row[6] #Individuals in building= estimated number of household units in building*average number of individuals per housing unit in quadrat
                cursor.updateRow(row)

#Create matching continuous 200 m fishnet---------------------
arcpy.Project_management(pop_count_variable_mesh, out_dataset=pop_count_variable_mesh_proj,
                         out_coor_system=arcpy.SpatialReference(3035)) # The mesh was produced based on data projected in LAEA (EPSG 3035), the European cs
variable_mesh_ext = arcpy.Describe(pop_count_variable_mesh_proj).Extent
with arcpy.EnvManager(outputCoordinateSystem=arcpy.SpatialReference(3035)):
    arcpy.CreateFishnet_management(
        out_feature_class=pop_count_200m_fishnet,
        origin_coord="{0} {1}".format(variable_mesh_ext.XMin, variable_mesh_ext.YMin),
        y_axis_coord="{0} {1}".format(variable_mesh_ext.XMin, variable_mesh_ext.YMax),
        cell_width=200,
        cell_height=200,
        corner_coord="{0} {1}".format(variable_mesh_ext.XMax, variable_mesh_ext.YMax),
        geometry_type="POLYGON")

#Reproject to Lambert ---------------------
arcpy.Project_management(pop_count_200m_fishnet, out_dataset=pop_count_200m_fishnet_lambert,
                         out_coor_system=sr_template)

#Compute total population in each 200-m quadrat, summing pop and housing units for all buildings ---------------------
# create a list of fields to sum
f_tosum = ['NB_LOGTS_EST', 'ind_est']

# create the field mapping object
fms_pop = arcpy.FieldMappings()

# populate the field mapping object with the fields from both feature classes
fms_pop.addTable(pop_count_200m_fishnet_lambert)
fms_pop.addTable(buildings_popvariable)

# loop through the field names to sum
for fieldName in f_tosum:
    # get the field map index of this field and get the field map
    fieldIndex = fms_pop.findFieldMapIndex(fieldName)
    fieldMap = fms_pop.getFieldMap(fieldIndex)
    # update the field map with the new merge rule
    fieldMap.mergeRule = 'Sum'
    # replace with the updated field map
    fms_pop.replaceFieldMap(fieldIndex, fieldMap)

arcpy.SpatialJoin_analysis(target_features=pop_count_200m_fishnet_lambert,
                           join_features=buildings_popvariable,
                           out_feature_class=pop_count_200m_fishnet_buildingsjoin,
                           join_operation='JOIN_ONE_TO_MANY',
                           join_type='KEEP_ALL',
                           field_mapping=fms_pop,
                           match_option='LARGEST_OVERLAP'
                           )

with arcpy.da.UpdateCursor(pop_count_200m_fishnet_buildingsjoin, [f_tosum]) as cursor:
    for row in cursor:
        if row[0] == None:
            row[0] = 0
        if row[1] == None:
            row[1] = 0
        cursor.updateRow(row)

#Convert 200-m mesh to raster------------------------------
arcpy.PolygonToRaster_conversion(in_features=pop_count_200m_fishnet_buildingsjoin,
                                 value_field='ind_est',
                                 out_rasterdataset=pop_ras_200m,
                                 cell_assignment='CELL_CENTER',
                                 cellsize=200)



###################### EXTRA stuff ##################
# for yr in lc_filedict:
#     lcyr_gdb = os.path.join(lcav_dir, 'lc{}_tiles.gdb'.format(yr))
#     if not arcpy.Exists(lcyr_gdb):
#         arcpy.CreateFileGDB_management(out_folder_path=os.path.split(lcyr_gdb)[0],
#                                        out_name=os.path.split(lcyr_gdb)[1])
#
#     lcext = arcpy.Describe(lc_filedict[yr]).Extent
#     lc_bbox = [lcext.XMin, lcext.YMin, lcext.XMax, lcext.YMax]
#     lc_tiles_bblist = divbb(bbox=lc_bbox,
#                      res=arcpy.Describe(lc_filedict[yr]).meanCellWidth,
#                      divratio=6)
#     if len(getfilelist(lcyr_gdb, gdbf=True)) < len(lc_tiles_bblist):
#         x=1
#         for tile_bb in lc_tiles_bblist:
#             out_tile = os.path.join(lcyr_gdb, 'lc{0}_{1}'.format(yr, x))
#             if x > 34:
#                 if not arcpy.Exists(out_tile):
#                     print('Processing {}...'.format(out_tile))
#                     #arcpy.env.extent = ' '.join(map(str, tile_bb))
#                     arcpy.Clip_management(in_raster=lc_filedict[yr],
#                                           rectangle=' '.join(map(str, tile_bb)),
#                                           out_raster=out_tile,
#                                           maintain_clipping_extent='NO_MAINTAIN_EXTENT')
#                     arcpy.ClearEnvironment('extent')
#                 else:
#                     print('{} already exists...'.format(out_tile))
#             x += 1