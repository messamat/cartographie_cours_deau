import arcpy.analysis

from setup_classement import  *

anci_dir = os.path.join(datdir, "données_auxiliaires")#Ancillary data directory

#Preprocessing geodatabase
pregdb = os.path.join(resdir, "preprocessing_ancillary_data.gdb")

#------------------------------------------- INPUTS --------------------------------------------------------------------
#Departments and communes - admin express
admin_dir = os.path.join(anci_dir, "admin_express", "ADMIN-EXPRESS_3-2__SHP_WGS84G_FRA_2023-05-03",
                    "1_DONNEES_LIVRAISON_2023-05-03", "ADE_3-2_SHP_WGS84G_FRA")
deps = os.path.join(admin_dir, "DEPARTEMENT.shp")
coms = os.path.join(admin_dir, "COMMUNE.shp")

#Catchments (cats) - Bassin versant topographique BD Topage
cats_raw = os.path.join(anci_dir, "topage", "BassinVersantTopographique_FXX.gpkg", "main.BassinVersantTopographique_FXX")

#Coordinate system template
fr_sr = arcpy.Describe(cats_raw).SpatialReference

#Hydrobasins level 8 and 9
hybas08 = os.path.join(anci_dir, "hydrosheds", "hybas_lake_eu_lev08_v1c.shp")
hybas09 = os.path.join(anci_dir, "hydrosheds", "hybas_lake_eu_lev09_v1c.shp")

#BCAE
bcae_dir = os.path.join(anci_dir, "bcae_20231106")
bcae_filelist = getfilelist(bcae_dir, "BCAE_lineaire.shp$")

#BDTOPO 2015
bdtopo2015_dir = os.path.join(anci_dir, "bdtopo151")
bdtopo2015_ce_filelist = getfilelist(bdtopo2015_dir, "TRONCON_COURS_EAU.SHP$")

#BD Carthage
bdcarthage = os.path.join(anci_dir, "carthage", "TRONCON_HYDROGRAPHIQUE.shp")

#RHT
rht = os.path.join(anci_dir, "rht_2020", "rht_lbt93.shp")

#------------- OUTPUTS ---------------------------------------------------------------------------------------------
cats_hybasjoin = os.path.join(pregdb, 'BV_hybas0809_join') #Bassins versants topographiques (catchments) spatially jointed to HydroBASINS levels 8 and 9
cats_hybasdeps = os.path.join(pregdb, 'BV_hybas0809_depsinters') #BV joined to HydroBASINS and intersected with Departements

bcae_fr = os.path.join(pregdb, 'bcae_fr')
bdtopo2015_fr = os.path.join(pregdb, 'bdtopo2015_fr')

#------------- PREPARE UNITS OF ANALYSIS -------------------------------------------------------------------------------
#Spatial join with HydroBASINS level 8-----------------------
if not arcpy.Exists(cats_hybasjoin):
    hybas08proj = arcpy.Project_management(hybas08, os.path.join(pregdb, 'hybas08proj'),
                                           out_coor_system=cats_raw)
    hybas09proj = arcpy.Project_management(hybas09, os.path.join(pregdb, 'hybas09proj'),
                                           out_coor_system=cats_raw)

    cats_hybasjoin08 = os.path.join(pregdb, 'BV_hybas08_join')
    arcpy.MakeFeatureLayer_management(cats_raw, out_layer='cats_raw_lyr',
                                      where_clause='NOT CdBH=12')
    arcpy.SpatialJoin_analysis(target_features='cats_raw_lyr', join_features=hybas08proj, out_feature_class=cats_hybasjoin08,
                               join_operation="JOIN_ONE_TO_ONE", join_type="KEEP_ALL", match_option="LARGEST_OVERLAP")
    arcpy.AlterField_management(cats_hybasjoin08, field="PFAF_ID",
                                new_field_name="PFAF_ID08", new_field_alias="PFAF_ID08")
    arcpy.AlterField_management(cats_hybasjoin08, field="ORDER_",
                                new_field_name="ORDER08", new_field_alias="ORDER08")

    arcpy.SpatialJoin_analysis(target_features=cats_hybasjoin08, join_features=hybas09proj, out_feature_class=cats_hybasjoin,
                               join_operation="JOIN_ONE_TO_ONE", join_type="KEEP_ALL", match_option="LARGEST_OVERLAP")
    arcpy.AlterField_management(cats_hybasjoin, field="PFAF_ID",
                                new_field_name="PFAF_ID09", new_field_alias="PFAF_ID09")
    arcpy.AlterField_management(cats_hybasjoin, field="ORDER_",
                                new_field_name="ORDER09", new_field_alias="ORDER09")

    fieldstodelete = [f.name for f in arcpy.ListFields(cats_hybasjoin)
                      if (f.name not in ([f2.name for f2 in arcpy.ListFields(cats_raw)]
                                         +['PFAF_ID08', 'PFAF_ID09','ORDER08', 'ORDER09', 'OBJECTID',
                                           'geom_Length', 'geom_Area']))]
    arcpy.DeleteField_management(cats_hybasjoin, drop_field=fieldstodelete)

    arcpy.management.Delete(os.path.join(pregdb, 'hybas08proj'))
    arcpy.management.Delete(os.path.join(pregdb, 'hybas09proj'))

#Intersect with Departement----------------------
if not arcpy.Exists(cats_hybasdeps):
    arcpy.analysis.Intersect(in_features=[cats_hybasjoin, deps], out_feature_class=cats_hybasdeps,
                             join_attributes='ALL')

#Compute area-------------------------------------
if 'POLY_AREA' not in [f.name for f in arcpy.ListFields(cats_hybasdeps)]:
    arcpy.management.AddGeometryAttributes(cats_hybasdeps, Geometry_Properties="AREA", Area_Unit='SQUARE_METERS')

#------------- PREPARE NETWORKS ----------------------------------------------------------------------------------------
#Merge BCAE
if not arcpy.Exists(bcae_fr):
    arcpy.management.Merge(inputs=bcae_filelist, output=bcae_fr)

#Merge BD Topo 2015
if not arcpy.Exists(bdtopo2015_fr):
    arcpy.management.Merge(inputs=bdtopo2015_ce_filelist, output=bdtopo2015_fr)

#------------- INTERSECT WITH UNITS OF ANALYSIS ------------------------------------------------------------------------
for net in [bdtopo2015_fr, bcae_fr, bdcarthage, rht]:
    root_name = os.path.split(os.path.splitext(net)[0])[1]
    out_net = os.path.join(pregdb, "{}_bvinters".format(root_name))
    out_tab = os.path.join(resdir, "{}.csv".format(root_name))
    if not arcpy.Exists(out_net):
        print("Processing {}...".format(out_net))
        arcpy.analysis.Intersect([net, cats_hybasdeps], out_feature_class=out_net, join_attributes="ALL")
    if not arcpy.Exists(out_tab):
        print("Exporting {}...".format(out_tab))
        arcpy.CopyRows_management(out_net, out_tab)