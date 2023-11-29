import arcpy.analysis

from setup_classement import *

overwrite=False

anci_dir = os.path.join(datdir, "données_auxiliaires")#Ancillary data directory

#Preprocessing geodatabase
pregdb = os.path.join(resdir, "preprocessing_ancillary_data.gdb")

#Scratch gdb
tempgdb = os.path.join(resdir, "scratch.gdb")
if not arcpy.Exists(tempgdb):
    arcpy.CreateFileGDB_management(out_folder_path=os.path.split(tempgdb)[0],
                                   out_name=os.path.split(tempgdb)[1])
#------------------------------------------- INPUTS --------------------------------------------------------------------
#Departments and communes - admin express
admin_dir = os.path.join(anci_dir, "admin_express", "ADMIN-EXPRESS_3-2__SHP_WGS84G_FRA_2023-05-03",
                         "1_DONNEES_LIVRAISON_2023-05-03", "ADE_3-2_SHP_WGS84G_FRA")
deps = os.path.join(admin_dir, "DEPARTEMENT.shp")
coms = os.path.join(admin_dir, "COMMUNE.shp")

#Catchments (cats) - Bassin versant topographique BD Topage
cats_raw = os.path.join(anci_dir, "bdtopage", "BassinVersantTopographique_FXX.gpkg", "main.BassinVersantTopographique_FXX")

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

#DDT network
outputs_gdb = os.path.join(resdir, 'analysis_outputs.gdb')
ce_net = os.path.join(outputs_gdb, 'carto_loi_eau_fr')

#------------- OUTPUTS ---------------------------------------------------------------------------------------------
cats_hybasjoin = os.path.join(pregdb, 'BV_hybas0809_join') #Bassins versants topographiques (catchments) spatially jointed to HydroBASINS levels 8 and 9
cats_hybasdeps = os.path.join(pregdb, 'BV_hybas0809_depsinters') #BV joined to HydroBASINS and intersected with Departements
cats_hybasdeps_tab = os.path.join(resdir, 'BV_hybas0809_depsinters.csv')

bcae_fr = os.path.join(pregdb, 'bcae_fr')
bdtopo2015_fr = os.path.join(pregdb, 'bdtopo2015_fr')
bdcarthage_nodupli = os.path.join(pregdb, 'bdcarthage_nodupli')

ddtnet_bdtopo_inters = os.path.join(pregdb, 'carto_loi_eau_bdtopo_inters')
ddtnet_carthage_inters = os.path.join(pregdb, 'carto_loi_eau_carthage_inters')
ddtnet_bdtopo_inters_tab = os.path.join(resdir, 'carto_loi_eau_bdtopo_inters_tab.csv')
ddtnet_carthage_inters_tab = os.path.join(resdir, 'carto_loi_eau_carthage_inters_tab.csv')

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
    arcpy.management.AddGeometryAttributes(cats_hybasdeps, Geometry_Properties="AREA", Area_Unit='SQUARE_KILOMETERS')

#Create stable UID
if 'UID_BV' not in [f.name for f in arcpy.ListFields(cats_hybasdeps)]:
    arcpy.AddField_management(in_table=cats_hybasdeps, field_name='UID_BV', field_type='LONG')
    with arcpy.da.UpdateCursor(cats_hybasdeps, ['UID_BV', 'OID@']) as cursor:
        for row in cursor:
            row[0] = row[1]
            cursor.updateRow(row)

#Export table
if not arcpy.Exists(cats_hybasdeps_tab):
    arcpy.CopyRows_management(in_rows=cats_hybasdeps,
                              out_table=cats_hybasdeps_tab)

#------------- PREPARE NETWORKS ----------------------------------------------------------------------------------------
#Merge BCAE
if not arcpy.Exists(bcae_fr):
    arcpy.management.Merge(inputs=bcae_filelist, output=bcae_fr)

#Merge BD Topo 2015
if not arcpy.Exists(bdtopo2015_fr):
    arcpy.management.Merge(inputs=bdtopo2015_ce_filelist, output=bdtopo2015_fr)
    arcpy.management.DeleteIdentical(in_dataset=bdtopo2015_fr, fields='Shape', xy_tolerance='1 meter')

#Remove duplicates in Carthage as well
if not arcpy.Exists(bdcarthage_nodupli):
    arcpy.CopyFeatures_management(bdcarthage, bdcarthage_nodupli)
    arcpy.management.DeleteIdentical(in_dataset=bdtopo2015_fr, fields='Shape', xy_tolerance='1 meter')

##------------- Identify identical geometries between the DDT network and BD Topo or Carthage --------------------------
# in_target_ft=ce_net
# in_join_ft=bdtopo2015_fr
# out_ft=ddtnet_bdtopo_inters
# tempgdb=tempgdb
# in_join_ft_suffix='bdtopo'
# in_join_ft_fieldstojoin='ID'
# tolerance='1 Meters'
def inters_linetoline(in_target_ft,
                      in_join_ft,
                      out_ft,
                      tempgdb,
                      in_join_ft_suffix,
                      in_join_ft_fieldstojoin,
                      in_target_ft_fieldstojoin,
                      tolerance):

    #Make sure in_join_ft_fieldstojoin is a list
    if type(in_join_ft_fieldstojoin) == str:
        in_join_ft_fieldstojoin = [in_join_ft_fieldstojoin]
    if type(in_target_ft_fieldstojoin) == str:
        in_target_ft_fieldstojoin = [in_target_ft_fieldstojoin]

    #Get field type for join fields
    in_join_ft_fdict = {f.name:f.type for f in arcpy.ListFields(in_join_ft)}
    ftojoin_dict = {}
    for ftojoin in in_join_ft_fieldstojoin:
        if ftojoin in in_join_ft_fdict:
            ftojoin_dict[ftojoin] = in_join_ft_fdict[ftojoin]

    #Get field type for target fields
    in_target_ft_fdict = {f.name:f.type for f in arcpy.ListFields(in_target_ft)}
    target_ftojoin_dict = {}
    for ftojoin in in_target_ft_fieldstojoin:
        if ftojoin in in_target_ft_fdict :
            target_ftojoin_dict[ftojoin] = in_target_ft_fdict[ftojoin]

    #Project the dataset
    in_join_ft_proj = os.path.join(tempgdb, 'join_ft_proj_{}'.format(in_join_ft_suffix))
    if not arcpy.Exists(in_join_ft_proj):
        if (arcpy.Describe(in_join_ft).SpatialReference != arcpy.Describe(in_target_ft).SpatialReference):
            arcpy.Project_management(in_dataset=in_join_ft,
                                     out_dataset=in_join_ft_proj,
                                     out_coor_system=in_target_ft)
        else:
            arcpy.CopyFeatures_management(in_features=in_join_ft, out_feature_class=in_join_ft_proj)

    #Add a specifically named length field
    lenf = "length_{}".format(in_join_ft_suffix)
    if lenf not in [f.name for f in arcpy.ListFields(in_join_ft_proj)]:
        arcpy.AddGeometryAttributes_management(in_join_ft_proj, Geometry_Properties='LENGTH',
                                               Length_Unit='METERS')
        arcpy.AlterField_management(in_join_ft_proj, 'LENGTH', new_field_name=lenf, new_field_alias=lenf)

    #Make a buffer of the desired tolerance around each line in the join feature class
    in_join_ft_buf = os.path.join(tempgdb, 'join_ft_buf_{}'.format(in_join_ft_suffix))
    if not arcpy.Exists(in_join_ft_buf):
        arcpy.Buffer_analysis(in_features=in_join_ft_proj,
                              out_feature_class=in_join_ft_buf,
                              buffer_distance_or_field=tolerance,
                              line_side='FULL',
                              line_end_type='FLAT',
                              dissolve_option='NONE'
                              )

    #Make sure to uniquely name the fields based on the join feature class to avoid duplicates with target feature class
    for f in arcpy.ListFields(in_join_ft_buf):
        if ((f.name not in [arcpy.Describe(in_join_ft_buf).OIDFieldName, 'Shape', 'Shape_Area', 'Shape_Length', lenf]) and
            (not re.findall("_{}$".format(in_join_ft_suffix), f.name))):
            new_fname = "{0}_{1}".format(f.name, in_join_ft_suffix)
            print('Editing {0} to {1}...'.format(f.name, new_fname))
            arcpy.AlterField_management(in_join_ft_buf, field=f.name,
                                        new_field_name=new_fname,
                                        new_field_alias=new_fname)
    for fname in in_join_ft_fieldstojoin:
        if fname in ftojoin_dict:
            new_fname = "{0}_{1}".format(fname, in_join_ft_suffix)
            ftojoin_dict[new_fname] = ftojoin_dict.pop(fname)
    ftojoin_dict[lenf]='DOUBLE'

    #Intersect the target feature class with the buffer of the join feature class
    arcpy.Intersect_analysis(in_features=[ce_net, in_join_ft_buf],
                             out_feature_class=out_ft,
                             join_attributes='ONLY_FID',
                             output_type='INPUT')

    # ldict = {row[0]:row[1] for row in arcpy.da.SearchCursor(in_join_ft_proj, ['OID@', 'length_carthage'])}
    # arcpy.AddField_management(out_ft, 'length_carthage', 'DOUBLE')
    # with arcpy.da.UpdateCursor(out_ft, ['FID_join_ft_buf_carthage', 'length_carthage']) as cursor:
    #     for row in cursor:
    #         if row[0] in ldict:
    #             row[1] = ldict[row[0]]
    #             cursor.updateRow(row)

    #Join the desired fields from the join feature class to the target feature class
    def join_field_fromdict(in_tab, out_tab, field_dict):
        ftojoin_valid = {}
        for fname, ftype in field_dict.items():
            if ftype=='String':
                ftojoin_valid[fname]='TEXT'
            if ftype=='Integer':
                ftojoin_valid[fname] = 'LONG'
            else:
                ftojoin_valid[fname]=ftype.upper()

        if lenf in field_dict:
            ftojoin_valid[lenf] = 'DOUBLE'

        for fname, ftype in ftojoin_valid.items():
            if fname not in [f.name for f in arcpy.ListFields(out_tab)]:
                arcpy.AddField_management(out_tab, field_name=fname, field_type=ftype)

        fcontent_dict = {row[0]: row[1:] for row in
                         arcpy.da.SearchCursor(in_tab, ["OID@"] + list(ftojoin_valid.keys())
                                               )}

        with arcpy.da.UpdateCursor(out_tab, (["FID_{}".format(os.path.basename(in_tab))] +
                                            list(ftojoin_valid.keys()))) as cursor:
            for row in cursor:
                if row[0] in fcontent_dict:
                    row[1:] = list(fcontent_dict[row[0]])
                    cursor.updateRow(row)

    join_field_fromdict(in_tab = in_join_ft_buf,
                        out_tab = out_ft,
                        field_dict = ftojoin_dict
                        )

    join_field_fromdict(in_tab = in_target_ft,
                        out_tab = out_ft,
                        field_dict = target_ftojoin_dict
                        )

if not arcpy.Exists(ddtnet_bdtopo_inters):
    print("Intersecting bdtopo")
    inters_linetoline(
        in_target_ft=ce_net,
        in_join_ft=bdtopo2015_fr,
        out_ft=ddtnet_bdtopo_inters,
        tempgdb=tempgdb,
        in_join_ft_suffix='bdtopo',
        in_join_ft_fieldstojoin=['ID'],
        in_target_ft_fieldstojoin=['UID_CE'],
        tolerance='5 Meters'
    )
if not arcpy.Exists(ddtnet_bdtopo_inters_tab):
    CopyRows_pd(in_table=ddtnet_bdtopo_inters,
                out_table=ddtnet_bdtopo_inters_tab,
                fields_to_copy={'UID_CE': 'UID_CE',
                                'geom_Length': 'length_inters_bdtopo',
                                'ID_bdtopo': 'ID_bdtopo',
                                'length_bdtopo': 'length_bdtopo'})

if not arcpy.Exists(ddtnet_carthage_inters):
    print("Intersecting carthage")
    inters_linetoline(in_target_ft=ce_net,
                      in_join_ft=bdcarthage,
                      out_ft=ddtnet_carthage_inters,
                      tempgdb=tempgdb,
                      in_join_ft_suffix='carthage',
                      in_join_ft_fieldstojoin=['CODE_HYDRO', "ID_BDCARTH"],
                      in_target_ft_fieldstojoin=['UID_CE'],
                      tolerance='5 Meters'
                      )
if not arcpy.Exists(ddtnet_carthage_inters_tab):
    CopyRows_pd(in_table=ddtnet_carthage_inters,
                out_table=ddtnet_carthage_inters_tab,
                fields_to_copy={'UID_CE': 'UID_CE',
                                'geom_Length': 'length_inters_carthage',
                                'ID_BDCARTH_carthage': 'ID_carthage',
                                'length_carthage': 'length_carthage'})

#------------- IDENTIFY LONGITUDINAL DISCONNECTIONS IN DDT NETWORK  ----------------------------------------------------
#Start trying with BD Alti
#(Maybe smooth elevation)
#Merge all lines between confluences
#Extract elevation at very node along line
#Run a regression across nodes
#Establish direction of line based on regression

#Then identify all lines with a dangle point that is also a start point
#Check how many have a dangle point that is an end point - random sample

#Check adjacent lines:
#if end point is connected to an endpoint (or multiple endpoints), and start point is connected to a start point



#------------- INTERSECT ALL NETWORKS WITH UNITS OF ANALYSIS -----------------------------------------------------------
for net in [ce_net, bdtopo2015_fr, bcae_fr, bdcarthage, rht]:
    root_name = os.path.split(os.path.splitext(net)[0])[1]
    out_net = os.path.join(pregdb, "{}_bvinters".format(root_name))
    out_tab = os.path.join(resdir, "{}_bvinters.csv".format(root_name))
    if (not arcpy.Exists(out_net)) or overwrite:
        print("Processing {}...".format(out_net))
        arcpy.analysis.Intersect([net, cats_hybasdeps], out_feature_class=out_net, join_attributes="ALL")
    if (not arcpy.Exists(out_tab)) or overwrite:
        print("Exporting {}...".format(out_tab))
        arcpy.CopyRows_management(out_net, out_tab)


######Extra stuff
    # # create the field mapping object
    # fms_spjoin = arcpy.FieldMappings()
    # # populate the field mapping object with the fields from both feature classes
    # fms_spjoin.addTable(in_target_ft)
    # fms_spjoin.addTable(in_join_ft_buf)
    # # loop through the field names to sum
    # for field in fms_spjoin.fields:
    #     if field.name not in (
    #             [f.name for f in arcpy.ListFields(in_target_ft)] + [lenf] +
    #             ["{0}_{1}".format(ftojoin, in_join_ft_suffix) for ftojoin in in_join_ft_fieldstojoin]):
    #         fms_spjoin.removeFieldMap(fms_spjoin.findFieldMapIndex(field.name))
    #
    # arcpy.SpatialJoin_analysis(target_features=in_target_ft,
    #                            join_features=in_join_ft_buf,
    #                            out_feature_class=out_ft,
    #                            join_operation='JOIN_ONE_TO_ONE',
    #                            join_type='KEEP_ALL',
    #                            field_mapping=fms_spjoin,
    #                            match_option='WITHIN'
    #                            )
