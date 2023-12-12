import os.path

import arcpy.analysis

from setup_classement import *
from sklearn.linear_model import HuberRegressor, LinearRegression

overwrite=False

anci_dir = os.path.join(datdir, "données_auxiliaires")#Ancillary data directory

#Preprocessing geodatabase
pregdb = os.path.join(resdir, "preprocessing_ancillary_data.gdb")

#Scratch gdb
tempgdb = os.path.join(resdir, "scratch.gdb")
if not arcpy.Exists(tempgdb):
    arcpy.CreateFileGDB_management(out_folder_path=os.path.split(tempgdb)[0],
                                   out_name=os.path.split(tempgdb)[1])

#DDT network
outputs_gdb = os.path.join(resdir, 'analysis_outputs.gdb')
ddt_net = os.path.join(outputs_gdb, 'carto_loi_eau_fr')

#ddt net after linkage to BD topo and BD carthage for auxilliary information (intermittence status, name, etc.)
ddt_net_refnetsattri_tab = os.path.join(resdir, 'carto_loi_eau_refnetsattris.csv')

#BDTOPO
bdtopo2015_fr = os.path.join(pregdb, 'bdtopo2015_fr')

#Large hydrographic units
BH = os.path.join(anci_dir, 'bdtopage', 'BassinHydrographique_FXX.shp')

# BD alti
bdalti_mosaic = os.path.join(pregdb, "bdalti_25m_mosaic")
#------------------- OUTPUTS ------------------------------------------------------------------------------------------
ddt_net_endpts = os.path.join(pregdb, 'carto_loi_eau_fr_endpts')
ddt_net_endpts_inters = os.path.join(pregdb, 'carto_loi_eau_fr_endpts_inters')
ddt_net_isolatedsegs = os.path.join(pregdb, 'carto_loi_eau_fr_isolated')

ddt_net_artif = os.path.join(pregdb, 'bdtopo_artif')
ddt_net_noartif = os.path.join(pregdb, 'bdtopo_noartif')
ddt_net_noartif_bh = os.path.join(pregdb, 'bdtopo_noartif_bh')


bdtopo_endpts = os.path.join(pregdb, 'bdtopo_fr_endpts')
bdtopo_endpts_inters = os.path.join(pregdb, 'bdtopo_fr_endpts_inters')
bdtopo_isolatedsegs = os.path.join(pregdb, 'bdtopo_fr_isolated')

bdtopo_artif = os.path.join(pregdb, 'bdtopo_artif')
bdtopo_noartif = os.path.join(pregdb, 'bdtopo_noartif')
bdtopo_noartif_bh = os.path.join(pregdb, 'bdtopo_noartif_bh')

#------------- IDENTIFY DISCONNECTIONS --------------------------------------------------------------------------------
#Convert lines to end points
if not arcpy.Exists(ddt_net_endpts):
    arcpy.FeatureVerticesToPoints_management(in_features=ddt_net,
                                             out_feature_class=ddt_net_endpts,
                                             point_location="BOTH_ENDS")
#Intersect end points
if not arcpy.Exits(ddt_net_endpts_inters):
    arcpy.analysis.Intersect(in_features=[ddt_net_endpts,ddt_net_endpts],
                             out_feature_class=ddt_net_endpts_inters,
                             join_attributes='ALL')

#Identify isolated river segmenets (i.e., segmenets surrounded by segments with another category. e.g., nce surrounded by ce)
if not arcpy.Exists(ddt_net_isolatedsegs):
    #Identify lines with one connected end point and one connected start point with dictionary
    lines_inters_dict = {}
    #[f.name for f in arcpy.ListFields(ddt_net_endpts_inters)]
    with arcpy.da.SearchCursor(ddt_net_endpts_inters,
                               ['OID@', 'FID_carto_loi_eau_fr_endpts', 'UID_CE', 'UID_CE_1',
                                'type_stand', 'type_stand_1']) as cursor:
        for row in cursor:
            if row[2] != row[3]:
                lines_inters_dict[row[0]] = row[1:]
    lines_inters_pd = pd.DataFrame.from_dict(lines_inters_dict, orient='index')
    del lines_inters_dict
    lines_inters_pd.columns = ['FID_carto_loi_eau_fr_endpts', 'UID_CE', 'UID_CE_1', 'type_stand', 'type_stand_1']
    #Count number of intersecting points per end point
    lines_inters_pd_pts_interscount = lines_inters_pd.groupby(by=['FID_carto_loi_eau_fr_endpts', 'UID_CE']).size()
    #Count number of end points with intersecting points per segment
    lines_inters_pd_ptswinters = lines_inters_pd_pts_interscount.groupby(by='UID_CE').size()
    #Examine only those with two end points intersecting with other points (excluding first and last order streams)
    #And segments with a total of intersecting points of 3 and under, including those either without confluence
    #with upstream confluence or downstream bifurcation
    sel1 = (lines_inters_pd_pts_interscount.groupby('UID_CE').sum() <3)
    sel2 = (lines_inters_pd_ptswinters == 2)

    simple_lines_pd = lines_inters_pd.loc[
                          (lines_inters_pd.UID_CE.isin(list(sel1[sel1].index))) &
                          (lines_inters_pd.UID_CE.isin(list(sel2[sel2].index))),
                          :]

    #For those lines, select those whose type of both ends differ from their own
    simple_lines_pd['dif_i'] = np.where(simple_lines_pd.type_stand!=simple_lines_pd.type_stand_1, 0, 1)
    sel3 = (simple_lines_pd.groupby('UID_CE').sum('dif_i').dif_i == 0)
    simple_lines_pd_isolated_ID = sel3[sel3].index

    #Copy ddt_net the delete all features but those to inspect, add fields with types of the surrounding segments
    arcpy.CopyFeatures_management(ddt_net, ddt_net_isolatedsegs)
    with arcpy.da.UpdateCursor(ddt_net_isolatedsegs, ["UID_CE"]) as cursor:
        for row in cursor:
            if row[0] not in simple_lines_pd_isolated_ID:
                cursor.deleteRow()

#Manually go through them to check which ones are "unjustified"
#If line segment only partly contained in standing water body
6071, 7628, 8646, 9834, 9916,


#Check percentage of NCE for which surrounding ones are CE

#------------- TRY TO ENHANCE DDT AND BDTOPO NETWORKS WITH FLOW DIRECTION --------------------------------------------------------
#Start trying with BD Alti
#(Maybe smooth elevation)
##(maybe also try removing all ARTIF from BD Topo, but inconsistently labeled)
def get_lr_coef(df, x, y):
    return (LinearRegression().fit(df[[x]], df[[y]]).coef_[0][0])
def identify_strahler1(in_net, in_dem, suffix, out_gdb, out_net, prefix='carto_loi_eau_fr'):
    ddt_net_dissolved = os.path.join(out_gdb, '{0}_noartif_dissolved_{1}'.format(prefix, suffix))
    ddt_net_conflupts = os.path.join(out_gdb, '{0}_noartif_conflupts_{1}'.format(prefix, suffix))
    ddt_net_conflusplit = os.path.join(out_gdb, '{0}_noartif_conflusplit_{1}'.format(prefix, suffix))
    ddt_net_nodes25m = os.path.join(out_gdb, '{0}_nodes25m_{1}'.format(prefix, suffix))
    ddt_nodes_elv = os.path.join(out_gdb, '{0}_nodes25m_elv_{1}'.format(prefix, suffix))
    ddt_net_conflusplit_directed = os.path.join(out_gdb, '{0}_noartif_conflusplit_directed_{1}'.format(prefix, suffix))
    ddt_net_conflusplit_dangle = os.path.join(out_gdb, '{0}_noartif_conflusplit_dangle_{1}'.format(prefix, suffix))
    ddt_net_conflusplit_start = os.path.join(out_gdb, '{0}_noartif_conflusplit_start_{1}'.format(prefix, suffix))
    ddt_net_conflusplit_danglestart_inters = os.path.join(out_gdb, '{0}_noartif_conflusplit_danglestartinters_{1}'.format(prefix, suffix))

    if not arcpy.Exists(ddt_net_dissolved):
        print('Dissolving the original network...')
        arcpy.Dissolve_management(in_net,
                                  out_feature_class=ddt_net_dissolved,
                                  multi_part='SINGLE_PART',
                                  unsplit_lines='UNSPLIT_LINES')

    if not arcpy.Exists(ddt_net_conflupts):
        print('Get confluence points to make sure that all lines are split at confluences...')
        arcpy.Intersect_analysis(in_features=ddt_net_dissolved,
                                 out_feature_class=ddt_net_conflupts,
                                 join_attributes='ONLY_FID',
                                 output_type='POINT')

    if not arcpy.Exists(ddt_net_conflusplit):
        print('Split dissolved network at confluence points...')
        arcpy.SplitLineAtPoint_management(in_features=ddt_net_dissolved,
                                          point_features=ddt_net_conflupts,
                                          out_feature_class=ddt_net_conflusplit,
                                          search_radius=0.1)

    if not arcpy.Exists(ddt_net_nodes25m):
        print('Generate points along lines...')
        arcpy.GeneratePointsAlongLines_management(Input_Features=ddt_net_conflusplit,
                                                  Output_Feature_Class=ddt_net_nodes25m,
                                                  Point_Placement="DISTANCE",
                                                  Distance="25 meters",
                                                  Include_End_Points="END_POINTS")
                                                  #Add_Chainage_Fields="ADD_CHAINAGE")

    if not arcpy.Exists(ddt_nodes_elv):
        print('Get elevation for each point...')
        Sample(in_rasters=in_dem,
               in_location_data=ddt_net_nodes25m,
               out_table=ddt_nodes_elv,
               resampling_type='NEAREST',
               unique_id_field=arcpy.Describe(ddt_net_nodes25m).OIDFieldName,
               statistics_type='MEDIAN',
               layout='ROW_WISE',
               generate_feature_class='TABLE')

        colstoread = [f.name for f in arcpy.ListFields(ddt_nodes_elv) if f.type!="Geometry"] #List the fields you want to include. I want all columns except the geometry
        #Read elevation table and join OBJECTID of line associated with the node
        ddt_nodes_elv_pd = pd.DataFrame(data=arcpy.da.SearchCursor(ddt_nodes_elv, colstoread), columns=colstoread).\
            merge(pd.DataFrame(data=arcpy.da.SearchCursor(ddt_net_nodes25m, ['OID@', 'ORIG_FID']),
                                        columns=[os.path.basename(ddt_net_nodes25m), 'ORIG_FID']),
                 on=os.path.basename(ddt_net_nodes25m))


        print("Compute linear regression for each line...")
        coefs_dict = (ddt_nodes_elv_pd.loc[~ddt_nodes_elv_pd.bdalti_25m_mosaic_Band_1.isna()].groupby('ORIG_FID').
                      apply(get_lr_coef, x='OBJECTID', y='bdalti_25m_mosaic_Band_1').
                      to_dict())

        elv_coef_fname = 'elv_coef'
        if not elv_coef_fname in [f.name for f in arcpy.ListFields(ddt_net_conflusplit)]:
            arcpy.AddField_management(ddt_net_conflusplit, elv_coef_fname, 'DOUBLE')
        with arcpy.da.UpdateCursor(ddt_net_conflusplit, ['OID@', elv_coef_fname]) as cursor:
            for row in cursor:
                if row[0] in coefs_dict:
                    row[1] = coefs_dict[row[0]]
                    cursor.updateRow(row)

    if not arcpy.Exists(ddt_net_conflusplit_directed):
        print("Flip lines...")
        arcpy.CopyFeatures_management(ddt_net_conflusplit, ddt_net_conflusplit_directed)
        arcpy.MakeFeatureLayer_management(ddt_net_conflusplit_directed, 'segs_to_flip_lyr',
                                          where_clause="elv_coef > 0.1")
        arcpy.edit.FlipLine(in_features='segs_to_flip_lyr')

    print("Assign first order to lines with dangle point and start point coinciding")
    if not arcpy.Exists(ddt_net_conflusplit_dangle):
        arcpy.FeatureVerticesToPoints_management(in_features=ddt_net_conflusplit_directed,
                                                 out_feature_class=ddt_net_conflusplit_dangle,
                                                 point_location='DANGLE')

    if not arcpy.Exists(ddt_net_conflusplit_start):
        arcpy.FeatureVerticesToPoints_management(in_features=ddt_net_conflusplit_directed,
                                                 out_feature_class=ddt_net_conflusplit_start,
                                                 point_location='START')

    if not arcpy.Exists(ddt_net_conflusplit_danglestart_inters):
        arcpy.Intersect_analysis(in_features=[ddt_net_conflusplit_dangle, ddt_net_conflusplit_start],
                                 out_feature_class=ddt_net_conflusplit_danglestart_inters,
                                 join_attributes='ALL')

    first_order_list = [row[0] for row in arcpy.da.SearchCursor(ddt_net_conflusplit_danglestart_inters, ['ORIG_FID'])]
    if 'strahler' not in [f.name for f in arcpy.ListFields(ddt_net_conflusplit_directed)]:
        arcpy.AddField_management(ddt_net_conflusplit_directed, 'strahler', 'SHORT')
    with arcpy.da.UpdateCursor(ddt_net_conflusplit_directed, ['OID@', 'strahler']) as cursor:
        for row in cursor:
            if row[0] in first_order_list:
                row[1] = 1
                cursor.updateRow(row)

    #Remove lines with dangle points under 10 m
    # arcpy.MakeFeatureLayer_management(in_net, 'ddt_confusplit_sublyr',
    #                                   where_clause="orig_layer='D73_Savoie'")
    if not arcpy.Exists(out_net):
        print("Re-assign Strahler order attribute to original network...")
        arcpy.SpatialJoin_analysis(target_features=in_net,
                                   join_features=ddt_net_conflusplit_directed,
                                   out_feature_class=out_net,
                                   join_operation='JOIN_ONE_TO_ONE',
                                   join_type='KEEP_ALL',
                                   match_option='LARGEST_OVERLAP'
                                   )

#Remove all truly artificial watercourses that may disrupt topology
if not arcpy.Exists(ddt_net_noartif):
    arcpy.CopyFeatures_management(ddt_net, ddt_net_artif)
    arcpy.CopyFeatures_management(ddt_net, ddt_net_noartif)

    ddt_net_refnetsattri_pd = pd.read_csv(ddt_net_refnetsattri_tab)
    regex = r'\b)|(\b'.join([r"(\bcanal", "foss[eé]", "roubine", "craste", "d[ée]rivation",
                      "bief", "aber", r"hydraulique\b)"])
    ddt_net_refnetsattri_pd['NOM_bdtopo_orig'] = ddt_net_refnetsattri_pd['NOM_bdtopo_orig'].str.lower()
    ddt_net_refnetsattri_pd['TOPONYME1_carthage'] = ddt_net_refnetsattri_pd['TOPONYME1_carthage'].str.lower()
    artif_ids = list(
        ddt_net_refnetsattri_pd.loc[(ddt_net_refnetsattri_pd['NOM_bdtopo_orig'].str.contains(regex)|
                                     ddt_net_refnetsattri_pd['TOPONYME1_carthage'].str.contains(regex)), 'UID_CE']
    )

    with arcpy.da.UpdateCursor(ddt_net_artif, ['UID_CE']) as cursor:
        for row in cursor:
            if not row[0] in artif_ids:
                cursor.deleteRow()

    with arcpy.da.UpdateCursor(ddt_net_noartif, ['UID_CE']) as cursor:
        for row in cursor:
            if row[0] in artif_ids:
                cursor.deleteRow()

#Spatial join to large basin region
if not arcpy.Exists(ddt_net_noartif_bh):
    arcpy.SpatialJoin_analysis(ddt_net_noartif,
                               join_features=BH,
                               out_feature_class=ddt_net_noartif_bh,
                               join_operation='JOIN_ONE_TO_ONE',
                               join_type='KEEP_ALL',
                               match_option='LARGEST_OVERLAP')

#Merge all lines between confluences
bh_numset = {row[0] for row in arcpy.da.SearchCursor(ddt_net_noartif_bh, 'CdBH')}
for bh_num in bh_numset:
    if bh_num is not None:
        print(bh_num)
        arcpy.MakeFeatureLayer_management(ddt_net_noartif_bh, 'ddt_net_noartif_sub',
                                          where_clause="CdBH='{}'".format(bh_num))
        ddt_net_noartif_wstrahler1 = os.path.join(pregdb, 'carto_loi_eau_noartif_strahler1_{}'.format(bh_num))
        identify_strahler1(in_net='ddt_net_noartif_sub',
                           in_dem=bdalti_mosaic,
                           suffix=bh_num,
                           out_gdb=pregdb,
                           out_net=ddt_net_noartif_wstrahler1)

################################ DO THE SAME FOR BDTOPO ################################################################
# Remove all truly artificial watercourses that may disrupt topology
if not arcpy.Exists(bdtopo_noartif):
    arcpy.CopyFeatures_management(bdtopo2015_fr, bdtopo_artif)
    arcpy.CopyFeatures_management(bdtopo2015_fr, bdtopo_noartif)

    bdtopo_attri_pd = pd.DataFrame(data=arcpy.da.SearchCursor(bdtopo2015_fr, ['ID', 'NOM']),
                                               columns=['ID', 'NOM'])
    regex = r'\b)|(\b'.join([r"(\bcanal", "foss[eé]", "roubine", "craste", "d[ée]rivation",
                             "bief", "aber", r"hydraulique\b)"])
    bdtopo_attri_pd['NOM'] = bdtopo_attri_pd['NOM'].str.lower()
    artif_ids = list(
        bdtopo_attri_pd.loc[(bdtopo_attri_pd['NOM'].str.contains(regex)), 'ID']
    )

    with arcpy.da.UpdateCursor(bdtopo_artif, ['ID']) as cursor:
        for row in cursor:
            if not row[0] in artif_ids:
                cursor.deleteRow()

    with arcpy.da.UpdateCursor(bdtopo_noartif, ['ID']) as cursor:
        for row in cursor:
            if row[0] in artif_ids:
                cursor.deleteRow()

# Spatial join to large basin region
if not arcpy.Exists(bdtopo_noartif_bh):
    arcpy.SpatialJoin_analysis(bdtopo_noartif,
                               join_features=BH,
                               out_feature_class=bdtopo_noartif_bh,
                               join_operation='JOIN_ONE_TO_ONE',
                               join_type='KEEP_ALL',
                               match_option='LARGEST_OVERLAP')

# Merge all lines between confluences
bh_numset = {row[0] for row in arcpy.da.SearchCursor(bdtopo_noartif_bh, 'CdBH')}
for bh_num in bh_numset:
    if bh_num is not None:
        print(bh_num)
        arcpy.MakeFeatureLayer_management(bdtopo_noartif_bh, 'bdtopo_noartif_sub',
                                          where_clause="CdBH='{}'".format(bh_num))
        bdtopo_noartif_wstrahler1 = os.path.join(pregdb, 'bdtopo_noartif_strahler1_{}'.format(bh_num))
        identify_strahler1(in_net='bdtopo_noartif_sub',
                           in_dem=bdalti_mosaic,
                           suffix=bh_num,
                           out_gdb=pregdb,
                           out_net=bdtopo_noartif_wstrahler1,
                           prefix='bdtopo')
#Check how many have a dangle point that is an end point - random sample

#Check adjacent lines:
#if end point is connected to an endpoint (or multiple endpoints), and start point is connected to a start point