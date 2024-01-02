import os
import arcpy.analysis
import collections

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

ddt_net_artif = os.path.join(pregdb, 'carto_loi_eau_fr_artif')
ddt_net_noartif = os.path.join(pregdb, 'carto_loi_eau_fr_noartif')
ddt_net_noartif_bh = os.path.join(pregdb, 'carto_loi_eau_fr_noartif_bh')

bdtopo_endpts = os.path.join(pregdb, 'bdtopo_fr_endpts')
bdtopo_endpts_inters = os.path.join(pregdb, 'bdtopo_fr_endpts_inters')
bdtopo_isolatedsegs = os.path.join(pregdb, 'bdtopo_fr_isolated')

bdtopo_artif = os.path.join(pregdb, 'bdtopo_artif')
bdtopo_noartif = os.path.join(pregdb, 'bdtopo_noartif')
bdtopo_noartif_bh = os.path.join(pregdb, 'bdtopo_noartif_bh')

ddt_net_wstrahler = os.path.join(pregdb, 'carto_loi_eau_noartif_strahler_fr')
ddt_net_wstrahler_tab = os.path.join(resdir, 'carto_loi_eau_noartif_strahler_fr.csv')

bdtopo_wstrahler = os.path.join(pregdb, 'bdtopo_noartif_strahler_fr')
bdtopo_wstrahler_tab = os.path.join(resdir, 'bdtopo_noartif_strahler_fr.csv')

#------------- IDENTIFY DISCONNECTIONS --------------------------------------------------------------------------------
#Convert lines to end points
if not arcpy.Exists(ddt_net_endpts):
    arcpy.FeatureVerticesToPoints_management(in_features=ddt_net,
                                             out_feature_class=ddt_net_endpts,
                                             point_location="BOTH_ENDS")
#Intersect end points
if not arcpy.Exists(ddt_net_endpts_inters):
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
#6071, 7628, 8646, 9834, 9916,


#Check percentage of NCE for which surrounding ones are CE

#------------- TRY TO ENHANCE DDT AND BDTOPO NETWORKS WITH FLOW DIRECTION --------------------------------------------------------
#Start trying with BD Alti
#(Maybe smooth elevation)
##(maybe also try removing all ARTIF from BD Topo, but inconsistently labeled)
def get_lr_coef(df, x, y):
    return (LinearRegression().fit(df[[x]], df[[y]]).coef_[0][0])

def assign_strahler_splitline(in_lines, out_gdb, prefix, suffix, in_strahler):
        for vertex_direction in ['START', 'END']:
            qualifying_endpts_selfinters_ids = {1}
            while len(qualifying_endpts_selfinters_ids) > 0:
                print("Assign same Strahler order to lines connected by one end exclusively to another line of order "
                      "{0} based on {1} points".format(in_strahler, vertex_direction))

                in_lines_oneendpts = os.path.join(
                    out_gdb,'{0}_noartif_conflusplit_directed_{1}_{2}'.format(prefix, vertex_direction, suffix))
                arcpy.FeatureVerticesToPoints_management(in_features=in_lines,
                                                         out_feature_class=in_lines_oneendpts,
                                                         point_location=vertex_direction)

                in_lines_bothendpts = os.path.join(
                    out_gdb,'{0}_noartif_conflusplit_directed_bothendpts_{1}'.format(prefix, suffix))
                arcpy.FeatureVerticesToPoints_management(in_features=in_lines,
                                                         out_feature_class=in_lines_bothendpts,
                                                         point_location='BOTH_ENDS')

                in_lines_bothendpts_selfinters = os.path.join(
                    out_gdb,'{0}_noartif_conflusplit_directed_bothendpts_selfinters_{1}'.format(prefix, suffix))
                arcpy.Intersect_analysis(in_features=[in_lines_oneendpts,
                                                      in_lines_bothendpts],
                                         out_feature_class=in_lines_bothendpts_selfinters,
                                         join_attributes='ALL')

                split_order1_pts_selfinters_1 = defaultdict(list)
                split_order1_pts_selfinters_undefined = defaultdict(list)
                with (arcpy.da.SearchCursor(in_lines_bothendpts_selfinters,
                                           ['ORIG_FID', 'ORIG_FID_1', 'strahler', 'strahler_1']) as cursor):
                    for row in cursor:
                        if ((row[0] != row[1]) and (row[1] not in split_order1_pts_selfinters_1[row[0]])
                            and (row[2] is None) and (row[3]==in_strahler)):
                            #If intersecting point is from different line, and other line is not already in dict
                            #and strahler order is undefined for this line and one for this other line
                            split_order1_pts_selfinters_1[row[0]].append(row[1])
                        if ((row[0] != row[1]) and (row[1] not in split_order1_pts_selfinters_1[row[0]])
                                and (row[2] is None) and (row[3] is None)):
                            split_order1_pts_selfinters_undefined[row[0]].append(row[1])

                #Identify those "other 1st order lines" which intersect with only one undefined line
                #(to avoid first order streams intersecting already-second order streams)
                split_order1_pts_selfinters_rev = defaultdict(list)
                for k, v in split_order1_pts_selfinters_1.items():
                    for i in v:
                        split_order1_pts_selfinters_rev[i].append(k)
                endpts_selfinters_unique_match_ids = {k for k,v in split_order1_pts_selfinters_rev.items() if len(v)==1}

                qualifying_endpts_selfinters_ids = {k for k, v in split_order1_pts_selfinters_1.items()
                                                    if (len(v) == 1) and #Make sure there aren't two first-order streams intersecting the main line
                                                    (len(split_order1_pts_selfinters_undefined[k])==0) and #And make sure that the undefined line does not intersect with another undefined line on that end
                                                    v[0] in endpts_selfinters_unique_match_ids}

                with arcpy.da.UpdateCursor(in_lines, ['OID@', 'strahler']) as cursor:
                    for row in cursor:
                        if row[0] in qualifying_endpts_selfinters_ids:
                            row[1] = in_strahler
                            cursor.updateRow(row)

# in_lines = ddt_net_conflusplit_directed
# in_loop_ids_tab_path = loop_ids_tab_path
# strahler_ini = 1
def iterate_strahler(in_lines, out_gdb, prefix, suffix, in_loop_ids_tab_path, strahler_ini):
    print('Assigning strahler order {}'.format(strahler_ini+1))
    for vertex_direction in ['START', 'END']:
        in_lines_oneendpts = os.path.join(
            out_gdb, '{0}_noartif_conflusplit_directed_{1}_{2}'.format(prefix, vertex_direction, suffix))
        arcpy.FeatureVerticesToPoints_management(in_features=in_lines,
                                                 out_feature_class=in_lines_oneendpts,
                                                 point_location=vertex_direction)

        in_lines_bothendpts = os.path.join(
            out_gdb, '{0}_noartif_conflusplit_directed_bothendpts_{1}'.format(prefix, suffix))
        arcpy.FeatureVerticesToPoints_management(in_features=in_lines,
                                                 out_feature_class=in_lines_bothendpts,
                                                 point_location='BOTH_ENDS')

        in_lines_bothendpts_selfinters = os.path.join(
            out_gdb, '{0}_noartif_conflusplit_directed_bothendpts_selfinters_{1}'.format(prefix, suffix))
        arcpy.Intersect_analysis(in_features=[in_lines_oneendpts,
                                              in_lines_bothendpts],
                                 out_feature_class=in_lines_bothendpts_selfinters,
                                 join_attributes='ALL')

        #Unique IDs of loops
        loop_ids_list = set().union(*
            [{v_2 for v_2 in v_1.split('_')} for v_1 in pd.read_csv(in_loop_ids_tab_path)['OBJECTID']]
                                    )
        #Get the set of lines intersecting each line with strahler_ini, excluding loops
        split_order1_pts_selfinters = defaultdict(set)
        with (arcpy.da.SearchCursor(in_lines_bothendpts_selfinters,
                                    ['ORIG_FID', 'ORIG_FID_1', 'strahler', 'strahler_1']) as cursor):
            for row in cursor:
                if ((row[0] != row[1]) and (row[1] not in split_order1_pts_selfinters[row[0]])
                        and (row[2] is None) and (row[3] ==strahler_ini)) and (not row[1] in loop_ids_list):
                    # If intersecting point is from different line, and other line is not already in dict
                    # and strahler order is undefined for this line and strahler_ini for this other line
                    # and this other line is not a flagged loop line
                    split_order1_pts_selfinters[row[0]].add(row[1])

        #If intersecting at least two lines with strahler_ini, assign strahler_ini + 1
        with arcpy.da.UpdateCursor(in_lines, ['OID@', 'strahler']) as cursor:
            for row in cursor:
                if row[0] in split_order1_pts_selfinters:
                    # if len(split_order1_pts_selfinters[row[0]]) == 1:
                    #     row[1] = strahler_ini
                    if len(split_order1_pts_selfinters[row[0]]) > 1:
                        row[1] = strahler_ini + 1
                    cursor.updateRow(row)

#Deal with simple loops (two lines sharing start and end points) ---------------------------
#Make a subset of undefined lines
#Unsplit lines
#Identify loops
#Find lines with an assigned strahler order that connect with two undefined on one end - dictionary key-value
#Reverse dictionary to identify loops with only one undefined line connected to it
#Exclude those with more than one undefined line connected to it
#Exclude lines in those loops that have a dangle point -> assign higher strahler order of line in the loop without dangle and line connected to that loop
#Else (those with one undefined line connected to it and no dangle points): assign higher connected strahler order
#Assign that strahler order to undefined connected line
def identify_strahler(in_net, in_dem, suffix, out_gdb, out_net, prefix):
    net_dissolved = os.path.join(out_gdb, '{0}_noartif_dissolved_{1}'.format(prefix, suffix))
    net_conflupts = os.path.join(out_gdb, '{0}_noartif_conflupts_{1}'.format(prefix, suffix))
    net_conflusplit = os.path.join(out_gdb, '{0}_noartif_conflusplit_{1}'.format(prefix, suffix))
    net_nodes25m = os.path.join(out_gdb, '{0}_nodes25m_{1}'.format(prefix, suffix))
    net_nodes_elv = os.path.join(out_gdb, '{0}_nodes25m_elv_{1}'.format(prefix, suffix))
    net_conflusplit_directed_ini = os.path.join(out_gdb, '{0}_noartif_conflusplit_directed_ini_{1}'.format(prefix, suffix))
    net_conflusplit_dangle = os.path.join(out_gdb, '{0}_noartif_conflusplit_dangle_{1}'.format(prefix, suffix))
    net_conflusplit_start = os.path.join(out_gdb, '{0}_noartif_conflusplit_start_{1}'.format(prefix, suffix))
    net_conflusplit_danglestart_inters = os.path.join(out_gdb, '{0}_noartif_conflusplit_danglestartinters_{1}'.format(prefix, suffix))
    net_conflusplit_directed = os.path.join(out_gdb, '{0}_noartif_conflusplit_directed_{1}'.format(prefix, suffix))

    if not arcpy.Exists(net_dissolved):
        print('Dissolving the original network...')
        arcpy.Dissolve_management(in_net,
                                  out_feature_class=net_dissolved,
                                  multi_part='SINGLE_PART',
                                  unsplit_lines='UNSPLIT_LINES')

    if not arcpy.Exists(net_conflupts):
        print('Get confluence points to make sure that all lines are split at confluences...')
        arcpy.Intersect_analysis(in_features=net_dissolved,
                                 out_feature_class=net_conflupts,
                                 join_attributes='ONLY_FID',
                                 output_type='POINT')

    if not arcpy.Exists(net_conflusplit):
        print('Split dissolved network at confluence points...')
        arcpy.SplitLineAtPoint_management(in_features=net_dissolved,
                                          point_features=net_conflupts,
                                          out_feature_class=net_conflusplit,
                                          search_radius=0.1)

    if not arcpy.Exists(net_nodes25m):
        print('Generate points along lines...')
        arcpy.GeneratePointsAlongLines_management(Input_Features=net_conflusplit,
                                                  Output_Feature_Class=net_nodes25m,
                                                  Point_Placement="DISTANCE",
                                                  Distance="25 meters",
                                                  Include_End_Points="END_POINTS")
                                                  #Add_Chainage_Fields="ADD_CHAINAGE")

    if not arcpy.Exists(net_nodes_elv):
        print('Get elevation for each point...')
        Sample(in_rasters=in_dem,
               in_location_data=net_nodes25m,
               out_table=net_nodes_elv,
               resampling_type='NEAREST',
               unique_id_field=arcpy.Describe(net_nodes25m).OIDFieldName,
               statistics_type='MEDIAN',
               layout='ROW_WISE',
               generate_feature_class='TABLE')

        colstoread = [f.name for f in arcpy.ListFields(net_nodes_elv) if f.type!="Geometry"] #List the fields you want to include. I want all columns except the geometry
        #Read elevation table and join OBJECTID of line associated with the node
        net_nodes_elv_pd = pd.DataFrame(data=arcpy.da.SearchCursor(net_nodes_elv, colstoread), columns=colstoread).\
            merge(pd.DataFrame(data=arcpy.da.SearchCursor(net_nodes25m, ['OID@', 'ORIG_FID']),
                                        columns=[os.path.basename(net_nodes25m), 'ORIG_FID']),
                 on=os.path.basename(net_nodes25m))


        print("Compute linear regression for each line...")
        coefs_dict = (net_nodes_elv_pd.loc[~net_nodes_elv_pd.bdalti_25m_mosaic_Band_1.isna()].groupby('ORIG_FID').
                      apply(get_lr_coef, x='OBJECTID', y='bdalti_25m_mosaic_Band_1').
                      to_dict())

        elv_coef_fname = 'elv_coef'
        if not elv_coef_fname in [f.name for f in arcpy.ListFields(net_conflusplit)]:
            arcpy.AddField_management(net_conflusplit, elv_coef_fname, 'DOUBLE')
        with arcpy.da.UpdateCursor(net_conflusplit, ['OID@', elv_coef_fname]) as cursor:
            for row in cursor:
                if row[0] in coefs_dict:
                    row[1] = coefs_dict[row[0]]
                    cursor.updateRow(row)

    #Identify loops ------------------------------------------------------------------------------------------------------
    loop_ids_tab_path = os.path.join(resdir, '{0}_loop_ids_tab_{1}.csv'.format(prefix, suffix))
    if not arcpy.Exists(loop_ids_tab_path):
        print("Identifying simple loops...")
        net_conflusplit_unsplit = os.path.join(out_gdb, '{0}_noartif_conflusplit_unsplit_{1}'.format(prefix, suffix))
        arcpy.UnsplitLine_management(in_features=net_conflusplit,
                                     out_feature_class=net_conflusplit_unsplit,
                                     statistics_fields=[[arcpy.Describe(net_conflusplit).OIDFieldName,
                                                          'CONCATENATE']],
                                     concatenation_separator='_')

        net_conflusplit_unsplit_bothendpts = os.path.join(
            out_gdb, '{0}_noartif_conflusplit_unsplit_bothendpts_{1}'.format(prefix, suffix))
        arcpy.FeatureVerticesToPoints_management(in_features=net_conflusplit_unsplit,
                                                 out_feature_class=net_conflusplit_unsplit_bothendpts,
                                                 point_location='BOTH_ENDS')

        #Intersect end points
        net_conflusplit_unsplit_bothendpts_selfinters = os.path.join(
            out_gdb, '{0}_noartif_conflusplit_unsplit_bothendpts_selfinters_{1}'.format(prefix, suffix))
        arcpy.Intersect_analysis(in_features=[net_conflusplit_unsplit_bothendpts,
                                              net_conflusplit_unsplit_bothendpts],
                                 out_feature_class=net_conflusplit_unsplit_bothendpts_selfinters,
                                 join_attributes='ALL')

        #Intersect end points with lines
        net_conflusplit_unsplit_segs_bothendpts_inters = os.path.join(
            out_gdb, '{0}_noartif_conflusplit_unsplit_segs_bothendpts_inters_{1}'.format(prefix, suffix))
        arcpy.Intersect_analysis(in_features=[net_conflusplit_unsplit,
                                              net_conflusplit_unsplit_bothendpts],
                                 out_feature_class=net_conflusplit_unsplit_segs_bothendpts_inters,
                                 join_attributes='ALL')
        #If two end points intersect, then the lines form a loop
        bothendpts_selfinters_dict = defaultdict(list)
        with arcpy.da.SearchCursor(net_conflusplit_unsplit_bothendpts_selfinters,
                                   ['CONCATENATE_OBJECTID', 'CONCATENATE_OBJECTID_1']) as cursor:
            for row in cursor:
                if row[0] != row[1]:
                    bothendpts_selfinters_dict[row[0]].append(row[1])
        loop_ids = {}
        for k_1,v_1 in bothendpts_selfinters_dict.items():
            dupli_value = [k_2 for k_2, v_2 in collections.Counter(v_1).items() if v_2 == 2]
            #if len(dupli_value) > 1: #One occurrence of double-loop, ignore for now
            #   print(dupli_value)
            if len(dupli_value) == 1 and (dupli_value[0] not in loop_ids): #(dupli_value not in loop_ids) ensures that loops are included only once rather than once for each participating line
                loop_ids[k_1] = dupli_value[0]

        pd.DataFrame.from_dict(loop_ids, orient='index').\
            reset_index().\
            rename({'index':'OBJECTID', 0:'OBJECTID_1'}, axis='columns').\
            to_csv(loop_ids_tab_path, index='False')

    #Flip lines and assign Strahler order-------------------------------------------------------------------
    if not arcpy.Exists(net_conflusplit_directed):
        print("Flip lines...")
        arcpy.CopyFeatures_management(net_conflusplit, net_conflusplit_directed_ini)
        arcpy.MakeFeatureLayer_management(net_conflusplit_directed_ini, 'segs_to_flip_lyr',
                                          where_clause="elv_coef > 0.1")
        arcpy.edit.FlipLine(in_features='segs_to_flip_lyr')

        print("Assign first order to lines with dangle point and start point coinciding")
        if not arcpy.Exists(net_conflusplit_dangle):
            arcpy.FeatureVerticesToPoints_management(in_features=net_conflusplit_directed_ini,
                                                     out_feature_class=net_conflusplit_dangle,
                                                     point_location='DANGLE')

        if not arcpy.Exists(net_conflusplit_start):
            arcpy.FeatureVerticesToPoints_management(in_features=net_conflusplit_directed_ini,
                                                     out_feature_class=net_conflusplit_start,
                                                     point_location='START')

        if not arcpy.Exists(net_conflusplit_danglestart_inters):
            arcpy.Intersect_analysis(in_features=[net_conflusplit_dangle, net_conflusplit_start],
                                     out_feature_class=net_conflusplit_danglestart_inters,
                                     join_attributes='ALL')

        first_order_list = [row[0] for row in arcpy.da.SearchCursor(net_conflusplit_danglestart_inters, ['ORIG_FID'])]
        if 'strahler' not in [f.name for f in arcpy.ListFields(net_conflusplit_directed_ini)]:
            arcpy.AddField_management(net_conflusplit_directed_ini, 'strahler', 'SHORT')
        with arcpy.da.UpdateCursor(net_conflusplit_directed_ini, ['OID@', 'strahler']) as cursor:
            for row in cursor:
                if row[0] in first_order_list:
                    row[1] = 1
                    cursor.updateRow(row)

        arcpy.Delete_management(net_conflusplit_dangle)
        arcpy.Delete_management(net_conflusplit_start)

        #Some segments had not been well dissolved, so only part of the reach is order 1 and the rest is undefined
        #So: if a segment intersects only one strahler order 1, and does so from a start or end point, make it strahler order 1
        #Iterate multiple times until there are not such cases because there are sometimes multiple segment pieces one after the other
        arcpy.CopyFeatures_management(net_conflusplit_directed_ini, net_conflusplit_directed)

        for s_o in range(1,4):
            assign_strahler_splitline(in_lines=net_conflusplit_directed,
                                      out_gdb=out_gdb,
                                      prefix=prefix,
                                      suffix=suffix,
                                      in_strahler=s_o)

            #TO BE CONTINUED
            #For some reason, sometimes does not always process everything
            assign_strahler_splitline(in_lines=net_conflusplit_directed,
                                      out_gdb=out_gdb,
                                      prefix=prefix,
                                      suffix=suffix,
                                      in_strahler=s_o)

            iterate_strahler(in_lines = net_conflusplit_directed,
                             out_gdb = out_gdb,
                             prefix = prefix,
                             suffix = suffix,
                             in_loop_ids_tab_path = loop_ids_tab_path,
                             strahler_ini=s_o)

        #Re-run full loop from strahler order two to max.

        #Extend order to undefined one after loop

    #Remove lines with dangle points under 10 m
    # arcpy.MakeFeatureLayer_management(in_net, 'ddt_confusplit_sublyr',
    #                                   where_clause="orig_layer='D73_Savoie'")
    if not arcpy.Exists(out_net):
        print("Re-assign Strahler order attribute to original network...")
        arcpy.SpatialJoin_analysis(target_features=in_net,
                                   join_features=net_conflusplit_directed,
                                   out_feature_class=out_net,
                                   join_operation='JOIN_ONE_TO_ONE',
                                   join_type='KEEP_ALL',
                                   match_option='LARGEST_OVERLAP'
                                   )

#########################################################################################################################
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
    print('Join ddt net to hydrographic basins')
    arcpy.SpatialJoin_analysis(ddt_net_noartif,
                               join_features=BH,
                               out_feature_class=ddt_net_noartif_bh,
                               join_operation='JOIN_ONE_TO_ONE',
                               join_type='KEEP_ALL',
                               match_option='LARGEST_OVERLAP')

#Merge all lines between confluences and assign strahler order
bh_numset = {row[0] for row in arcpy.da.SearchCursor(ddt_net_noartif_bh, 'CdBH')}
for bh_num in bh_numset:
    if bh_num is not None:
        print(bh_num)
        arcpy.MakeFeatureLayer_management(ddt_net_noartif_bh, 'ddt_net_noartif_sub',
                                          where_clause="CdBH='{}'".format(bh_num))
        ddt_net_noartif_wstrahler = os.path.join(pregdb, 'carto_loi_eau_noartif_strahler_{}'.format(bh_num))
        identify_strahler(in_net='ddt_net_noartif_sub',
                          in_dem=bdalti_mosaic,
                          prefix='carto_loi_eau_fr',
                          suffix=bh_num,
                          out_gdb=pregdb,
                          out_net=ddt_net_noartif_wstrahler)

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

# Merge all lines between confluences and assign strahler order
bh_numset = {row[0] for row in arcpy.da.SearchCursor(bdtopo_noartif_bh, 'CdBH')}
for bh_num in bh_numset:
    if bh_num is not None:
        print(bh_num)
        arcpy.MakeFeatureLayer_management(bdtopo_noartif_bh, 'bdtopo_noartif_sub',
                                          where_clause="CdBH='{}'".format(bh_num))
        bdtopo_noartif_wstrahler = os.path.join(pregdb, 'bdtopo_noartif_strahler_{}'.format(bh_num))
        identify_strahler(in_net='bdtopo_noartif_sub',
                           in_dem=bdalti_mosaic,
                           suffix=bh_num,
                           out_gdb=pregdb,
                           out_net=bdtopo_noartif_wstrahler,
                           prefix='bdtopo')
#Check how many have a dangle point that is an end point - random sample


################################ Merge and export #####################################################################
if not arcpy.Exists(ddt_net_wstrahler):
    print('Merging all regions of ddt net with Strahler order')
    ddt_net_wstrahler_list = getfilelist(dir=pregdb, repattern='carto_loi_eau_noartif_strahler_[0-9]{2}$',
                                         nongdbf=False, gdbf=True, fullpath=True)
    arcpy.Merge_management(ddt_net_wstrahler_list, ddt_net_wstrahler)
    ftodelete_list = [f.name for f in arcpy.ListFields(ddt_net_wstrahler) if f.name not in
                      [arcpy.Describe(ddt_net_wstrahler).OIDFieldName, 'Shape', 'UID_CE', 'strahler', 'geom_Length',
                       'Shape_Length']]
    arcpy.DeleteField_management(ddt_net_wstrahler, drop_field=ftodelete_list)

if not arcpy.Exists(ddt_net_wstrahler_tab):
    print('Exporting attribute table to csv')
    arcpy.CopyRows_management(ddt_net_wstrahler, ddt_net_wstrahler_tab)

if not arcpy.Exists(bdtopo_wstrahler):
    print('Merging all regions of bd topo with Strahler order')
    bdtopo_wstrahler_list = getfilelist(dir=pregdb, repattern='bdtopo_noartif_strahler_[0-9]{2}$',
                                        nongdbf=False, gdbf=True, fullpath=True)
    arcpy.Merge_management(bdtopo_wstrahler_list, bdtopo_wstrahler)
    ftodelete_list = [f.name for f in arcpy.ListFields(bdtopo_wstrahler) if f.name not in
                      [arcpy.Describe(bdtopo_wstrahler).OIDFieldName, 'Shape', 'ID', 'strahler', 'Shape_Length']]
    arcpy.DeleteField_management(bdtopo_wstrahler, drop_field=ftodelete_list)
if not arcpy.Exists(bdtopo_wstrahler_tab):
    print('Exporting attribute table to csv')
    arcpy.CopyRows_management(bdtopo_wstrahler, bdtopo_wstrahler_tab)
