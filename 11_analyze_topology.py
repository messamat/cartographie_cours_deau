import os
import arcpy
import collections
import time

from setup_classement import *
from sklearn.linear_model import HuberRegressor, LinearRegression

overwrite = False

anci_dir = os.path.join(datdir, "données_auxiliaires")  # Ancillary data directory

# Preprocessing geodatabase
pregdb = os.path.join(resdir, "preprocessing_ancillary_data.gdb")

# Scratch gdb
temp_gdb = os.path.join(resdir, "scratch.gdb")
if not arcpy.Exists(temp_gdb):
    arcpy.CreateFileGDB_management(out_folder_path=os.path.split(temp_gdb)[0],
                                   out_name=os.path.split(temp_gdb)[1])

# DDT network
outputs_gdb = os.path.join(resdir, 'analysis_outputs.gdb')
ddt_net = os.path.join(outputs_gdb, 'carto_loi_eau_fr')

# ddt net after linkage to BD topo and BD carthage for auxilliary information (intermittence status, name, etc.)
ddt_net_refnetsattri_tab = os.path.join(resdir, 'carto_loi_eau_refnetsattris.csv')

# BDTOPO
bdtopo2015_fr = os.path.join(pregdb, 'bdtopo2015_fr')

# Large hydrographic units
BH = os.path.join(anci_dir, 'bdtopage', 'BassinHydrographique_FXX.shp')

# BD alti
bdalti_mosaic = os.path.join(pregdb, "bdalti_25m_mosaic")
# ------------------- OUTPUTS ------------------------------------------------------------------------------------------
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

# ------------- TRY TO ENHANCE DDT AND BDTOPO NETWORKS WITH FLOW DIRECTION --------------------------------------------------------
# Remove lines that are fully contained within other ones without deleting both lines if perfect overlap between two lines
def remove_contained_lines(in_net, idfn, precision_digits=2, delete_lenfield=True, check_type=True, temp_gdb='memory'):
    selfinters_lines = os.path.join(temp_gdb, 'selfinters_lines')

    # Self intersect lines
    if ('ORIG_LENGTH' in [f.name for f in arcpy.ListFields(in_net)]) and delete_lenfield:
        arcpy.DeleteField_management(in_net, drop_field='ORIG_LENGTH')

    arcpy.AddGeometryAttributes_management(in_net,
                                           Geometry_Properties="LENGTH_GEODESIC",
                                           Length_Unit='meters')
    arcpy.AlterField_management(in_net,
                                field='LENGTH_GEO',
                                new_field_name='ORIG_LENGTH',
                                new_field_alias='ORIG_LENGTH')

    arcpy.Intersect_analysis(in_features=[in_net, in_net],
                             out_feature_class=selfinters_lines,
                             join_attributes='ALL',
                             output_type='LINE')
    arcpy.AddGeometryAttributes_management(selfinters_lines, "LENGTH_GEODESIC", Length_Unit='meters')

    ids_not_to_delete = set()
    ids_to_delete = set()
    # Let's say that idfn is the id of the source line, and idfn_1 the id of the target line
    cursor_cols = [idfn, '{}_1'.format(idfn),'ORIG_LENGTH', 'LENGTH_GEO', 'ORIG_LENGTH_1']
    if check_type:
        cursor_cols.extend(['type_stand', 'type_stand_1'])

    with arcpy.da.UpdateCursor(selfinters_lines, cursor_cols) as cursor:
        for row in cursor:
            # Only analyze intersections between two different lines
            if row[0] != row[1]:
                # Continue if the source line was flagged to not delete
                if row[0] in ids_not_to_delete:
                    continue
                # if the source line is already flagged to delete, then add target line to those not to delete
                elif row[0] in ids_to_delete:
                    ids_not_to_delete.add(row[1])
                # if the length of the intersection of the source line with the target line is the same length as the original source line
                elif (round(row[2], precision_digits) == round(row[3], precision_digits)):
                    #If not checking whether there is a difference in line status, delete source line
                    if not check_type:
                        ids_to_delete.add(row[0])
                        ids_not_to_delete.add(row[1])
                    else:
                        # if the target line is also fully overlapping, and it does not have watercourses status while the source line does
                        #  delete the target line
                        if (round(row[4], precision_digits) == round(row[3], precision_digits)) \
                            and (row[5] in ["Cours d'eau", "Non cours d'eau"]) \
                            and (row[6] not in ["Cours d'eau", "Non cours d'eau"]):
                            ids_to_delete.add(row[1])
                            ids_not_to_delete.add(row[0])
                        else: #otherwise, delete the source line
                            ids_to_delete.add(row[0])
                            ids_not_to_delete.add(row[1])

    with arcpy.da.UpdateCursor(in_net, [idfn]) as cursor:
        for row in cursor:
            if row[0] in ids_to_delete:
                cursor.deleteRow()

    #arcpy.DeleteField_management(in_net, 'ORIG_LENGTH')
    arcpy.Delete_management(selfinters_lines)

    return [ids_to_delete, ids_not_to_delete]

def get_lr_coef(df, x, y):
    return LinearRegression().fit(df[[x]], df[[y]]).coef_[0][0]

def assign_strahler_splitline(in_net, prefix, suffix, idfield, in_strahler, temp_gdb='memory'):
    print("-------- Assigning same Strahler order to lines connected by one end "
          "exclusively to another line of order {}".format(in_strahler))

    for vertex_direction in ['START', 'END']:
        qualifying_endpts_selfinters_ids = {1}
        while len(qualifying_endpts_selfinters_ids) > 0:
            in_net_oneendpts = os.path.join(
                temp_gdb, '{0}_noartif_conflusplit_directed_{1}_{2}'.format(prefix, vertex_direction, suffix))
            arcpy.MakeFeatureLayer_management(in_net,  'in_net_null','strahler IS NULL')
            arcpy.FeatureVerticesToPoints_management(in_features='in_net_null',
                                                     out_feature_class=in_net_oneendpts,
                                                     point_location=vertex_direction)

            in_net_bothendpts = os.path.join(
                temp_gdb, '{0}_noartif_conflusplit_directed_bothendpts_{1}'.format(prefix, suffix))
            arcpy.FeatureVerticesToPoints_management(in_features=in_net,
                                                     out_feature_class=in_net_bothendpts,
                                                     point_location='BOTH_ENDS')

            in_net_bothendpts_selfinters = os.path.join(
                temp_gdb, '{0}_noartif_conflusplit_directed_bothendpts_selfinters_{1}'.format(prefix, suffix))
            arcpy.Intersect_analysis(in_features=[in_net_oneendpts,
                                                  in_net_bothendpts],
                                     out_feature_class=in_net_bothendpts_selfinters,
                                     join_attributes='ALL')

            #Identify links to undefined lines
            split_order_pts_selfinters = defaultdict(list)
            split_order_pts_selfinters_undefined = defaultdict(set)
            split_higherorder_pts_selfinters = set()
            with arcpy.da.SearchCursor(in_net_bothendpts_selfinters,
                                        [idfield, '{}_1'.format(idfield), 'strahler', 'strahler_1']) as cursor:
                for row in cursor:
                    #Identify undefined lines connected to another line of the order of interest:
                    if ((row[0] != row[1]) and (row[1] not in split_order_pts_selfinters[row[0]])
                            and (row[0] not in split_higherorder_pts_selfinters)
                            and (row[2] is None) and (row[3] is not None)):
                        if row[3] == in_strahler:
                            # If intersecting point is from different line, and other line is not already in dict of defined lines
                            # and strahler order is undefined for this line and in_strahler for this other line
                            split_order_pts_selfinters[row[0]].append(row[1])
                        #If other line is of higher order, exclude undefined line from further analysis
                        elif row[3] > in_strahler:
                            split_higherorder_pts_selfinters.add(row[0])

                    #Identify undefined lines connected to another undefined line
                    if ((row[0] != row[1]) and (row[1] not in split_order_pts_selfinters[row[0]])
                            and (row[2] is None) and (row[3] is None)):
                        # If intersecting point is from different line, and other line is not already in dict of undefined lines
                        # and strahler order is undefined for this line and in_strahler for this other line
                        split_order_pts_selfinters_undefined[row[0]].add(row[1])

            # Identify those "other defined lines" which intersect with only one undefined line
            # (e.g., to avoid first order streams intersecting already-second order streams)
            split_order_pts_selfinters_rev = defaultdict(list)
            for k, v in split_order_pts_selfinters.items():
                if k not in split_higherorder_pts_selfinters:
                    for i in v:
                        split_order_pts_selfinters_rev[i].append(k)
            endpts_selfinters_unique_match_ids = {k for k, v in split_order_pts_selfinters_rev.items() if len(v) == 1}

            qualifying_endpts_selfinters_ids = {k for k, v in split_order_pts_selfinters.items()
                                                if (len(v) == 1) and  # Make sure there aren't two first-order streams intersecting the main line
                                                (len(split_order_pts_selfinters_undefined[k]) == 0) and  # and that the undefined line does not intersect with another undefined line on that end
                                                v[0] in endpts_selfinters_unique_match_ids}

            with arcpy.da.UpdateCursor(in_net, [idfield, 'strahler']) as cursor:
                for row in cursor:
                    if row[0] in qualifying_endpts_selfinters_ids:
                        row[1] = in_strahler
                        cursor.updateRow(row)

    if temp_gdb == 'memory':
        arcpy.Delete_management('memory')

# Deal with simple loops (two lines sharing start and end points) ---------------------------
def inner_identify_loops(in_net, ref_fn, suffix, temp_gdb='memory'):
    net_conflusplit_unsplit_bothendpts = os.path.join(
        temp_gdb, 'identify_loops_bothendpts_{0}'.format(suffix))
    arcpy.FeatureVerticesToPoints_management(in_features=in_net,
                                             out_feature_class=net_conflusplit_unsplit_bothendpts,
                                             point_location='BOTH_ENDS')

    # Intersect end points
    net_conflusplit_unsplit_bothendpts_selfinters = os.path.join(
        temp_gdb, 'identify_loops_bothendpts_selfinters_{0}'.format(suffix))
    arcpy.Intersect_analysis(in_features=[net_conflusplit_unsplit_bothendpts,
                                          net_conflusplit_unsplit_bothendpts],
                             out_feature_class=net_conflusplit_unsplit_bothendpts_selfinters,
                             join_attributes='ALL')

    # If their two end points intersect, then the lines form a loop
    bothendpts_selfinters_dict = defaultdict(list)
    with arcpy.da.SearchCursor(net_conflusplit_unsplit_bothendpts_selfinters,
                               [ref_fn,'{}_1'.format(ref_fn)]) as cursor:
        for row in cursor:
            bothendpts_selfinters_dict[row[0]].append(row[1])

    loop_ids = {}
    #First include self loops
    for k_1, v_1 in bothendpts_selfinters_dict.items():
        #If line loops unto itself
        self_loop = [k_2 for k_2, v_2 in collections.Counter(v_1).items() if ((v_2 == 4) and (k_1 == k_2))]
        if len(self_loop) > 0 and (self_loop[0] not in loop_ids):
            loop_ids[self_loop[0]] = None

    #Then include multi-line loops
    for k_1, v_1 in bothendpts_selfinters_dict.items():
        if k_1 not in loop_ids:
            dupli_value = [k_2 for k_2, v_2 in collections.Counter(v_1).items() if ((v_2 == 2) and (k_1 != k_2))]
            # if len(dupli_value) > 1: #One occurrence of double-loop, ignore for now
            #   print(dupli_value)
            if len(dupli_value) == 1 and (dupli_value[0] not in loop_ids):
                # len(dupli_value) == 1 means that we focus on simple loops
                # (dupli_value not in loop_ids) ensures that loops are included only once rather than once for each participating line
                loop_ids[k_1] = dupli_value[0]


    loop_ids_pd = pd.DataFrame.from_dict(loop_ids, orient='index'). \
        reset_index(). \
        rename({'index': ref_fn, 0: '{}_1'.format(ref_fn)}, axis='columns')

    if temp_gdb == 'memory':
        arcpy.Delete_management('memory')

    return loop_ids_pd

# Identify loops" subset network to remove lines with strahler order, identify simple loops, unsplit lines
# re-identify loops, output panda data frame
def identify_loops(in_net, idfield, subset_net, prefix, suffix, temp_gdb='memory', verbose=False):
    print("-------- Identifying simple loops...")

    if subset_net:
        if verbose:
            print("------------ Subsetting network")
        subset_net = arcpy.MakeFeatureLayer_management(in_net, out_layer='net_nostrahler',
                                                       where_clause='strahler IS NULL')
        in_net = subset_net

    # Identify first set of loops before unsplitting
    if verbose:
        print("------------ Identifying loops before unsplitting")
    #oidfn = arcpy.Describe(in_net).OIDFieldName
    if not 'CONCATENATE_{}'.format(idfield) in [f.name for f in arcpy.ListFields(in_net)]:
        arcpy.AddField_management(in_net, 'CONCATENATE_{}'.format(idfield), 'LONG')
        with arcpy.da.UpdateCursor(in_net, ['CONCATENATE_{}'.format(idfield), idfield]) as cursor:
            for row in cursor:
                row[0] = row[1]
                cursor.updateRow(row)
            del row

    loop_ids_beforeunsplit = inner_identify_loops(in_net=in_net,
                                                  ref_fn='CONCATENATE_{}'.format(idfield),
                                                  temp_gdb=temp_gdb,
                                                  suffix=suffix)

    # Unsplit lines, then re-identify loops
    if verbose:
        print("------------ Identifying loops after unsplitting")
    net_conflusplit_unsplit = os.path.join(temp_gdb,
                                           '{0}_noartif_conflusplit_unsplit_{1}'.format(prefix, suffix))
    arcpy.UnsplitLine_management(in_features=in_net,
                                 out_feature_class=net_conflusplit_unsplit,
                                 statistics_fields=[[idfield, 'CONCATENATE']],
                                 concatenation_separator='_')
    loop_ids_afterunsplit = inner_identify_loops(in_net=net_conflusplit_unsplit,
                                                 ref_fn='CONCATENATE_{}'.format(idfield),
                                                 temp_gdb=temp_gdb,
                                                 suffix=suffix
                                                 )

    #Concatenate the loops before and after unsplitting.
    loop_ids = pd.concat([
        loop_ids_beforeunsplit.astype("Int64").astype(str).replace('<NA>',None), #Make sure that the unsplit IDs are string with None passed through
        loop_ids_afterunsplit],
                         axis=0)

    #Drop full duplicates
    loop_ids.drop_duplicates(inplace=True)

    # reset index to have a unique loop ID
    loop_ids.reset_index(drop=True, inplace=True)

    # IGNORE - #Flatten partial duplicates forming large loops if remaining duplicates
    # loop_ids_na = loop_ids[loop_ids['CONCATENATE_{}_1'.format(idfield)].isnull()]
    # loop_ids_nona = loop_ids[~(loop_ids['CONCATENATE_{}_1'.format(idfield)].isnull())]
    #
    # loop_ids_nona_dupli_1 = loop_ids_nona[loop_ids_nona.duplicated(subset='CONCATENATE_{}_1'.format(idfield), keep=False)]
    # loop_ids_nona_dupli = loop_ids_nona[loop_ids_nona.duplicated(subset='CONCATENATE_{}'.format(idfield), keep=False)]
    #
    # if len(loop_ids_nona_dupli_1)==0 and len(loop_ids_nona_dupli)==0:
    #    #continue process
    # else:
    #     raise Exception('there are remaining duplicates in the loops. the function needs to be edited to deal with it')
    #
    # (loop_ids_duplis.groupby('CONCATENATE_{}_1'.format(idfield)).agg(
    #     {'CONCATENATE_{}'.format(idfield):'_'.join}).
    #  reset_index())

    # Put loop ids in a key column, and OBJECTID of all constitutive lines in another column
    loop_ids_df = pd.concat([
        (loop_ids['CONCATENATE_{}'.format(idfield)].str.split('_', expand=True)),
        # Make sure that the unsplit IDs are string with None passed through
        (loop_ids['CONCATENATE_{}_1'.format(idfield)].str.split('_', expand=True))],
        axis=1). \
        reset_index(). \
        rename(columns={'index':'loop_id'}). \
        melt(id_vars='loop_id')

    loop_ids_df_format = loop_ids_df.loc[~loop_ids_df['value'].isnull(), ['loop_id', 'value']]
    loop_ids_df_format.rename(mapper={'value': idfield}, axis=1, inplace=True)

    #Remove loops with identical constitutive lines (due to self-looping unsplit lines having been added)
    loop_ids_nodupli = set()
    loop_values = set()
    for k,v in loop_ids_df_format.groupby('loop_id')[idfield].apply(frozenset).to_dict().items():
        if v not in loop_values:
            loop_ids_nodupli.add(k)
            loop_values.add(v)

    loop_ids_df_format_nodupli =  loop_ids_df_format[loop_ids_df_format['loop_id'].isin(loop_ids_nodupli)]

    if temp_gdb == 'memory':
        arcpy.Delete_management('memory')

    return loop_ids_df_format_nodupli

def assign_strahler_toloops(in_net, idfield, in_loop_ids_df, prefix, suffix, temp_gdb='memory'):
    print('-------- Assigning Strahler order to simple loops with only one connected undefined line - assign order to that line...')

    # For loops connected to a single undefined line, and other lines with defined strahler order
    # Identify the highest strahler order connected to it
    # If only one line of that order, assign that order to all lines in loop and the connected undefined line
    # If two lines of that order, assign that order +1 to all lines in loop and the connected undefined line
    in_net_bothends = os.path.join(temp_gdb, '{0}_noartif_conflusplit_bothendpts_{1}'.format(prefix, suffix))

    arcpy.FeatureVerticesToPoints_management(in_features=in_net,
                                             out_feature_class=in_net_bothends,
                                             point_location='BOTH_ENDS')

    net_conflusplit_unsplit_segs_bothendpts_inters = os.path.join(
        temp_gdb, '{0}_noartif_conflusplit_segs_bothendpts_inters_{1}'.format(prefix, suffix))
    arcpy.Intersect_analysis(in_features=[in_net,
                                          in_net_bothends],
                             out_feature_class=net_conflusplit_unsplit_segs_bothendpts_inters,
                             join_attributes='ALL')

    # Create a nested dictionary for which each outer key is a unique loop_id
    # for each loop_id, the outer value is another dictionary which contains the list of connected lines IDs (inner values) by strahler order (inner keys)
    loop_connections_dict = defaultdict(lambda: defaultdict(set))
    uids_inloops = np.array(in_loop_ids_df[idfield])
    with arcpy.da.SearchCursor(net_conflusplit_unsplit_segs_bothendpts_inters,
                               [idfield, '{}_1'.format(idfield), 'strahler_1']) as cursor:
        for row in cursor:
            if row[0] in uids_inloops:
                if row[1] not in uids_inloops:  # for now, only deal with simple loops (not series of loops)
                    sel_loopid = in_loop_ids_df.loc[in_loop_ids_df[idfield] == row[0], 'loop_id'].values[0]
                    loop_connections_dict[sel_loopid][row[2]].add(row[1])

    loop_connections_df = pd.DataFrame.from_dict(loop_connections_dict, orient='index'). \
        reset_index(). \
        rename(columns={'index': 'loop_id'})

    # Identify loops that intersect with a single undefined line or only defined lines
    loop_connections_1NA = loop_connections_df[loop_connections_df[np.NaN].str.len() == 1]. \
        melt(id_vars='loop_id')
    loop_connections_1NA['variable'] = loop_connections_1NA['variable'].astype('Int64')

    # Identify the maximum strahler order of the connected lines with a defined order
    max_connected_strahler = loop_connections_1NA.loc[~loop_connections_1NA['value'].isnull()]. \
        groupby('loop_id'). \
        max('variable')

    # Join
    loop_connections_1NA_maxstrahler = pd.merge(
        loop_connections_1NA,
        max_connected_strahler[~max_connected_strahler['variable'].isna()],
        how='inner', on=['loop_id', 'variable'])

    # Count and assign new order
    loop_connections_1NA_maxstrahler['new_strahler'] = loop_connections_1NA_maxstrahler['variable']
    loop_connections_1NA_maxstrahler.loc[
        loop_connections_1NA_maxstrahler['value'].str.len() > 1, 'new_strahler'] += 1

    # Prepare dictionaries with key as the line tempID and value as the new strahler order
    new_strahler_dict = pd.merge(in_loop_ids_df,
                                 loop_connections_1NA_maxstrahler[['loop_id', 'new_strahler']],
                                 how='inner', on='loop_id')[[idfield, 'new_strahler']]. \
        set_index(idfield)['new_strahler'].to_dict()  # Convert to dict

    new_strahler_connectedNAs = loop_connections_1NA.loc[loop_connections_1NA['variable'].isna()]
    new_strahler_connectedNAs.loc[:, 'value'] = new_strahler_connectedNAs.value.map(lambda x: list(x)[0])
    new_strahler_connectedNAs = pd.merge(new_strahler_connectedNAs,
                                         loop_connections_1NA_maxstrahler[['loop_id', 'new_strahler']],
                                         how='inner', on='loop_id')[['value', 'new_strahler']]. \
        set_index('value')['new_strahler'].to_dict()

    new_strahler_dict = new_strahler_dict | new_strahler_connectedNAs

    with arcpy.da.UpdateCursor(in_net, [idfield, 'strahler']) as cursor:
        for row in cursor:
            if row[0] in new_strahler_dict:
                row[1] = new_strahler_dict[row[0]]
                cursor.updateRow(row)

    if temp_gdb == 'memory':
        arcpy.Delete_management('memory')

    return new_strahler_dict

def iterate_strahler(in_net, idfield, in_loop_ids_df, prefix, suffix, strahler_ini, temp_gdb='memory'):
    print('-------- Assigning strahler order {} at confluences'.format(strahler_ini + 1))
    for vertex_direction in ['START', 'END']:
        in_net_oneendpts = os.path.join(
            temp_gdb, '{0}_noartif_conflusplit_directed_{1}_{2}'.format(prefix, vertex_direction, suffix))
        arcpy.FeatureVerticesToPoints_management(in_features=in_net,
                                                 out_feature_class=in_net_oneendpts,
                                                 point_location=vertex_direction)

        in_net_bothendpts = os.path.join(
            temp_gdb, '{0}_noartif_conflusplit_directed_bothendpts_{1}'.format(prefix, suffix))
        arcpy.FeatureVerticesToPoints_management(in_features=in_net,
                                                 out_feature_class=in_net_bothendpts,
                                                 point_location='BOTH_ENDS')

        in_net_bothendpts_selfinters = os.path.join(
            temp_gdb, '{0}_noartif_conflusplit_directed_bothendpts_selfinters_{1}'.format(prefix, suffix))
        arcpy.Intersect_analysis(in_features=[in_net_oneendpts,
                                              in_net_bothendpts],
                                 out_feature_class=in_net_bothendpts_selfinters,
                                 join_attributes='ALL')

        # Unique IDs of loops
        loop_ids_list = in_loop_ids_df[idfield].values

        # Get the set of lines intersecting each line with strahler_ini, excluding loops
        split_order_pts_selfinters = defaultdict(set)
        with arcpy.da.SearchCursor(in_net_bothendpts_selfinters,
                                    [idfield, '{}_1'.format(idfield), 'strahler', 'strahler_1']) as cursor:
            for row in cursor:
                if (((row[0] != row[1]) and (row[1] not in split_order_pts_selfinters[row[0]])
                    and (row[2] is None) and (row[3] == strahler_ini))
                        and (not row[0] in loop_ids_list)): #(not row[1] in loop_ids_list) not needed because dealing with loops include downstream segments
                    # If intersecting point is from different line, and other line is not already in dict
                    # and strahler order is undefined for this line and strahler_ini for this other line
                    # and this other line is not a flagged loop line
                    split_order_pts_selfinters[row[0]].add(row[1])

        # If intersecting at least two lines with strahler_ini, assign strahler_ini + 1
        with arcpy.da.UpdateCursor(in_net, [idfield, 'strahler']) as cursor:
            for row in cursor:
                if row[0] in split_order_pts_selfinters:
                    if len(split_order_pts_selfinters[row[0]]) > 1:
                        row[1] = strahler_ini + 1
                    cursor.updateRow(row)

    if temp_gdb == 'memory':
        arcpy.Delete_management('memory')

def run_strahler_routing_loop(in_net, idfield, prefix, suffix, so_min, so_max, temp_gdb='in_memory'):
    for so in range(so_min, so_max + 1):
        # loop_ids_temptab = os.path.join(resdir, '{0}_loop_ids_temptab_{1}.csv'.format(prefix, suffix))
        print('------ Expanding stream order {} downstream...'.format(so))
        assigned_loops = {9999: 999}
        while len(assigned_loops.keys()) > 0:
            # Expand strahler down until next confluence
            # Some segments had not been well dissolved, so only part of the reach is assigned the order and the rest is undefined
            # So: if a segment intersects only one segment of the given strahler order, and does so from a start or end point, make it of that order
            start = time.time()
            assign_strahler_splitline(in_net=in_net,
                                      idfield=idfield,
                                      temp_gdb=temp_gdb,
                                      prefix=prefix,
                                      suffix=suffix,
                                      in_strahler=so)
            print(time.time()-start)

            # Identify loops with undefined strahler order
            loop_ids_df = identify_loops(in_net=in_net,
                                         idfield=idfield,
                                         temp_gdb=temp_gdb,
                                         subset_net=True,
                                         prefix=prefix,
                                         suffix=suffix)
            loop_ids_df = loop_ids_df.astype(float).astype(
                int)  # First convert to float to avoid potential ValueError: invalid literal for int() with base 10: '4e+04'

            # If there are loops: assign them as Strahler order, and that of the downstream segment
            if len(loop_ids_df) > 0:
                assigned_loops = assign_strahler_toloops(in_net=in_net,
                                                         in_loop_ids_df=loop_ids_df,
                                                         idfield=idfield,
                                                         temp_gdb=temp_gdb,
                                                         prefix=prefix,
                                                         suffix=suffix)
                print('-------- Assigned Strahlher order to {0} segments associated'
                      ' with simple loops for stream order {1}'.format(
                    len(assigned_loops.keys()), so))
            else:
                assigned_loops = {}

        # Identify remaining loops with undefined strahler order so as not to use them in defining new Strahler order
        loop_ids_df = identify_loops(in_net=in_net,
                                     idfield=idfield,
                                     temp_gdb=temp_gdb,
                                     subset_net=True,
                                     prefix=prefix,
                                     suffix=suffix)
        loop_ids_df = loop_ids_df.astype(float).astype(int)

        iterate_strahler(in_net=in_net,
                         temp_gdb=temp_gdb,
                         in_loop_ids_df=loop_ids_df,
                         idfield=idfield,
                         prefix=prefix,
                         suffix=suffix,
                         strahler_ini=so)

    if temp_gdb == 'memory':
        arcpy.Delete_management('memory')

def assign_strahler_to_stubs(in_net, in_stublist, idfield, prefix, suffix, max_stub_length=500, temp_gdb='memory'):
    in_net_bothends = os.path.join(temp_gdb, '{0}_noartif_conflusplit_bothendpts_{1}'.format(prefix, suffix))

    arcpy.FeatureVerticesToPoints_management(in_features=in_net,
                                             out_feature_class=in_net_bothends,
                                             point_location='BOTH_ENDS')

    net_conflusplit_unsplit_segs_bothendpts_inters = os.path.join(
        temp_gdb, '{0}_noartif_conflusplit_segs_bothendpts_inters_{1}'.format(prefix, suffix))
    arcpy.Intersect_analysis(in_features=[in_net_bothends,
                                          in_net],
                             out_feature_class=net_conflusplit_unsplit_segs_bothendpts_inters,
                             join_attributes='ALL')

    # if id in stublist, length under 500 m, no strahler order, and connected line has strahler order.
    # When multiple stubs are connected to a line with a stream order -- only assign stream order to one stub
    stubs_so1_list = set()
    connected_to_stubs_list = set()
    with arcpy.da.SearchCursor(net_conflusplit_unsplit_segs_bothendpts_inters,
                                [idfield, 'LENGTH_GEO', 'strahler', 'strahler_1', "{}_1".format(idfield)]) as cursor:
        for row in cursor:
            if ((row[0] in in_stublist) and (row[1] <= max_stub_length)
                    and (row[2] is None) and (row[3] is not None) and
                    (row[4] not in connected_to_stubs_list)):
                stubs_so1_list.add(row[0])
                connected_to_stubs_list.add(row[4])

    # Assign strahler order 1
    with arcpy.da.UpdateCursor(in_net, [idfield, 'strahler']) as cursor:
        for row in cursor:
            if row[0] in stubs_so1_list:
                row[1] = 1
                cursor.updateRow(row)

    if temp_gdb == 'memory':
        arcpy.Delete_management('memory')

    return stubs_so1_list

def identify_strahler(in_net, in_dem, out_gdb, out_net, prefix, suffix, temp_gdb='memory', overwrite=True):
    tempIDfn = 'tempID'
    net_noduplisubsegs = os.path.join(temp_gdb, '{0}_noduplisubsegs_{1}'.format(prefix, suffix))
    net_dissolved = os.path.join(temp_gdb, '{0}_dissolved_{1}'.format(prefix, suffix))
    net_conflupts = os.path.join(temp_gdb, '{0}_conflupts_{1}'.format(prefix, suffix))
    net_conflusplit = os.path.join(temp_gdb, '{0}_conflusplit_{1}'.format(prefix, suffix))
    net_nodes25m = os.path.join(temp_gdb, '{0}_nodes25m_{1}'.format(prefix, suffix))
    net_nodes_elv = os.path.join(temp_gdb, '{0}_nodes25m_elv_{1}'.format(prefix, suffix))
    net_conflusplit_directed_ini = os.path.join(out_gdb,
                                                '{0}_conflusplit_directed_ini_{1}'.format(prefix, suffix))
    net_conflusplit_dangle = os.path.join(temp_gdb, '{0}_conflusplit_dangle_{1}'.format(prefix, suffix))
    net_conflusplit_start = os.path.join(temp_gdb, '{0}_conflusplit_start_{1}'.format(prefix, suffix))
    net_conflusplit_danglestart_inters = os.path.join(temp_gdb,
                                                      '{0}_conflusplit_danglestartinters_{1}'.format(prefix,
                                                                                                             suffix))
    net_conflusplit_directed = os.path.join(out_gdb, '{0}_conflusplit_directed_{1}'.format(prefix, suffix))

    if (not arcpy.Exists(net_noduplisubsegs)) or overwrite:
        print("---- Deleting overlapping segment sections")
        arcpy.Intersect_analysis(in_features=[in_net, in_net],
                                 out_feature_class=net_noduplisubsegs,
                                 join_attributes='ALL',
                                 output_type='LINE')

        arcpy.DeleteIdentical_management(net_noduplisubsegs,
                                         fields=arcpy.Describe(net_noduplisubsegs).shapeFieldName)

    if (not arcpy.Exists(net_dissolved)) or overwrite:
        print('---- Dissolving the  network')
        arcpy.Dissolve_management(net_noduplisubsegs,
                                  out_feature_class=net_dissolved,
                                  multi_part='SINGLE_PART',
                                  unsplit_lines='UNSPLIT_LINES')

    if (not arcpy.Exists(net_conflupts)) or overwrite:
        print('---- Getting confluence points to make sure that all lines are split at confluences')
        arcpy.Intersect_analysis(in_features=net_dissolved,
                                 out_feature_class=net_conflupts,
                                 join_attributes='ONLY_FID',
                                 output_type='POINT')

    if (not arcpy.Exists(net_conflusplit)) or overwrite:
        print('---- Splitting dissolved network at confluence points')
        arcpy.SplitLineAtPoint_management(in_features=net_dissolved,
                                          point_features=net_conflupts,
                                          out_feature_class=net_conflusplit,
                                          search_radius=0.1)

    print('---- Identifying first-order streams')
    if (not arcpy.Exists(net_nodes25m)) or overwrite:
        print('-------- Generating points along lines')
        arcpy.GeneratePointsAlongLines_management(Input_Features=net_conflusplit,
                                                  Output_Feature_Class=net_nodes25m,
                                                  Point_Placement="DISTANCE",
                                                  Distance="25 meters",
                                                  Include_End_Points="END_POINTS")

    if (not arcpy.Exists(net_nodes_elv)) or overwrite:
        print('-------- Getting elevation for each point')
        oidfn = arcpy.Describe(net_nodes25m).OIDFieldName
        Sample(in_rasters=in_dem,
               in_location_data=net_nodes25m,
               out_table=net_nodes_elv,
               resampling_type='NEAREST',
               unique_id_field=oidfn,
               statistics_type='MEDIAN',
               layout='ROW_WISE',
               generate_feature_class='TABLE')

        colstoread = [f.name for f in arcpy.ListFields(net_nodes_elv) if
                      f.type != "Geometry"]  # List the fields you want to include. I want all columns except the geometry
        # Read elevation table and join OBJECTID of line associated with the node
        net_nodes_elv_pd = pd.DataFrame(data=arcpy.da.SearchCursor(net_nodes_elv, colstoread), columns=colstoread). \
            merge(pd.DataFrame(data=arcpy.da.SearchCursor(net_nodes25m, ['OID@', 'ORIG_FID']),
                               columns=[os.path.basename(net_nodes25m), 'ORIG_FID']),
                  on=os.path.basename(net_nodes25m))

        print("-------- Computing robust linear regression of elevation for each line...")
        coefs_dict = (net_nodes_elv_pd.loc[~net_nodes_elv_pd.bdalti_25m_mosaic_Band_1.isna()].groupby('ORIG_FID').
                      apply(get_lr_coef, x=oidfn, y='bdalti_25m_mosaic_Band_1').
                      to_dict())

        elv_coef_fname = 'elv_coef'
        if not elv_coef_fname in [f.name for f in arcpy.ListFields(net_conflusplit)]:
            arcpy.AddField_management(net_conflusplit, elv_coef_fname, 'DOUBLE')
        with arcpy.da.UpdateCursor(net_conflusplit, ['OID@', elv_coef_fname]) as cursor:
            for row in cursor:
                if row[0] in coefs_dict:
                    row[1] = coefs_dict[row[0]]
                    cursor.updateRow(row)

    # Flip lines and assign first Strahler order to out-most segments --------------------------------------------------
    if (not arcpy.Exists(net_conflusplit_directed_ini)) or overwrite:
        arcpy.CopyFeatures_management(net_conflusplit, net_conflusplit_directed_ini)

        # Create temporary ID field
        if not tempIDfn in [f.name for f in arcpy.ListFields(net_conflusplit_directed_ini)]:
            arcpy.AddField_management(net_conflusplit_directed_ini, tempIDfn, 'LONG')
            arcpy.CalculateField_management(
                in_table=net_conflusplit_directed_ini,
                field=tempIDfn,
                expression="!{}!".format(arcpy.Describe(net_conflusplit_directed_ini).OIDFieldName)
            )

        print("-------- Flipping lines oriented against the slope")
        arcpy.MakeFeatureLayer_management(net_conflusplit_directed_ini, 'segs_to_flip_lyr',
                                          where_clause="elv_coef > 0.1")
        arcpy.edit.FlipLine(in_features='segs_to_flip_lyr')

        print("-------- Identifying lines with dangle points")
        arcpy.FeatureVerticesToPoints_management(in_features=net_conflusplit_directed_ini,
                                                 out_feature_class=net_conflusplit_dangle,
                                                 point_location='DANGLE')
        stub_list = {row[0] for row in arcpy.da.SearchCursor(net_conflusplit_dangle, [tempIDfn])}
        print("-------- Deleted {} stubs: lines with dangle points under 10 m".format(len(stub_list)))
        geofn = arcpy.Describe(net_conflusplit_directed_ini).shapeFieldName
        with arcpy.da.UpdateCursor(net_conflusplit_directed_ini, [tempIDfn, '{}_Length'.format(geofn)]) as cursor:
            for row in cursor:
                if (row[0] in stub_list) and (row[1] < 10):
                    cursor.deleteRow()

        arcpy.FeatureVerticesToPoints_management(in_features=net_conflusplit_directed_ini,
                                                 out_feature_class=net_conflusplit_start,
                                                 point_location='START')

        arcpy.Intersect_analysis(in_features=[net_conflusplit_dangle, net_conflusplit_start],
                                 out_feature_class=net_conflusplit_danglestart_inters,
                                 join_attributes='ALL')

        print("-------- Assigning first order to lines with dangle point and start point coinciding")
        first_order_list = [row[0] for row in arcpy.da.SearchCursor(net_conflusplit_danglestart_inters, [tempIDfn])]
        if 'strahler' not in [f.name for f in arcpy.ListFields(net_conflusplit_directed_ini)]:
            arcpy.AddField_management(net_conflusplit_directed_ini, 'strahler', 'SHORT')
        with arcpy.da.UpdateCursor(net_conflusplit_directed_ini, [tempIDfn, 'strahler']) as cursor:
            for row in cursor:
                if row[0] in first_order_list:
                    row[1] = 1
                    cursor.updateRow(row)

    # Define Straler order for the rest of the network -----------------------------------------------------------------
    if (not arcpy.Exists(net_conflusplit_directed)) or overwrite:
        arcpy.CopyFeatures_management(net_conflusplit_directed_ini, net_conflusplit_directed)

        #Create length field
        arcpy.AddGeometryAttributes_management(net_conflusplit_directed, Geometry_Properties='LENGTH_GEODESIC')

        print("---- Running Strahler routing loop")
        #Run first loop iterating through all stream orders and dealing with simple loops
        #Modified in_net in place
        run_strahler_routing_loop(in_net=net_conflusplit_directed,
                                  temp_gdb='memory',
                                  idfield=tempIDfn,
                                  prefix=prefix,
                                  suffix=suffix,
                                  so_min=1,
                                  so_max=5)

        # Assign order 1 to lines with dangle points
        print("---- Assigning Strahler to stubs")
        stubs_so1_assigned = assign_strahler_to_stubs(in_net=net_conflusplit_directed,
                                                      in_stublist=stub_list,
                                                      temp_gdb='memory',
                                                      idfield=tempIDfn,
                                                      max_stub_length=500,
                                                      prefix=prefix,
                                                      suffix=suffix)

        print("-------- Assigning Strahler order 1 to {} stubs connected to lines with defined order".format(
            len(stubs_so1_assigned)))
        for so in reversed(range(5)):
            assign_strahler_splitline(in_net=net_conflusplit_directed,
                                      idfield=tempIDfn,
                                      temp_gdb='memory',
                                      prefix=prefix,
                                      suffix=suffix,
                                      in_strahler=so+1)

        run_strahler_routing_loop(in_net=net_conflusplit_directed,
                                  temp_gdb='memory',
                                  idfield=tempIDfn,
                                  prefix=prefix,
                                  suffix=suffix,
                                  so_min=1,
                                  so_max=5)

    # Remove lines with dangle points under 10 m
    # arcpy.MakeFeatureLayer_management(in_net, 'ddt_confusplit_sublyr',
    #                                   where_clause="orig_layer='D73_Savoie'")
    if (not arcpy.Exists(out_net)) or overwrite:
        print("---- Re-assigning Strahler order attribute to original network...")
        arcpy.SpatialJoin_analysis(target_features=in_net,
                                   join_features=net_conflusplit_directed,
                                   out_feature_class=out_net,
                                   join_operation='JOIN_ONE_TO_ONE',
                                   join_type='KEEP_ALL',
                                   match_option='LARGEST_OVERLAP'
                                   )

def enhance_network_topology(in_net, idfn, in_dem, out_net, temp_gdb, prefix, suffix, overwrite=False):
    start_bh = time.time()
    net_integrate = os.path.join(temp_gdb, '{0}_integrate_{1}'.format(prefix, suffix))

    #Turn off metadata logging for performance
    if float(arcpy.GetInstallInfo()['Version']) >= 3:
        if arcpy.GetLogMetadata():
            arcpy.SetLogMetadata(False)

    # Remove all lines which are fully contained within another line
    total_deleted_lines = 0
    if (not arcpy.Exists(net_integrate)) or overwrite:
        arcpy.CopyFeatures_management(in_net, net_integrate)

        ids_not_to_delete = {9999}
        while len(ids_not_to_delete) > 0:
            removal_output = remove_contained_lines(in_net=net_integrate,
                                                    temp_gdb=temp_gdb,
                                                    idfn=idfn,
                                                    check_type=False)
            ids_not_to_delete = removal_output[1]
            total_deleted_lines += len(removal_output[0])

        # Merge vertices that are within 10 cm from each other
        # (to deal with extremely small disconnections and nearly overlapping lines)
        arcpy.Integrate_analysis(in_features=net_integrate,
                                 cluster_tolerance='0.1')

        # Integrating joined lines that were within a few centimeters from each other created more overlapping
        # features -- delete them.
        ids_not_to_delete = {9999}
        while len(ids_not_to_delete) > 0:
            removal_output = remove_contained_lines(in_net=net_integrate,
                                                    temp_gdb=temp_gdb,
                                                    idfn=idfn,
                                                    check_type=False)
            ids_not_to_delete = removal_output[1]
            total_deleted_lines += len(removal_output[0])

    print('---- Deleted {} lines which were fully overlapping others'.format(total_deleted_lines))

    # Assign strahler order
    identify_strahler(in_net=net_integrate,
                      in_dem=in_dem,
                      prefix=prefix,
                      suffix=suffix,
                      out_gdb=temp_gdb,
                      temp_gdb='memory',
                      out_net=out_net,
                      overwrite=True)

    end_bh = time.time()
    print("RUNNING STRAHLER ROUTINE TOOK {0} s FOR RECORD {1}".format(
        round(end_bh - start_bh), suffix))

############################ RUN ANALYSIS FOR DDT NETWORKS #############################################################
# Remove all truly artificial watercourses that may disrupt topology
if not arcpy.Exists(ddt_net_noartif):
    arcpy.CopyFeatures_management(ddt_net, ddt_net_artif)
    arcpy.CopyFeatures_management(ddt_net, ddt_net_noartif)

    ddt_net_refnetsattri_pd = pd.read_csv(ddt_net_refnetsattri_tab)
    regex = r'\b)|(\b'.join([r"(\bcanal", "foss[eé]", "roubine", "craste", "d[ée]rivation",
                             "bief", "aber", r"hydraulique\b)"])
    ddt_net_refnetsattri_pd['NOM_bdtopo_orig'] = ddt_net_refnetsattri_pd['NOM_bdtopo_orig'].str.lower()
    ddt_net_refnetsattri_pd['TOPONYME1_carthage'] = ddt_net_refnetsattri_pd['TOPONYME1_carthage'].str.lower()
    artif_ids = list(
        ddt_net_refnetsattri_pd.loc[(ddt_net_refnetsattri_pd['NOM_bdtopo_orig'].str.contains(regex) |
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

# Spatial join to large basin region
if not arcpy.Exists(ddt_net_noartif_bh):
    print('Join ddt net to hydrographic basins')
    arcpy.SpatialJoin_analysis(ddt_net_noartif,
                               join_features=BH,
                               out_feature_class=ddt_net_noartif_bh,
                               join_operation='JOIN_ONE_TO_ONE',
                               join_type='KEEP_ALL',
                               match_option='LARGEST_OVERLAP')

# Merge all lines between confluences and assign strahler order
bh_numset = {row[0] for row in arcpy.da.SearchCursor(ddt_net_noartif_bh, 'CdBH')}
for bh_num in ['04', '05', '06']:#bh_numset:
    if bh_num is not None:
        print("PROCESSING HYDROGAPHIC BASIN {}".format(bh_num))
        temp_gdb = os.path.join(resdir, "scratch_{}.gdb".format(bh_num))
        if not arcpy.Exists(temp_gdb):
            arcpy.CreateFileGDB_management(out_folder_path=os.path.split(temp_gdb)[0],
                                           out_name=os.path.split(temp_gdb)[1])

        ddt_net_noartif_sub = os.path.join(temp_gdb, 'ddt_net_noartif_sub_{}'.format(bh_num))
        if not arcpy.Exists(ddt_net_noartif_sub):
            arcpy.MakeFeatureLayer_management(ddt_net_noartif_bh, 'ddt_net_noartif_sub',
                                              where_clause="CdBH='{}'".format(bh_num))
            arcpy.CopyFeatures_management('ddt_net_noartif_sub', ddt_net_noartif_sub)
            arcpy.Delete_management('ddt_net_noartif_sub')

        net_wstrahler = os.path.join(pregdb, 'ddt_net_noartif_strahler_{}'.format(bh_num))
        enhance_network_topology(in_net=ddt_net_noartif_sub,
                                 idfn='UID_CE',
                                 in_dem=bdalti_mosaic,
                                 temp_gdb=temp_gdb,
                                 out_net=net_wstrahler,
                                 prefix='ddt_net_noartif',
                                 suffix=bh_num,
                                 overwrite=True)

############################ RUN ANALYSIS FOR BDTOPO ###################################################################
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
for bh_num in ['02', '03', '04']:#bh_numset:
    if bh_num is not None:
        print("PROCESSING HYDROGAPHIC BASIN {}".format(bh_num))
        temp_gdb = os.path.join(resdir, "scratch_{}.gdb".format(bh_num))
        if not arcpy.Exists(temp_gdb):
            arcpy.CreateFileGDB_management(out_folder_path=os.path.split(temp_gdb)[0],
                                           out_name=os.path.split(temp_gdb)[1])

        bdtopo_noartif_sub = os.path.join(temp_gdb, 'bdtopo_noartif_sub_{}'.format(bh_num))
        arcpy.MakeFeatureLayer_management(bdtopo_noartif_bh, 'bdtopo_noartif_sub',
                                          where_clause="CdBH='{}'".format(bh_num))
        arcpy.CopyFeatures_management('bdtopo_noartif_sub', bdtopo_noartif_sub)
        arcpy.Delete_management('bdtopo_noartif_sub')

        bdtopo_wstrahler = os.path.join(pregdb, 'bdtopo_noartif_strahler_{}'.format(bh_num))
        enhance_network_topology(in_net=bdtopo_noartif_sub,
                                 idfn='ID',
                                 in_dem=bdalti_mosaic,
                                 temp_gdb=temp_gdb,
                                 out_net=bdtopo_wstrahler,
                                 prefix='bdtopo_noartif',
                                 suffix=bh_num,
                                 overwrite=True)


# Check how many have a dangle point that is an end point - random sample

################################ Merge and export #####################################################################
######################REMERGE ARTIFICIAL LINES#########################################################################
if not arcpy.Exists(ddt_net_wstrahler):
    print('Merging all regions of ddt net with Strahler order')
    ddt_net_wstrahler_list = getfilelist(dir=pregdb, repattern='carto_loi_eau_noartif_strahler_[0-9]{2}$',
                                         nongdbf=False, gdbf=True, fullpath=True)
    arcpy.Merge_management(ddt_net_wstrahler_list, ddt_net_wstrahler)
    ftodelete_list = [f.name for f in arcpy.ListFields(ddt_net_wstrahler) if f.name not in
                      [arcpy.Describe(ddt_net_wstrahler).OIDFieldName, 'Shape', 'UID_CE', 'strahler',
                       arcpy.Describe(ddt_net_wstrahler).shapeFieldName]]
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
                      [arcpy.Describe(bdtopo_wstrahler).OIDFieldName, 'Shape', 'ID', 'strahler',
                       arcpy.Describe(bdtopo_wstrahler).shapeFieldName]]
    arcpy.DeleteField_management(bdtopo_wstrahler, drop_field=ftodelete_list)

if not arcpy.Exists(bdtopo_wstrahler_tab):
    print('Exporting attribute table to csv')
    arcpy.CopyRows_management(bdtopo_wstrahler, bdtopo_wstrahler_tab)

# ------------- IDENTIFY DISCONNECTIONS --------------------------------------------------------------------------------
# Convert lines to end points
if not arcpy.Exists(ddt_net_endpts):
    arcpy.FeatureVerticesToPoints_management(in_features=ddt_net,
                                             out_feature_class=ddt_net_endpts,
                                             point_location="BOTH_ENDS")
# Intersect end points
if not arcpy.Exists(ddt_net_endpts_inters):
    arcpy.analysis.Intersect(in_features=[ddt_net_endpts, ddt_net_endpts],
                             out_feature_class=ddt_net_endpts_inters,
                             join_attributes='ALL')

# Identify isolated river segmenets (i.e., segmenets surrounded by segments with another category. e.g., nce surrounded by ce)
if not arcpy.Exists(ddt_net_isolatedsegs):
    # Identify lines with one connected end point and one connected start point with dictionary
    lines_inters_dict = {}
    # [f.name for f in arcpy.ListFields(ddt_net_endpts_inters)]
    with arcpy.da.SearchCursor(ddt_net_endpts_inters,
                               ['OID@', 'FID_carto_loi_eau_fr_endpts', 'UID_CE', 'UID_CE_1',
                                'type_stand', 'type_stand_1']) as cursor:
        for row in cursor:
            if row[2] != row[3]:
                lines_inters_dict[row[0]] = row[1:]
    lines_inters_pd = pd.DataFrame.from_dict(lines_inters_dict, orient='index')
    del lines_inters_dict
    lines_inters_pd.columns = ['FID_carto_loi_eau_fr_endpts', 'UID_CE', 'UID_CE_1', 'type_stand', 'type_stand_1']
    # Count number of intersecting points per end point
    lines_inters_pd_pts_interscount = lines_inters_pd.groupby(by=['FID_carto_loi_eau_fr_endpts', 'UID_CE']).size()
    # Count number of end points with intersecting points per segment
    lines_inters_pd_ptswinters = lines_inters_pd_pts_interscount.groupby(by='UID_CE').size()
    # Examine only those with two end points intersecting with other points (excluding first and last order streams)
    # And segments with a total of intersecting points of 3 and under, including those either without confluence
    # with upstream confluence or downstream bifurcation
    sel1 = (lines_inters_pd_pts_interscount.groupby('UID_CE').sum() < 3)
    sel2 = (lines_inters_pd_ptswinters == 2)

    simple_lines_pd = lines_inters_pd.loc[
                      (lines_inters_pd.UID_CE.isin(list(sel1[sel1].index))) &
                      (lines_inters_pd.UID_CE.isin(list(sel2[sel2].index))),
                      :]

    # For those lines, select those whose type of both ends differ from their own
    simple_lines_pd['dif_i'] = np.where(simple_lines_pd.type_stand != simple_lines_pd.type_stand_1, 0, 1)
    sel3 = (simple_lines_pd.groupby('UID_CE').sum('dif_i').dif_i == 0)
    simple_lines_pd_isolated_ID = sel3[sel3].index

    # Copy ddt_net the delete all features but those to inspect, add fields with types of the surrounding segments
    arcpy.CopyFeatures_management(ddt_net, ddt_net_isolatedsegs)
    with arcpy.da.UpdateCursor(ddt_net_isolatedsegs, ["UID_CE"]) as cursor:
        for row in cursor:
            if row[0] not in simple_lines_pd_isolated_ID:
                cursor.deleteRow()

# Manually go through them to check which ones are "unjustified"
# If line segment only partly contained in standing water body
# 6071, 7628, 8646, 9834, 9916,


# Check percentage of NCE for which surrounding ones are CE
# Could expand analysis to all departments with BDTOPO based on ID
