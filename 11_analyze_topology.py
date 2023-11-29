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

#DDT network
outputs_gdb = os.path.join(resdir, 'analysis_outputs.gdb')
ddt_net = os.path.join(outputs_gdb, 'carto_loi_eau_fr')

# BD alti
bdalti_mosaic = os.path.join(pregdb, "bdalti_25m_mosaic")


#------------------- OUTPUTS ------------------------------------------------------------------------------------------
ddt_net_endpts = os.path.join(pregdb, 'carto_loi_eau_fr_endpts')
ddt_net_endpts_inters = os.path.join(pregdb, 'carto_loi_eau_fr_endpts_inters')
ddt_net_isolatedsegs = os.path.join(pregdb, 'carto_loi_eau_fr_isolated')

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

#Check percentage of NCE for which surrounding ones are CE

#------------- TRY TO ENHANCE DDT NETWORK WITH FLOW DIRECTION --------------------------------------------------------
#Start trying with BD Alti
#(Maybe smooth elevation)

#Remove all those for which BD Topo is ARTIF and

#Merge all lines between confluences
ddt_net_dissolved = os.path.join(pregdb, 'carto_loi_eau_fr_dissolved')
arcpy.Dissolve_management(ddt_net,
                          out_feature_class=ddt_net_dissolved,
                          multi_part='SINGLE_PART',
                          unsplit_lines='UNSPLIT_LINES')
#arcpy.Intersect_analysis()
#arcpy.SplitLineAtPoint_management()

#Extract elevation at very node along line
#Run a regression across nodes
#Establish direction of line based on regression

#Then identify all lines with a dangle point that is also a start point
#Check how many have a dangle point that is an end point - random sample

#Check adjacent lines:
#if end point is connected to an endpoint (or multiple endpoints), and start point is connected to a start point