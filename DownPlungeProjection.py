#script to project map data into a cross section
# Ralph Haugerud, U.S. Geological Survey
# Code not subject to U.S. copyright

import arcpy, sys, os.path, math
import numpy as np

versionString = 'DownPlungeProjection.py, version of 11 February 2020'

'''
gdb    # GeMS database
dem    # coextensive with area to be projected x,y,z units are same as those of GeologicMap 
numberOfCells  # sampling distance, as cellSize multiplier
xsLine # feature class, ONE arc selected. We assume arc is STRAIGHT
token  # 'CS'+token = prefix for new feature classes
       # 'CrossSection'+token = name of new feature dataset
xsPlunge # plunge down which we project data. = 0 if no plunge. In degrees
bufferWidth # half-width, in map units, area around xsLine to be projected
'''
debug = False
debug2 = False

gdb = sys.argv[1]
dem = sys.argv[2]
numberOfCells = float(sys.argv[3])
xsLine = sys.argv[4]
token = sys.argv[5]
xsPlunge = float(sys.argv[6])
bufferWidth = float(sys.argv[7])
ForceExit = sys.argv[8]

tanPlunge = math.tan(math.radians(xsPlunge))
cosPlunge = math.cos(math.radians(xsPlunge))
sampleDistance = numberOfCells * arcpy.Describe(dem).meanCellWidth

#################################
def addMsgAndPrint(msg, severity=0): 
    # prints msg to screen and adds msg to the geoprocessor (in case this is run as a tool) 
    #print msg 
    try: 
        for string in msg.split('\n'): 
            # Add appropriate geoprocessing message 
            if severity == 0: 
                arcpy.AddMessage(string) 
            elif severity == 1: 
                arcpy.AddWarning(string) 
            elif severity == 2: 
                arcpy.AddError(string) 
    except: 
        pass

def testAndDelete(fc):
    if arcpy.Exists(fc):
        arcpy.Delete_management(fc)

def numberOfRows(aTable):
    return int(str(arcpy.GetCount_management(aTable)))

def forceExit():
    addMsgAndPrint('Forcing exit by raising ExecuteError')
    raise arcpy.ExecuteError

#################################

def transformAndSwap(ZName,transName,linkFeatures,tanPlunge):
    minY = 0; maxY = 0
    # transform
    if debug:
        addMsgAndPrint('    transforming')
    arcpy.TransformFeatures_edit(ZName,linkFeatures,'SIMILARITY')
    testAndDelete(transName)
    arcpy.Copy_management(ZName,transName)
    # swap coordinates
    if debug:
        addMsgAndPrint('    swapping coordinates')
    geomType = arcpy.Describe(ZName).shapeType
    #addMsgAndPrint('    '+geomType
    if geomType == 'Point': #Point
        with arcpy.da.UpdateCursor(ZName, ['SHAPE@X','SHAPE@Y','SHAPE@Z']) as cursor:
            for row in cursor:
                x = row[0]
                y = row[1]
                z = row[2]
                newY = (z + tanPlunge*y) * cosPlunge
                newZ = y*cosPlunge - z*cosPlunge
                if newY > maxY:
                    maxY = newY
                elif newY < minY:
                    minY = newY
                cursor.updateRow([x,newY,newZ])
    elif geomType in ('Polygon','Polyline'):
        with arcpy.da.UpdateCursor(ZName, ['SHAPE@','OBJECTID']) as cursor:
            for row in cursor:
                #addMsgAndPrint(str(row[0].partCount)+' parts')
                #addMsgAndPrint(str(row[0].pointCount)+ ' points')
                allParts = arcpy.Array()
                try:
                    for part in row[0]:
                        newPart = arcpy.Array()
                        for pnt in part:
                          if pnt:
                            x = pnt.X
                            y = pnt.Y
                            z = pnt.Z
                            newY = (z + tanPlunge*y) * cosPlunge
                            newZ = y
                            if newY > maxY:
                                maxY = newY
                            elif newY < minY:
                                minY = newY
                            pnt = arcpy.Point(x,newY,newZ)
                            newPart.add(pnt)
                        allParts.add(newPart)
                    if geomType == 'Polygon':
                        plg = arcpy.Polygon(allParts)
                    else:
                        plg = arcpy.Polyline(allParts)
                    #addMsgAndPrint(str(plg.pointCount)+' new points')
                    cursor.updateRow([plg,row[1]])
                except:
                    addMsgAndPrint('Problem. Skipping OBJECTID = '+str(row[1]))
    return minY, maxY

def unitVectorPerpSD(strike, dip):
    # returns the unit vector perpendicular to a plane defined by strike and dip
    # formulas after Allmendinger, 2015-2016
    # http://www.geo.cornell.edu/geology/faculty/RWA/structure-lab-manual/chapter-2.pdf    p. 38
    # as of 27 Sept 2019
    strikeRads = math.radians(strike)
    dipRads = math.radians(dip)
    cosA = math.sin(strikeRads) * math.sin(dipRads)
    cosB = -math.cos(strikeRads) * math.sin(dipRads)
    cosC = math.cos(dipRads)
    return (cosA,cosB,cosC)

def angle_between(v1_u,v2_u):
    # assumes v1_u and v2_u are unit vectors!!
    # After David Wolever,
    # https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python/13849249#13849249
    # as of 27 Sept 2019
    #   Returns the angle in radians between vectors 'v1' and 'v2'::
    #   >>> angle_between((1, 0, 0), (0, 1, 0))
    #   1.5707963267948966
    #   >>> angle_between((1, 0, 0), (1, 0, 0))
    #   0.0
    #   >>> angle_between((1, 0, 0), (-1, 0, 0))
    #   3.141592653589793
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def strikelineVector(strike):
    # assumes input strike is in degrees
    return (math.cos(math.radians(strike)), math.sin(math.radians(strike)), 0)

def mirrorAzimuth(azi):
    newAzi = 360-azi
    if newAzi > 180:
        newAzi = newAzi-180
    return newAzi

def dipsToRight(xsAzi,opAzi):
    #### this may not work for some orientations of cross section
    # xsAzi is azimuth of right-hand end of cross section
    # opAzi is azimuth of a point planar feature
    # returned value = True if dip is to right, False if dip is to left
    if opAzi < xsAzi or opAzi > xsAzi+180:
        return True
    else:
        return False

#################################

addMsgAndPrint(versionString)

# validate inputs
for x in xsLine, dem, gdb:
    addMsgAndPrint(x+'  '+str(arcpy.Exists(x)))
    if arcpy.Exists(x) == False:
        forceExit()

# in case DEM is in feet and GeologicMap is in meters, or vice versa...
##  Note that we assume DEM XY units = Z units
demSpatialReference = arcpy.Describe(dem).spatialReference
gMapSpatialReference = arcpy.Describe(gdb+'/GeologicMap').spatialReference
zFactor = demSpatialReference.metersPerUnit/gMapSpatialReference.metersPerUnit
       
# calculate azimuth of xs line
srf = arcpy.Describe(xsLine).spatialReference
for row in arcpy.da.SearchCursor(xsLine,["SHAPE@"]):
    lineGeom = row[0]

xsStart = arcpy.PointGeometry(lineGeom.firstPoint,srf)
xsEnd = arcpy.PointGeometry(lineGeom.lastPoint,srf)
xsAzimuth,xsLengthM = xsStart.angleAndDistanceTo(xsEnd, method = "PLANAR")

addMsgAndPrint("  section azimuth = "+str(xsAzimuth))
addMsgAndPrint("  section length = "+str(xsLengthM))
addMsgAndPrint("  sampleDistance = "+str(sampleDistance))
addMsgAndPrint("  tanPlunge = "+str(tanPlunge))
addMsgAndPrint("  cosPlunge = "+str(cosPlunge))

# calculate link features (FROM end of xs line to (0,0),
#    TO end of xs line to (lengthOfXsLine, 0) )
addMsgAndPrint('Creating link features')
fromPoint1 = lineGeom.firstPoint
toPoint1 = arcpy.Point(0.0,0.0)
fromPoint2 = lineGeom.lastPoint
toPoint2 = arcpy.Point(xsLengthM/srf.metersPerUnit,0.0)
link1 = arcpy.Polyline(arcpy.Array([fromPoint1,toPoint1]),srf)
link2 = arcpy.Polyline(arcpy.Array([fromPoint2,toPoint2]),srf)
linkFeatures = [link1,link2]

# get list of fc to be projected
addMsgAndPrint('Getting list of feature classes to be projected')
projectFCS = []
arcpy.env.workspace = gdb+'/GeologicMap'
fcRoot = arcpy.env.workspace+'/'
for f in ('Point','Line','Polygon'):
    fcs = arcpy.ListFeatureClasses(feature_type=f)
    for fc in fcs:
        #addMsgAndPrint(str(fc)
        if fc.find('edit_') < 0 and fc.find('errors_') < 0:
            projectFCS.append(fcRoot+fc)

# make new feature dataset
xsFDS = gdb+'/CrossSection'+token
if not arcpy.Exists(xsFDS):
    addMsgAndPrint('Making new feature dataset CrossSection' + token)
    arcpy.CreateFeatureDataset_management(gdb,'CrossSection'+token,srf)
else:
    addMsgAndPrint('Feature dataset CrossSection'+token+' already exists')

# buffer xsLine to get clip poly
addMsgAndPrint('Buffering xsLine to get clip poly')
clipPoly = xsFDS+'/'+token+'clipPoly'
testAndDelete(clipPoly)
arcpy.Buffer_analysis(xsLine,clipPoly,bufferWidth,'FULL','ROUND')

addMsgAndPrint('Processing feature classes')
arcpy.CheckOutExtension("3D")
allMinY = 0; allMaxY = 0
for fc in projectFCS:
  #if fc.find('OrientationPoints') > -1:  # if uncommented, restricts to OrientationPoints, for debug purposes only
    shortName = os.path.basename(fc)
    addMsgAndPrint('  '+fc)
    clippedName = xsFDS+'/clipped_'+shortName
    transName = xsFDS+'/trans_'+shortName
    ZName = xsFDS+'/CS'+token+shortName
    for nm in clippedName,ZName:
        testAndDelete(nm)
    # clip
    if debug:
        addMsgAndPrint('    clipping')
    arcpy.Clip_analysis(fc,clipPoly,clippedName)
    if numberOfRows(clippedName) == 0:
        if debug:
            addMsgAndPrint('    empty output')
        testAndDelete(clippedName)
    else:
        if debug:
            addMsgAndPrint('    clipped dataset has '+str(numberOfRows(clippedName))+' rows')
        # enable Z
        if debug:
            addMsgAndPrint('    adding Z to make '+os.path.basename(ZName))
        arcpy.InterpolateShape_3d(dem,clippedName,ZName,sampleDistance,zFactor)
        

        if os.path.basename(fc) == xsLine:  ### this isn't very robust!  
            if debug:
                addMsgAndPrint('    **Making line segments of surface mapunit polys**')
            SMUL = xsFDS+'/CS'+token+'SurfaceMapUnitLines'
            testAndDelete(SMUL)
            arcpy.Intersect_analysis([ZName,gdb+'/GeologicMap/MapUnitPolys'],SMUL)
            minY, maxY = transformAndSwap(SMUL,xsFDS+'/trans_SurfaceMapUnitLines',linkFeatures,tanPlunge)
            testAndDelete(xsFDS+'/trans_SurfaceMapUnitLines')
            arcpy.RecalculateFeatureClassExtent_management(SMUL)
        
            
            if debug:
                addMsgAndPrint('    **Making points where section line crosses contacts**')
                addMsgAndPrint(ZName)
                addMsgAndPrint('  hasZ = '+str(arcpy.Describe(ZName).hasZ))          

            SCAF1 = xsFDS+'/CS'+token+'SurfaceCAF_pts1'
            SCAF = xsFDS+'/CS'+token+'SurfaceCAF_pts'
            testAndDelete(SCAF1); testAndDelete(SCAF)
            arcpy.Intersect_analysis([ZName,gdb+'/GeologicMap/ContactsAndFaults'],SCAF1,'ALL','#','POINT')
            addMsgAndPrint(SCAF1+', # rows = '+str(numberOfRows(SCAF1)))
            addMsgAndPrint('  hasZ = '+str(arcpy.Describe(SCAF1).hasZ))
            arcpy.InterpolateShape_3d(dem,SCAF1,SCAF)
            addMsgAndPrint(SCAF+', # rows = '+str(numberOfRows(SCAF)))
            addMsgAndPrint('  hasZ = '+str(arcpy.Describe(SCAF).hasZ))
            minY, maxY = transformAndSwap(SCAF,xsFDS+'/trans_SurfaceCAF_pts',linkFeatures,tanPlunge)
            testAndDelete(SCAF1)
            testAndDelete(xsFDS+'/trans_SurfaceCAF_pts')
            arcpy.RecalculateFeatureClassExtent_management(SCAF)
            
            
        minY, maxY = transformAndSwap(ZName,transName,linkFeatures,tanPlunge)    
        if minY < allMinY:
            allMinY = minY
        if maxY > allMaxY:
            allMaxY = maxY
        arcpy.RecalculateFeatureClassExtent_management(ZName)
        if os.path.basename(fc) == xsLine:
            if debug:
                addMsgAndPrint('    **Making surface profile**')
            testAndDelete(xsFDS+'/CS'+token+'TopographicProfile')
            arcpy.Copy_management(ZName,xsFDS+'/CS'+token+'TopographicProfile')
            
        testAndDelete(clippedName)
        testAndDelete(transName)
testAndDelete(clipPoly)
arcpy.CheckInExtension('3D')

addMsgAndPrint('minY = '+str(minY))
addMsgAndPrint('maxY = '+str(maxY))

#  do OrientationPoints
orp = xsFDS+'/CS'+token+'OrientationPoints'

xsNormal = unitVectorPerpSD(xsAzimuth,90-xsPlunge)
xsStrikelineVector = strikelineVector(xsAzimuth)

if debug2:
    addMsgAndPrint('xsNormal = '+str(xsNormal))
if arcpy.Exists(orp):
    addMsgAndPrint('Calculating apparent dips and stuff')
    ## add fields
    arcpy.AddField_management(orp,'DistanceFromSection','FLOAT')
    arcpy.AddField_management(orp,'MapAzimuth','FLOAT')
    arcpy.AddField_management(orp,'IntersectAngle','FLOAT')
    arcpy.AddField_management(orp,'ApparentDip','FLOAT')
    arcpy.AddField_management(orp,'Rake','FLOAT')
    arcpy.CalculateField_management(orp,'MapAzimuth','!Azimuth!','Python')
    ## recalc some fields
    ## step through points
    fields = ['Type','Azimuth','Inclination','DistanceFromSection','IntersectAngle','SHAPE@Z','ApparentDip','Rake','OrientationPoints_ID']
    with arcpy.da.UpdateCursor(orp, fields) as cursor:
        for row in cursor:
            mapAzi = row[1]
            mapInc = row[2]
            if debug2:
                addMsgAndPrint(row[7])
                addMsgAndPrint('XS azi = '+str(xsAzimuth)+'  Map azi = '+str(mapAzi)+'  Map inc = '+str(mapInc))
            
            dfs = row[5]
            obliq = mapAzi - xsAzimuth
            if obliq < 0:
                obliq = obliq+360

            appDip = math.atan(math.sin(obliq) * math.tan(math.radians(row[2])))
            # next two lines flatten for down-plunge projectoion
            appDip = math.atan(math.tan(appDip) * cosPlunge)
            appDip = math.degrees(appDip)
            if obliq <= 180:
                plotAz = 270 - abs(appDip)
            else:
                plotAz = 270 + abs(appDip)
                
            if plotAz >= 360:
                plotAz = plotAz - 360.0

            attitudeNormal = unitVectorPerpSD(mapAzi,mapInc)
            
            # alternate apparent dip calculation
            intersection = np.cross(xsNormal,attitudeNormal)
            rake = math.degrees(angle_between(intersection,xsStrikelineVector))
            plotAz = rake+90
            if plotAz > 180:
                plotAz = plotAz - 180.0

            if dipsToRight(xsAzimuth,mapAzi):
                plotAz = mirrorAzimuth(plotAz)
                
            
            if debug2:
                addMsgAndPrint('app dip = '+str(appDip)+'  plotAzi = '+str(plotAz)+'  obliquity = '+str(obliq))
                if plotAz < 180:
                    addMsgAndPrint('   oops!')
                elif plotAz <= 270:
                    addMsgAndPrint('   dip to left')
                else:
                    addMsgAndPrint('   dip to right')
                 
            angleBetween = math.degrees(angle_between(xsNormal,attitudeNormal))
            # correct from vectors to axes
            if angleBetween > 90:
                angleBetween = 180 - angleBetween
            if debug2:
                addMsgAndPrint('angleBetween = '+str(angleBetween))
                addMsgAndPrint(' ')

            cursor.updateRow([row[0],plotAz,row[2],dfs,angleBetween,row[5],appDip,rake,row[7]])
    del cursor
    

# build framing lines
addMsgAndPrint('Adding framing lines')
## make sure fc is there and empty
cl = xsFDS+'/CS'+token+'CartographicLines'
if arcpy.Exists(cl):
    # delete all features
    arcpy.DeleteFeatures_management(cl)
    cursor = arcpy.da.InsertCursor(cl, ['Type','Label','SHAPE@'])
    ## z = 0
    addMsgAndPrint('  adding horizon line')
    lin = arcpy.Polyline(arcpy.Array([toPoint1,toPoint2]))
    cursor.insertRow(['horizon','',lin])

    leftX = toPoint1.X
    rightX = toPoint2.X
    ## left end and right end
    addMsgAndPrint('  adding end caps')
    minYindex = -10
    maxYindex = 6
    minY = minYindex * 1000
    maxY = maxYindex * 1000
    ticLength = 200
    lin = arcpy.Polyline(arcpy.Array([arcpy.Point(leftX,minY),arcpy.Point(leftX,maxY)]))
    cursor.insertRow(['end cap','',lin])
    lin = arcpy.Polyline(arcpy.Array([arcpy.Point(rightX,minY),arcpy.Point(rightX,maxY)]))
    cursor.insertRow(['end cap','',lin])
    
    #3 elevation tics
    addMsgAndPrint('  adding elevation tics')
    for yi in range(minYindex,maxYindex+1):
        y = yi*1000
        p1 = arcpy.Point(leftX-ticLength,y)
        p2 = arcpy.Point(leftX,y)
        lin = arcpy.Polyline(arcpy.Array([p1,p2]))
        cursor.insertRow(['left tic',str(y),lin])

        p1 = arcpy.Point(rightX+ticLength,y)
        p2 = arcpy.Point(rightX,y)
        lin = arcpy.Polyline(arcpy.Array([p1,p2]))
        cursor.insertRow(['right tic',str(y),lin])
        
        
    del cursor

addMsgAndPrint('DONE')
if ForceExit == 'true':
    forceExit()



