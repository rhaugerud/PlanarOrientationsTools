#script to project map data into a cross section
# Ralph Haugerud, U.S. Geological Survey
# Code not subject to U.S. copyright

import arcpy, sys, os.path, math
import numpy as np

versionString = 'DownPlungeProjection2.py, version of 7 May 2021'


#### Still to be done
## test for linear (not-planar) orientationpoints features (DONE) and do a different calculation (NOT DONE)
## use allMinY and allMaxY to calculate vertical extent of cartographic frame
## check that tool help matches changes I've made
## append a note about inputs to CrossSectionNotes.txt in host .gdb. Should include fds name, date, code version, domain-name, and projection axis

'''
gdb    # GeMS database
domLyr
axisTrendField
axisPlungeField
dem    # coextensive with area to be projected x,y,z units are same as those of GeologicMap 
numberOfCells  # sampling distance, as cellSize multiplier
token  # 'CS'+token = prefix for new feature classes
       # 'CrossSection'+token = name of new feature dataset

ASSUMES
    DEM and input geology have same spatial reference framework
    XY units = Z units
'''
debug = False
debug2 = False

gdb = sys.argv[1]
domLyr = sys.argv[2]
axisTrendField = sys.argv[3]
axisPlungeField = sys.argv[4]
dem = sys.argv[5]
numberOfCells = float(sys.argv[6])
token = sys.argv[7]
ForceExit = sys.argv[8]


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

def fieldNameList(aTable):
    fns = arcpy.ListFields(aTable)
    fns2 = []
    for fn in fns:
        fns2.append(fn.name)
    return fns2

def forceExit():
    addMsgAndPrint('Forcing exit by raising ExecuteError')
    raise arcpy.ExecuteError

#################################

def transformAndSwap(ZName,transName,linkFeatures,tanPlunge,xOffset):
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
                x = row[0] + xOffset
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
                            x = pnt.X + xOffset
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

def isPlanar(orpType):
    if orpType == None:
        return False
    isp = False
    for x in ('bedding','foliation','cleavage','joint','layering','plane','schistosity'):
        if orpType.lower().find(x) > -1:
            isp = True
    return isp

def unitVectorPerpSD(strike, dip):
    # returns the unit vector perpendicular to a plane defined by strike and dip
    # formulas after Allmendinger, 2015-2016
    # http://www.geo.cornell.edu/geology/faculty/RWA/structure-lab-manual/chapter-2.pdf    p. 38
    # as of 27 Sept 2019
    ## A,B,C = North, East, Down
    ## rest of my code uses East,North,Up as axes, thus the rearrangement in the return statement
    strikeRads = math.radians(strike)
    dipRads = math.radians(dip)
    cosA = math.sin(strikeRads) * math.sin(dipRads)  # north
    cosB = -math.cos(strikeRads) * math.sin(dipRads)# east
    cosC = math.cos(dipRads)  # down
    return (cosB,cosA,0-cosC)

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


def unitVector(trend,plunge):
    #assumes input trend and plunge are in degrees
    rPlunge = math.radians(plunge)
    rTrend = math.radians(90-trend)
    z = 0-math.sin(rPlunge)
    zMult = math.cos(rPlunge)
    x = zMult*math.cos(rTrend)
    y = zMult*math.sin(rTrend)
    #print 'unitVector length = ',math.sqrt(x*x+y*y+z*z)
    return x,y,z

def projectdownplunge(axisTrend,axisPlunge,strike,dip):

    xsAzimuth = (axisTrend + 90) % 360
    xsNormal = unitVector(axisTrend,axisPlunge)
    attitudeNormal = unitVectorPerpSD(strike,dip)
    xsStrikelineVector = unitVector((axisTrend+90) % 360,0)
 
    intersec = np.cross(xsNormal,attitudeNormal)
    # scale intersection to unit vector
    xLength = math.sqrt(intersec[0]*intersec[0]+intersec[1]*intersec[1]+intersec[2]*intersec[2])
    intersection = [intersec[0]/xLength,intersec[1]/xLength,intersec[2]/xLength]
    rake = math.degrees(angle_between(intersection,xsStrikelineVector))

    if debug2:
        print 'xsAzimuth = ',xsAzimuth
        print 'xsNormal = ',xsNormal
        print 'attitudeNormal = ',attitudeNormal
        print 'xsStrikelineVector = ',xsStrikelineVector
        print 'intersection = ',intersection
        print 'rake = ',rake
    print 'rake = ',rake
           
    if abs(axisTrend-strike)<= 90:  # dips to right
        print 'dips to right'
        plotAz = (90 - rake) % 360
    else:  # dips to left
        print 'dips to left'
        plotAz = (90 + rake) % 360
        
    obliq = math.degrees(angle_between(attitudeNormal, xsNormal))
    if obliq > 90:
        obliq = 180-obliq
    
    if plotAz > 0 and plotAz <= 180:
        plotAz = plotAz+180
         
    return plotAz, obliq,


#################################

addMsgAndPrint(versionString)

# validate inputs
for x in domLyr, dem, gdb:
    if arcpy.Exists(x) == False:
        addMsgAndPrint(x+' does not exist. Exiting')
        forceExit()

sampleDistance = numberOfCells * arcpy.Describe(dem).meanCellWidth
srf = arcpy.Describe(domLyr).spatialReference

##  Note that we assume DEM XY units = DEM Z units = Geology XY units
## zFactor = DEM Z metersPerUnit / DEM XY metersPerUnit
zFactor = 1

# make new feature dataset
xsFDS = gdb+'/CrossSection'+token
if not arcpy.Exists(xsFDS):
    addMsgAndPrint('Making new feature dataset CrossSection' + token)
    arcpy.CreateFeatureDataset_management(gdb,'CrossSection'+token,srf)
else:
    addMsgAndPrint('Feature dataset CrossSection'+token+' already exists')

#@# extract trend and plunge and length of xs line from domLyr, axisTrendField, and axisPlungeField
fields = ["SHAPE@TRUECENTROID",axisTrendField,axisPlungeField]
addMsgAndPrint('Layer '+domLyr+' has '+str(numberOfRows(domLyr))+' row(s)')
for row in arcpy.da.SearchCursor(domLyr,fields):
    [[centroidX,centroidY],axisTrend,axisPlunge] = row

if axisTrend == None or axisPlunge == None:
    addMsgAndPrint('Null value in '+axisTrendField+' or '+axisPlungeField+'.  Halting')
    forceExit()
    

tanPlunge = math.tan(math.radians(axisPlunge))
cosPlunge = math.cos(math.radians(axisPlunge))
cosTrendGeom = math.cos(math.radians(90-axisTrend))
sinTrendGeom = math.sin(math.radians(90-axisTrend))
xsAzimuth = (axisTrend + 90) % 360

if debug2:
    addMsgAndPrint('xsAzimuth = '+str(xsAzimuth)) 
    
#@# figure out length of xs line and get links
fromPoint1 = arcpy.Point(centroidX,centroidY)
toPoint1 = arcpy.Point(0,0)
fromPoint2 = arcpy.Point(centroidX + cosTrendGeom*10, centroidY + sinTrendGeom*10)
toPoint2 = arcpy.Point(0,10)

link1 = arcpy.Polyline(arcpy.Array([fromPoint1,toPoint1]),srf)
link2 = arcpy.Polyline(arcpy.Array([fromPoint2,toPoint2]),srf)
linkFeatures = [link1,link2]

shortName = os.path.basename(domLyr)
transName = xsFDS+'/trans_'+shortName
ZName = xsFDS+'/CS'+token+shortName
for nm in transName,ZName:
    testAndDelete(nm)

addMsgAndPrint('  adding Z values to domain boundary')
arcpy.CheckOutExtension("3D")

arcpy.InterpolateShape_3d(dem,domLyr,ZName,sampleDistance,zFactor)
addMsgAndPrint('  transforming domain boundary into XS space')
xOffset = 0.0
minY, maxY = transformAndSwap(ZName,transName,linkFeatures,tanPlunge, xOffset)

# scan transName for minX and MaxX
addMsgAndPrint('  getting xOffset')
for row in arcpy.da.SearchCursor(transName,['SHAPE@']):
    minX = maxX = row[0][0][0].X
    for part in row[0]:
        for pnt in part:
            if pnt.X < minX:
                minX = pnt.X
            if pnt.X > maxX:
                maxX = pnt.X
xOffset = 0.0 - minX
testAndDelete(transName)


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
            addMsgAndPrint('  '+str(fc))

clipPoly = domLyr

addMsgAndPrint('Processing feature classes')
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
            
        minY, maxY = transformAndSwap(ZName,transName,linkFeatures,tanPlunge,xOffset)    
        if minY < allMinY:
            allMinY = minY
        if maxY > allMaxY:
            allMaxY = maxY
        arcpy.RecalculateFeatureClassExtent_management(ZName)

        testAndDelete(clippedName)
        testAndDelete(transName)
arcpy.CheckInExtension('3D')

addMsgAndPrint('allMinY = '+str(allMinY))
addMsgAndPrint('allMaxY = '+str(allMaxY))

#  do OrientationPoints
addMsgAndPrint('Calculating apparent dips')
arcpy.env.workspace = xsFDS
fcs = arcpy.ListFeatureClasses(feature_type='Point')
addMsgAndPrint(arcpy.env.workspace)
addMsgAndPrint(str(fcs))
for fc in fcs:
    fcFields = fieldNameList(fc)
    addMsgAndPrint(fc+' has '+str(numberOfRows(fc))+' rows')
    if 'Azimuth' in fcFields and 'Inclination' in fcFields and numberOfRows(fc) > 0:
        orp = fc
        addMsgAndPrint('  calculating apparent dips for '+os.path.basename(orp))
        ## add fields
        for newField in ('DistanceFromSection','MapAzimuth','Obliquity'):
            if not newField in fcFields:
                arcpy.AddField_management(orp,newField,'FLOAT')
        arcpy.CalculateField_management(orp,'MapAzimuth','!Azimuth!','Python')
        ## step through points and calc some fields
        fields = ['Type','Azimuth','Inclination','DistanceFromSection','SHAPE@Z','Obliquity','Type']
        with arcpy.da.UpdateCursor(orp, fields) as cursor:
            for row in cursor:
                mapAzi = row[1]
                mapInc = row[2]
                orpType = row[6]
                if mapAzi <> None and mapInc <> None:
                    if isPlanar(orpType):
                        dfs = row[4]  # note that X-section coordinate axes are X=Right, Y=Up, Z=Away
                        plotAz, obliq = projectdownplunge(axisTrend,axisPlunge,mapAzi,mapInc)
                        if debug2:
                            addMsgAndPrint('strike = '+str(mapAzi)+', dip = '+str(mapInc))
                            addMsgAndPrint('plotAzimuth = '+str(plotAz)+'  obliquity = '+str(obliq))
                        #            0       1          2             3                     4           5       6
                        #fields = ['Type','Azimuth','Inclination','DistanceFromSection','SHAPE@Z','Obliquity','Type']
                        cursor.updateRow([row[0],plotAz,row[2],dfs,row[4],obliq,orpType])
                    else:  # is linear, need a different calculation
                        addMsgAndPrint(orpType+' is not recognized feature. Write some new projection code!')
        del cursor
    

# build framing lines
addMsgAndPrint('Adding framing lines')
## make sure fc is there and empty
cl = xsFDS+'/CS'+token+'CartographicLines'
if not arcpy.Exists(cl):
    arcpy.CreateFeatureclass_management(os.path.dirname(cl),os.path.basename(cl), "POLYLINE")
    for i in [['Type',254],
              ['Symbol',254],
	      ['Label',254],
              ['DataSourceID',50],
	      ['Notes',254],
              [os.path.basename(cl)+'_ID',50]]:
        arcpy.AddField_management(cl,i[0],'TEXT','','',i[1])
else: # cl exists
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



