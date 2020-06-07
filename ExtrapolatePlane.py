# In Arc, calculate the intersection of a plane (location, strike, dip)
#  with topography. Display the results
#
# Ralph Haugerud, U.S. Geological Survey
# Code not subject to U.S. copyright

versionString = 'ExtrapolatePlanes.py, version of 10 February 2020'

import math, os.path, arcpy, sys, shutil
from arcpy.sa import *

arcpy.CheckOutExtension("3D")
arcpy.CheckOutExtension("Spatial")
        

##=======================================
# first 4 functions copied from GeMS_utilityFunctions.py
def testAndDelete(fc):
    if arcpy.Exists(fc):
        arcpy.Delete_management(fc)

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

def forceExit():
    addMsgAndPrint('Forcing exit by raising ExecuteError')
    raise arcpy.ExecuteError

def numberOfRows(aTable):
    return int(str(arcpy.GetCount_management(aTable)))

def makeCoordsList(xyz,azi,inc,radius):
    addMsgAndPrint('  calculating control points')
    aziR = math.radians(azi)
    ddR = math.radians(azi+90)
    incR = math.radians(inc)

    x1 = xyz[0] + radius * math.sin(aziR)
    y1 = xyz[1] + radius * math.cos(aziR)
    z1 = xyz[2]

    x2 = xyz[0] + radius * math.sin(ddR) * math.cos(incR)
    y2 = xyz[1] + radius * math.cos(ddR) * math.cos(incR)
    z2 = xyz[2] - radius * math.sin(incR)

    x3 = xyz[0] - radius * math.sin(aziR)
    y3 = xyz[1] - radius * math.cos(aziR)
    z3 = xyz[2]

    x4 = xyz[0] - radius * math.sin(ddR) * math.cos(incR)
    y4 = xyz[1] - radius * math.cos(ddR) * math.cos(incR)
    z4 = xyz[2] + radius * math.sin(incR)

    coordsList = [[xyz[0],xyz[1],xyz[2]],
              [x1,y1,z1],
              [x2,y2,z2],
              [x3,y3,z3],
              [x4,y4,z4]]
    return coordsList

def makeControlPoints(dem,coordslist,scratchgdb):
    addMsgAndPrint('  making control-point feature class')
    sr = arcpy.Describe(dem).spatialReference
    scp = os.path.join(scratchgdb,'xxxControlPoints')
    testAndDelete(scp)
    arcpy.CreateFeatureclass_management(scratchgdb,'xxxControlPoints','POINT',
                                        '','DISABLED','ENABLED',sr)
    cursor = arcpy.da.InsertCursor(scp,['SHAPE@X','SHAPE@Y','SHAPE@Z'])
    for pt in coordslist:
        cursor.insertRow([pt[0],pt[1],pt[2]])
    del cursor
    return scp

def makeOutcropRaster(dem,controlPoints,inc,sgdb):
    addMsgAndPrint('  making predicted outcrop raster')
    # set up raster environment
    arcpy.env.snapRaster = dem
    arcpy.env.cellSize = arcpy.Describe(dem).meanCellWidth
    arcpy.env.extent = controlPoints
    # TREND
    addMsgAndPrint('    calculating plane with TREND')
    outRaster = os.path.join(sgdb,'xxxPlane')
    testAndDelete(outRaster)
    arcpy.Trend_3d(controlPoints, 'SHAPE', outRaster, '', 1, 'LINEAR')
    # subtract TREND plane from DEM
    addMsgAndPrint('    subtracting plane from DEM')
    intersectRaster1 = outRaster+'_i'
    testAndDelete(intersectRaster1)
    arcpy.Minus_3d(dem, outRaster, intersectRaster1)
    arcpy.CalculateStatistics_management(intersectRaster1)
    # classify results
    ###### need to trap for inclination = 90  ( tan(incR) = infinity )
    dH = math.tan(math.radians(inc)) * 2.0 * float(arcpy.env.cellSize)
    addMsgAndPrint('    classifying results, dH='+str(dH))
    outCon = Con(intersectRaster1, 0, 1, 'VALUE > '+str(dH)+' OR VALUE < -'+str(dH))
    testAndDelete(outRaster)
    testAndDelete(intersectRaster1)
    return outCon

#==================================
    
## define inputs
dem = sys.argv[1]
gridDec = float(sys.argv[2])
inputFC = sys.argv[3]
aziField = sys.argv[4]
incField = sys.argv[5]
radius = float(sys.argv[6])
scratchgdb = sys.argv[7]

inPts = os.path.join(scratchgdb,'inputPoints')

addMsgAndPrint('--')
addMsgAndPrint(versionString)

## copy inputFC to scratchgdb
addMsgAndPrint('Copying input features to feature class '+scratchgdb+'/inputPoints')
testAndDelete(inPts)
arcpy.CopyFeatures_management(inputFC,inPts)
#test for number of values
addMsgAndPrint('  '+str(numberOfRows(inPts))+' rows of input features')
if numberOfRows(inPts) > 5:
    addMsgAndPrint('Too many input features. Select fewer!')
    forceExit()

## if necessary, project inPts to dem srf
demSR = arcpy.Describe(dem).spatialReference
inPtsSR = arcpy.Describe(inPts).spatialReference
if demSR <> inPtsSR:
    addMsgAndPrint('Projecting structure points from WKID '+str(inPtsSR.factoryCode)+' to WKID '+str(demSR.factoryCode))
    # rename inPts
    inPts2 = inPts+'_2'
    testAndDelete(inPts2)
    arcpy.Rename_management(inPts,inPts2)
    # project inPts2
    arcpy.Project_management(inPts2,inPts,demSR)
    testAndDelete(inPts2)
    

## Add XYZ values
addMsgAndPrint('Getting XYZ values')
arcpy.AddXY_management(inPts)
arcpy.AddSurfaceInformation_3d(inPts,dem,'Z')

## for each orientation point
# extract xyz, azi, inc
# validate azi and inc values
with arcpy.da.SearchCursor(inPts,['OID@','POINT_X','POINT_Y','Z',aziField,incField]) as cursor:
    for row in cursor:
        oid = row[0]
        xyz = [row[1],row[2],row[3]]
        azi = row[4] + gridDec  ############### or is it minus?
        inc = row[5]
        addMsgAndPrint('OID='+str(oid)+' strike='+str(azi)+' dip='+str(inc))
        #addMsgAndPrint(str(xyz))
        coordsList = makeCoordsList(xyz, azi,inc,radius)
        controlPoints = makeControlPoints(dem,coordsList,scratchgdb)
        outCon = makeOutcropRaster(dem,controlPoints,inc,scratchgdb)
        testAndDelete(controlPoints)
        planeName = 'OID'+str(oid)+'_'+str(int(azi))+'_'+str(int(inc))
        saveRaster = os.path.join(scratchgdb,planeName)
        testAndDelete(saveRaster)
        outCon.save(saveRaster)
        ## import existing layer file, reset source, rename to source
        addMsgAndPrint('  inserting new layer')
        lyrFile = os.path.join(os.path.dirname(sys.argv[0]),'projectedBedding.lyr')
        newLyrFile = os.path.join(os.path.dirname(scratchgdb),planeName+'.lyr')
        # copy to newname in this workspace
        shutil.copy(lyrFile,newLyrFile)
        inLayer = arcpy.mapping.Layer(newLyrFile)
        inLayer.replaceDataSource(scratchgdb,'FILEGDB_WORKSPACE',planeName)
        inLayer.name = planeName
        mxd = arcpy.mapping.MapDocument("CURRENT")
        dataFrameList = arcpy.mapping.ListDataFrames(mxd)
        if len(dataFrameList) > 1:
            addMsgAndPrint('This script assumes the map document has only data frame')
            addMsgAndPrint('whereas this map document has '+str(len(dataFrameList)))
            addMsgAndPrint('ABORTING')
            forceExit()
        currentFrame = dataFrameList[0]

        arcpy.mapping.AddLayer(currentFrame,inLayer)
        # clean up
        del mxd
        #testAndDelete(newLyrFile)

testAndDelete(inPts)
addMsgAndPrint('DONE!')
addMsgAndPrint('=====')

forceExit()

print 'DONE!'




