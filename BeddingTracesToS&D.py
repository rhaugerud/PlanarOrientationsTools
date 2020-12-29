# calculate strike and dip from bedding trace
# Ralph Haugerud, U.S. Geological Survey
# Code not subject to U.S. copyright

versionString = 'BeddingTracesToS&D.py, version of 8 April 2020'
# 6 April 2020: Changed toolbox interface to use raster layer instead of raster dataset
#     changed scatter and spread fields to Scatter and Spread
# 29 December 2020:  Added optional Type and Symbol values


import arcpy, sys, os.path, math
import numpy as np
from datetime import datetime

densifyFactor = 3  # densification factor: input traces are densified by densifyFactor * DEM cellsize
sgdbName = 'xxx_S&D.gdb'

###### FIT FILTERS ######
#   Choices may not be entirely justifiable
#   Note that values of maxScatter and minSpread are multipliers for DEM cell size
minN = 15          # minimum number of observations, after densification
testOb = 1.2       # minimum acceptable oblateness
maxScatter = 3     # maximum acceptable value of sqrt(eigenvalue3/n) in units of DEM cell size
minSpread = 3      # minimum acceptable value of sqrt(eigenvalue2/n) in units of DEM cell size


##########Utility Functions####################
def addMsgAndPrint(msg, severity=0): 
    # prints msg to screen and adds msg to the geoprocessor (in case this is run as a tool) 
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
        if debug1: addMsgAndPrint('    deleting '+fc)
        arcpy.Delete_management(fc)
    return

def nRows(aTable):
    return int(arcpy.GetCount_management(aTable).getOutput(0))

def forceExit():
    addMsgAndPrint('Forcing exit by raising ExecuteError')
    raise arcpy.ExecuteError

def fieldNameList(aTable):
    fns = arcpy.ListFields(aTable)
    fns2 = []
    for fn in fns:
        fns2.append(fn.name)
    return fns2

############Script-specific Functions##################

def planeFit(xyz):
    # xyz is python list of point positions: [[x1,y1,z1],[x2,y2,z2],[x3,y3,z3],...]
    # must have at least 3 elements!
    # returns strike,dip,eigenVals,eigenVects, centroid, n
    ## convert to numpy array and transpose, so that columns are points
    xyz1 = np.array(xyz).T
    ## get mean
    mX = np.mean(xyz1[0,:])
    mY = np.mean(xyz1[1,:])
    mZ = np.mean(xyz1[2,:])
    mean_vector = np.array([[mX],[mY],[mZ]])
    ## compute scatter matrix (xyz2)
    xyz2 = np.zeros((3,3))
    for i in range(xyz1.shape[1]):
        xyz2 += (xyz1[:,i].reshape(3,1) - mean_vector).dot((xyz1[:,i].reshape(3,1) - mean_vector).T)
    ## calculate eigenvalues and eigenvectors
    eigVal,eigVec = np.linalg.eig(xyz2)
    ## make eigenvalue-eigenvector pairs and sort
    eigPairs = [(np.abs(eigVal[i]),eigVec[:,i]) for i in range(len(eigVal))]
    eigPairs.sort()
    eigPairs.reverse()
    # smallest eigenvector is normal to best-fit plane
    normal = eigPairs[2][1]
    # if normal is in upper hemisphere, take its opposite
    if normal[2] > 0:
        for i in 0,1,2:
            normal[i] = normal[i]*-1.0
    # calculate dip
    ## if normal[2] = 0, dip is 90
    ## if normal[2] = -1, dip is 0
    dip = math.degrees(math.acos(0-normal[2])) 
    # dip direction
    if dip == 0:
        strike = 0
    else:
        if normal[0] <> 0:
            dipDirectionC = ( math.degrees(math.atan2(normal[1],normal[0])) + 180 ) % 360
        else: # normal[0] == 0
            if normal[1] > 0:
                dipDirectionC = 270
            elif normal[1] <  0:
                dipDirectionC = 90
            else: # normal[0] = 0, normal[1] = 0, normal[2] = -1
                dipDirectionC = 0                
        #dipDirectionG = 90 - dipDirectionC
        #strike = dipDirectionG - 90
        #strike = ( strike + 360 ) % 360
        # above 3 lines add up to
        strike = ( 360 - dipDirectionC ) % 360
    eigenVects = [eigPairs[0][1],eigPairs[1][1],eigPairs[2][1]]
    eigenVals = [eigPairs[0][0],eigPairs[1][0],eigPairs[2][0]]
    if eigenVals[1] == 0:
        addMsgAndPrint('*** eigenVals[1] = 0')
        addMsgAndPrint('*** substituting value of 0.00001 to avoid divide-by-zero errors')
        eigenVals[1] = 0.00001
    if eigenVals[2] == 0:
        addMsgAndPrint('*** eigenVals[2] = 0')
        addMsgAndPrint('*** substituting value of 0.00001 to avoid divide-by-zero errors')
        eigenVals[2] = 0.00001
    centroid = [mX,mY,mZ]
    n = len(xyz)
    return strike,dip,eigenVals,eigenVects,centroid,n


def distance(xy1,xy2):
    return math.sqrt((xy1[0]-xy2[0])**2+(xy1[1]-xy2[1])**2)

def makeXyzList(pointclass,stackLines):
    xyzList = []
    xyzListStacked = []
    fields = ['SHAPE@']
    sumX = 0; sumY = 0; sumZ = 0; nRows = 0
    for row in arcpy.da.SearchCursor(pointclass, fields):
        #addMsgAndPrint(str(row[1])+','+str(row[2])+','+str(row[3]))
        #midPoints.append([row[1],row[2],row[3]])
        partXSum = 0; partYSum = 0; partZSum = 0
        nParts = 0
        for part in row[0]:
            xSum = 0; ySum = 0; zSum = 0; n = 0
            for pnt in part:
                xSum = xSum+pnt.X
                ySum = ySum+pnt.Y
                zSum = zSum+pnt.Z
                n = n+1
            partX = xSum/n; partY = ySum/n; partZ = zSum/n
            nParts = nParts+1
            partXSum = partXSum+partX; partYSum = partYSum+partY; partZSum = partZSum+partZ
        rowX = partXSum/nParts; rowY = partYSum/nParts; rowZ = partZSum/nParts
        for part in row[0]:
            for pnt in part:
                if stackLines:
                    xyzListStacked.append([pnt.X-rowX,pnt.Y-rowY,pnt.Z-rowZ])
                xyzList.append([pnt.X,pnt.Y,pnt.Z])
        nRows = nRows+1
        sumX = sumX+rowX; sumY = sumY+rowY; sumZ = sumZ+rowZ
    meanXYZ = [sumX/nRows,sumY/nRows,sumZ/nRows]
    if nRows == 1:  # find closest input point
        closest = xyzList[0]
        dMin = distance(meanXYZ,xyzList[0])
        for i in xyzList:
            d = distance(meanXYZ,i)
            if d < dMin:
                closest = i
                dMin = d
        meanXYZ = closest
        
    if stackLines:
        return xyzListStacked,meanXYZ
    else:
        return xyzList,meanXYZ

def gallo(ob):
    # calculates OrientationConfidenceDegrees from Oblateness
    # Gallo et al. (JGR:Solid Earth, 2018) suggest
    #    O = ln(ev2/ev3) = Oblateness
    #    O > 5 = OrientationConfidenceDegrees <= 5
    #    O > 1.5 = OrConfDeg <= 10
    #  Code below is my continuous extrapolation from their Figure 6
    #  Note that this function perhaps returns too-low values for ob < 1.2
    if ob >= 1.5:
        return 10**(1.30103 - 0.200687 *ob)
    elif ob < 1.5 and ob >= 1:
        return 10**(2.2 - 0.8 * ob)
    else:
        addMsgAndPrint('Oblateness <= 1. OCD not defined')
        forceExit()
    
#####################################

dem = sys.argv[1]
gridDec = float(sys.argv[2])      
bt = sys.argv[3]
orPts = sys.argv[4]
scratchFolder = sys.argv[5]
dsID = sys.argv[6]
if sys.argv[7].upper() == 'TRUE':
    stackLines = True
else:
    stackLines = False
if sys.argv[8] <> '#':
    calcType = sys.argv[8]
else:
    calcType = None
if sys.argv[9] <> '#':
    calcSymbol = sys.argv[9]
else:
    calcSymbol = None

dens = densifyFactor * arcpy.Describe(dem).meanCellHeight
maxScatter = maxScatter * arcpy.Describe(dem).meanCellHeight
minSpread = minSpread * arcpy.Describe(dem).meanCellHeight

debug1 = False

addMsgAndPrint(versionString)

## check for required fields in orPts
for xx in ('Type','Azimuth','Inclination','OrientationConfidenceDegrees',
           'OrientationSourceID','Notes','Oblateness','Spread','Scatter'):
    if not xx in fieldNameList(orPts):
        addMsgAndPrint('Please add field "'+xx+'" to '+orPts)
        addMsgAndPrint('    you may do this manually or')
        addMsgAndPrint('    you may use script tool Prep OrientationPoints')
        forceExit()

if 'Symbol' in fieldNameList(orPts):
    hasSymbol = True
else:
    hasSymbol = False

## make sgdb if it doesn't already exist
sgdb = scratchFolder +'/'+sgdbName
if not arcpy.Exists(sgdb):
    arcpy.CreateFileGDB_management(scratchFolder,sgdbName)

## Get list of OBJECTIDS
addMsgAndPrint('Getting list of OBJECTIDs')
iDs = "OIDs:"
for row in arcpy.da.SearchCursor(bt,["OBJECTID"]):
    iDs= iDs+' '+str(row[0])
if stackLines:
    iDs = iDs+' stacked'
if len(iDs) > 254:
    iDs = iDs[0:254]


## make copy of (selected arcs in) scratchgdb
xlbt = os.path.join(sgdb,'inputArcs')
arcpy.env.outputZFlag = 'Enabled'
addMsgAndPrint('Copying input arcs to feature class '+xlbt)
testAndDelete(xlbt)
arcpy.CopyFeatures_management(bt,xlbt)

## if necessary, project inPts to dem srf
demSR = arcpy.Describe(dem).spatialReference
xlbtSR = arcpy.Describe(xlbt).spatialReference
if demSR <> xlbtSR:
    addMsgAndPrint('Projecting input arcs from WKID '+str(xlbtSR.factoryCode)+' to WKID '+str(demSR.factoryCode))
    # rename inPts
    xlbt2 = xlbt+'_2'
    testAndDelete(xlbt2)
    arcpy.Rename_management(xlbt,xlbt2)
    # project inPts2
    arcpy.Project_management(xlbt2,xlbt,demSR)
    testAndDelete(xlbt2)

## add z values 
addMsgAndPrint('Adding Z values')
arcpy.CheckOutExtension('3D')
xlbtdz = xlbt+'_dz'
testAndDelete(xlbtdz)
# note that we densify
arcpy.InterpolateShape_3d(dem,xlbt,xlbtdz,dens,1,'BILINEAR','DENSIFY')
if nRows(xlbtdz) == 0:
    print '  Bedding trace(s) appear to be outside extent of DEM'
    forceExit()

##  Make point lists
addMsgAndPrint('Building point lists')
xyzListd,meanXYZd = makeXyzList(xlbtdz,stackLines)

## Fit a plane
addMsgAndPrint('Fitting a plane')
strike,dip,eigenVals,eigenVects,centroid,n = planeFit(xyzListd)
strike = strike + gridDec
# Fernandez, 2005; Woodcock, 1977, Gallo and others (2018)
M = math.log(eigenVals[0]/eigenVals[2])
K = math.log(eigenVals[0]/eigenVals[1])/math.log(eigenVals[1]/eigenVals[2])
oblateness = math.log(eigenVals[1]/eigenVals[2])
OCD = round(gallo(oblateness),1)
scatter = math.sqrt(eigenVals[2]/n)
spread = math.sqrt(eigenVals[1]/n)
    
addMsgAndPrint('Densified, n = '+str(n))
addMsgAndPrint('  strike/dip = {:.0f}'.format(strike)+'/{:.0f}'.format(dip))
addMsgAndPrint('  eigenvalues: {:.2f}'.format(eigenVals[0])+', {:.2f}'.format(eigenVals[1])+', {:.2f}'.format(eigenVals[2]))
addMsgAndPrint('  M = {:.1f}'.format(M)+', K = {:.1f}'.format(K))
addMsgAndPrint('  Scatter = {:.1f}'.format(scatter)+', Spread = {:.1f}'.format(spread))
addMsgAndPrint('  oblateness = {:.1f}'.format(oblateness)+'  OCD, 95% confidence = '+str(OCD))

        
if oblateness < testOb or n < minN or scatter > maxScatter or spread < minSpread:
    if n < minN:
        arcpy.AddWarning('****n = '+str(n)+' is less than testvalue of '+str(minN))
    if oblateness < testOb:
        arcpy.AddWarning('****Oblateness = {:.1f}'.format(oblateness)+' is less than testvalue of '+str(testOb))
    if scatter > maxScatter:
        arcpy.AddWarning('****scatter = {:.2f}'.format(scatter)+' is more than testvalue of '+str(maxScatter))    
    if spread < minSpread:
        arcpy.AddWarning('****spread = {:.2f}'.format(spread)+' is less than testvalue of '+str(minSpread))
    arcpy.AddWarning('****Data are insufficient to define a plane. Aborting.')
    forceExit()

## if necessary, project from demSR to xlbtSR
if demSR <> xlbtSR:
    addMsgAndPrint('Projecting mean location from WKID '+str(demSR.factoryCode)+' to WKID '+str(xlbtSR.factoryCode))
    newPoint = arcpy.PointGeometry(arcpy.Point(meanXYZd[0],meanXYZd[1]),demSR).projectAs(xlbtSR)

## add new feature to OrientationPoints
addMsgAndPrint('Adding result to '+orPts)
fields = ["SHAPE@XY","Type","Azimuth","Inclination","OrientationSourceID","LocationSourceID","OrientationConfidenceDegrees","Notes","Oblateness","Scatter","Spread"]
if hasSymbol:
    fields.append('Symbol')
cursor = arcpy.da.InsertCursor(orPts,fields)
xy = [newPoint[0].X,newPoint[0].Y]
if hasSymbol:
    cursor.insertRow([xy,calcType,round(strike,1),round(dip,1),dsID,dsID,OCD,iDs,round(oblateness,3),round(scatter,2),round(spread,1),calcSymbol])
else:
    cursor.insertRow([xy,calcType,round(strike,1),round(dip,1),dsID,dsID,OCD,iDs,round(oblateness,3),round(scatter,2),round(spread,1)])
arcpy.RefreshActiveView()

## write output to log file
orPtsDesc = arcpy.Describe(orPts)
logFileHost = arcpy.Describe(orPtsDesc.path)

if logFileHost.dataType == 'FeatureDataset':
    logPath = logFileHost.path
else:
    logPath = logFileHost.catalogPath
logfileName = os.path.join(logPath,os.path.basename(sys.argv[0])[:-3]+'-logfile.txt')
addMsgAndPrint('Writing results to '+logfileName)
logfile = open(logfileName,'a')
logfile.write(str(datetime.now())[:-7]+'\n')
logfile.write(iDs+'\n')
logfile.write('location = {:.0f}'.format(newPoint[0].X)+' {:.0f}'.format(newPoint[0].Y)+', srf = WKID '+str(xlbtSR.factoryCode)+'\n')
logfile.write('LocationSourceID = OrientationSourceID = '+dsID+'\n')
logfile.write('strike = {:.1f}'.format(strike)+', dip = {:.1f}'.format(dip)+',  OrientationConfidenceDegrees = '+str(OCD)+'\n')
logfile.write('n = '+str(n)+', Oblateness = {:.2f}'.format(oblateness))
logfile.write(', M = {:.1f}'.format(M)+', K = {:.1f}'.format(K))
logfile.write(', Scatter = {:.1f}'.format(scatter)+', Spread = {:.1f}'.format(spread)+'\n')
logfile.write('eigenvalues: '+str(eigenVals)+'\n')
logfile.write('eigenvectors: \n')
for x in eigenVects:
    logfile.write('  '+str(x)+'\n')
logfile.write('----\n')
logfile.close()

addMsgAndPrint("DONE!")
addMsgAndPrint('=====')
forceExit()




