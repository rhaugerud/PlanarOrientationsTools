# code to export a selection on OrientationPoints to a text file formatted
#  for Stereonet 10.0 (Allmendinger) input.

# This code is not subject to US Copyright

versionString = 'StereonetExport2.py, version of 25 May 2021'
import arcpy, os.path, time

debug = False
c = ', '

##########################
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
	
def forceExit():
    addMsgAndPrint('Forcing exit by raising ExecuteError')
    raise arcpy.ExecuteError

def testAndDelete(fc):
    if arcpy.Exists(fc):
        if debug: addMsgAndPrint('    deleting '+fc)
        try:
            arcpy.Delete_management(fc)
        except:
            addMsgAndPrint('failed to delete '+fc)
    return


##########################
op = sys.argv[1]  #  OrientationPoints layer
dl = sys.argv[2]  #  domain layer
typeField = sys.argv[3]
aziField = sys.argv[4]
incField = sys.argv[5]
outDir = sys.argv[6]
ofl = open(os.path.join(outDir,sys.argv[7]),'w')
comment3 = sys.argv[8]

addMsgAndPrint(versionString)

scratchName = 'xxx'+str(int(time.time()))[-5:]+'.gdb'

xxxgdb = os.path.join(outDir,scratchName)
testAndDelete(xxxgdb)
if not arcpy.Exists(xxxgdb):
    arcpy.CreateFileGDB_management(outDir,scratchName)
selDom = xxxgdb+'/selDom'
selOp1 = xxxgdb+'/selOp1'
selOp = xxxgdb+'/selOp'
prjOp = xxxgdb+'/prjOp'
if dl <> '#':
    arcpy.CopyFeatures_management(dl,selDom)
    arcpy.CopyFeatures_management(op,selOp1)
    arcpy.MakeFeatureLayer_management(selDom,'Domain')
    arcpy.MakeFeatureLayer_management(selOp1,'OrPts')
    arcpy.SelectLayerByLocation_management('OrPts','INTERSECT','Domain')
    arcpy.CopyFeatures_management('OrPts',selOp)
else:   # dl is blank
    arcpy.CopyFeatures_management(op, selOp)
spref = arcpy.SpatialReference(4326)
arcpy.Project_management(selOp,prjOp,spref,"")

mxd = arcpy.mapping.MapDocument("CURRENT")
opSource = op
for lyr in arcpy.mapping.ListLayers(mxd):
    if lyr.name == op:
        opSource = lyr.dataSource

#ofl.write('AZ\n')
ofl.write('# file written by '+versionString+'\n')
ofl.write('# data from '+opSource+'\n')
ofl.write('# '+comment3+'\n')

with arcpy.da.SearchCursor(prjOp,[aziField,incField,"SHAPE@X","SHAPE@Y",typeField]) as cursor:
    for row in cursor:
        if row[1] <> -9:  # skip if Inclination = NoData
            ofl.write(str(round(row[0],1)).rjust(5)+c+str(round(row[1],1)).rjust(4)+c+str(round(row[2],5))+c+str(round(row[3],5))+c+str(row[4])+'\n')
ofl.close()

del cursor
arcpy.Delete_management('Domain')
arcpy.Delete_management('OrPts')
arcpy.Compact_management(xxxgdb)
for fc in (selDom,selOp1,selOp,prjOp):
    testAndDelete(fc)
#testAndDelete(xxxgdb)

if arcpy.Exists(xxxgdb):
    addMsgAndPrint(scratchName+' is still here')

addMsgAndPrint('DONE!')
addMsgAndPrint('=====')
#forceExit()
