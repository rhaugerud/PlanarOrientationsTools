# code to export a selection on OrientationPoints to a text file formatted
#  for Stereonet 10.0 (Allmendinger) input.

# This code is not subject to US Copyright

versionString = 'StereonetExport.py, version of 11 February 2020'
import arcpy, os.path

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
        arcpy.Delete_management(fc)
    return


##########################
op = sys.argv[1]
typeField = sys.argv[2]
aziField = sys.argv[3]
incField = sys.argv[4]
ofl = open(os.path.join(sys.argv[5],sys.argv[6]),'w')
comment3 = sys.argv[7]

addMsgAndPrint(versionString)

xxxgdb = os.path.join(sys.argv[5],'xxxxScratch.gdb')
testAndDelete(xxxgdb)
arcpy.CreateFileGDB_management(sys.argv[5],'xxxxScratch.gdb')
selOp = xxxgdb+'/selOp'
prjOp = xxxgdb+'/prjOp'
arcpy.CopyFeatures_management(op,selOp)
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

testAndDelete(xxxgdb)

addMsgAndPrint('DONE!')
addMsgAndPrint('=====')
#forceExit()
