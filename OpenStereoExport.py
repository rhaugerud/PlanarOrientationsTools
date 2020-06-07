# code to export a selection on OrientationPoints to a text file formatted
#  for Openstereo input.
#
# Ralph Haugerud, U.S. Geological Survey
# Code not subject to U.S. copyright

versionString = 'OpenStereoExport.py, version of 11 February 2020'
import arcpy, os.path

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


##########################
op = sys.argv[1]
aziField = sys.argv[2]
incField = sys.argv[3]
ofl = open(os.path.join(sys.argv[4],sys.argv[5]),'w')
comment3 = sys.argv[6]

mxd = arcpy.mapping.MapDocument("CURRENT")
for lyr in arcpy.mapping.ListLayers(mxd):
    if lyr.name == op:
        opSource = lyr.dataSource

addMsgAndPrint(opSource)

ofl.write('# file written by '+versionString+'\n')
ofl.write('# data from '+op+'\n')
ofl.write('# '+comment3+'\n')

with arcpy.da.SearchCursor(op,[aziField,incField]) as cursor:
    for row in cursor:
        ofl.write(str(row[0])+' '+str(row[1])+'\n')

ofl.close()

addMsgAndPrint('DONE!')
addMsgAndPrint('=====')

