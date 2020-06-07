# script to ensure that necessary fields are present in
# point feature class used for BeddingTracesToS&D.py output.
# Ralph Haugerud, U.S. Geological Survey
# Code not subject to U.S. copyright

import arcpy, sys

versionString = 'PrepOrientationPoints.py, version of 1 June 2020'

defaultLength = 254
IDLength = 50
memoLength = 3000
booleanLength = 1

neededFields = [['Type','String',defaultLength],
                ['Azimuth','Single',''],
                ['Inclination','Single',''],
                ['OrientationConfidenceDegrees','Single'],
                ['OrientationSourceID','String',IDLength],
                ['Notes','String',defaultLength],
                ['Oblateness','Single',''],
                ['Scatter','Single',''],
                ['Spread','Single','']]

# get input feature class

op = sys.argv[1]

# cycle through fields
fields = arcpy.ListFields(op)
fieldNames = []
for f in fields:
    fieldNames.append(f.name.encode("ascii"))
#try:
for ff in neededFields:
    if not ff[0] in fieldNames:
        print 'adding field '+ff[0]
        arcpy.AddField_management(op,ff[0],ff[1],'','',ff[2])
    elif ff[1] <> fields[fieldNames.index(ff[0])].type:
        print 'field '+ff[0]+', type mismatch, '+ff[1]+' <> '+fields[fields.index(ff[0])].type
#except:
#    print 'Is there an open edit session?  If so, please save edits and stop editing!'
#    print 'Then try this again.'
    
