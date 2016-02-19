
import yoda
import glob
import os
import pprint

AllYodaFiles = glob.glob("*.yoda")

for file in AllYodaFiles:
  histosDict = yoda.readYODA( file )
  newDict = { key: value for key, value in histosDict.items()
                         if (key != '/_XSEC') and (key != '/_EVTCOUNT') }
  os.rename( file, file+".old" )
  yoda.writeYODA( newDict, file )
