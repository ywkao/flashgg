associated = ['WH','ZH','VBF']
thMap = {'bbH':'testBBH','THQ':'testTHQ','THW':'testTHW'}

from os import listdir
from os import system
for filename in listdir('./'):
  for thName in thMap:
    if thName in filename:
      print "renaming %s.root as %s.root"%(thName,thMap[thName])
      system("mv %s %s.root"%(filename,thMap[thName]))
  for assoc in associated:
    if ('%sToGG'%assoc in filename or '%sHToGG'%assoc in filename) and 'VH2HQQ' in filename:
      newfilename = filename.replace('VH2HQQ','%s2HQQ'%assoc)
      newfilename = newfilename.replace('VBF2HQQ','VBF')
      print "renaming",filename,"as",newfilename
      system("mv %s %s"%(filename,newfilename))
