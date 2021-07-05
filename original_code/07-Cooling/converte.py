import os

"""
for n in xrange(241): # 241
   nome="FieldValues"+format(n,"04d")
   if os.path.isfile(nome+".ppm"):
      cmdexe="convert -size 700x700 "+nome+".ppm "+nome+".jpg"
      r=os.system(cmdexe)
      print "r, cmdexe = ",r,cmdexe
"""
np=0
for n in xrange(241): # 241
   nome="FieldValues"+format(n,"04d")
   nuovo="FieldValues"+format(np,"04d")
   if os.path.isfile(nome+".ppm"):
      cmdexe="mv "+nome+".jpg "+nuovo+".jpg"
      r=os.system(cmdexe)
      print "r, cmdexe = ",r,cmdexe   
      np=np+1
