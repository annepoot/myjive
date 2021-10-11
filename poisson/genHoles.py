import re
import os
import numpy as np
import pylab
import random
import math

nh = 2   # number of holes
r  = 2.0 # hole radius

dx = 10
dy = 10
le = 0.3
le = 0.06
le = 0.1
le = 0.5

seed = 311
seed = 11
seed = 1

random.seed ( seed )

with open ( 'mesh.geo' ,'w' ) as geo:
  geo.write ('Point(1) = {0,0,0,' + str(le) + '};\n')
  geo.write ('Point(2) = {' + str(dx) + ',0,0,' + str(le) + '};\n')
  geo.write ('Point(3) = {' + str(dx) + ',' + str(dy) + ',0,' + str(le) + '};\n')
  geo.write ('Point(4) = {0,' + str(dy) + ',0,' + str(le) + '};\n')
  geo.write ('Line(1) = {1,2};\n');
  geo.write ('Line(2) = {2,3};\n');
  geo.write ('Line(3) = {3,4};\n');
  geo.write ('Line(4) = {4,1};\n');

  counter = 4

  hxs = []
  hys = []
  
  for i in range(nh):
    hx = random.uniform ( r, dx-r )
    hy = random.uniform ( r, dy-r )
    if len(hxs) > 0:
      while any(math.sqrt((hx-x)*(hx-x)+(hy-y)*(hy-y)) < 2*r for x,y in zip(hxs,hys)):
        hx = random.uniform ( r, dx-r )
        hy = random.uniform ( r, dy-r )
    hxs.append ( hx )
    hys.append ( hy )
    geo.write ('Point(' + str(counter+1) + ') = {' + str(hx) + ',' + str(hy) + ',0,1.0};\n'   )
    geo.write ('Point(' + str(counter+2) + ') = {' + str(hx-r) + ',' + str(hy) + ',0,' + str(le) +'};\n' )
    geo.write ('Point(' + str(counter+3) + ') = {' + str(hx+r) + ',' + str(hy) + ',0,' + str(le) +'};\n' )
    counter = counter + 3

  for i in range(nh):
    fst = 5 + 3*i
    snd = 6 + 3*i
    trd = 7 + 3*i
    geo.write ('Circle(' + str(counter+1) + ') = {' + str(snd) + ',' + str(fst) + ',' + str(trd) + '};\n')
    geo.write ('Circle(' + str(counter+2) + ') = {' + str(trd) + ',' + str(fst) + ',' + str(snd) + '};\n')
    counter = counter + 2

  geo.write ('Line Loop(' + str(counter+1) + ') = {3,4,1,2};\n')

  counter = counter + 1
  bulkid  = counter

  for i in range(nh):
    fst = 5 + 3*nh + 2*i
    snd = 6 + 3*nh + 2*i
    geo.write ('Line Loop(' + str(counter+1) +') = {' + str(fst) + ',' + str(snd) + '};\n')
    counter = counter + 1

  geo.write ('Plane Surface(' + str(counter+1) + ') = {' + str(bulkid))

  counter2 = counter
  counter = counter + 1

  for i in range(nh):
    geo.write (',' + str(counter2))
    counter2 = counter2 - 1

  geo.write ('};\n')
  geo.write ('Physical Surface(' + str(counter+1) + ') = {' + str(counter) + '};\n')

os.system('gmsh -2 mesh.geo')
