#!/usr/bin/python3

import math
import numpy as np
import sys
from eglib import *

lat_imp =  34.3853
lon_imp = 132.4533

lat_ba  = 34.413462
lon_ba  = 132.951369
az      = 262.0

tmax = 44.9
dist = 45000.0
vo   = dist / tmax
print(vo)
g    = 9.8
el   = 0.0

elr = el * D2R
azr = az * D2R

vox = vo * math.cos(elr) * math.cos(azr)
voy = vo * math.sin(elr) * math.cos(azr)

fo = open('bmbpath.txt','w')

x0 = 0.0
y0 = 8974.0
t  = 0.0
dt = 0.10

x = x0
y = y0

while (y > 0.0):
   y = y0 + voy * t - 0.50 * g * t * t
   x = x0 + vox * t
   er = abs(x / 100.0)

   (nlat,nlon) = sphericalDisplaceLatitudeLongitude(lat_ba,lon_ba,az,er)

   fo.write('%10.3f %15.5f %15.5f %15.5f\n' % (t,nlon,nlat,y))

   t = t + dt
fo.close()
