#!/usr/bin/python3

import math
import numpy as np
import re
import sys
from eglib import *

p = dict()

fi = open('man.inp','r')

for aLine in fi:
   aLine   = aLine.strip()
   parts   = aLine.split()
   val     = float(parts[0])
   var     = parts[1]
   p[var] = val

fi.close()

#  math came from "The math of saving the Enola Gay #SoME3" Dr. Jorge S. Diaz (YouTube)

#  8 -> 9 miles (13 -> 14 km) was considered safe distance based on Trinity test

rsafe   = 13.0    # km

# params

oldway = 0
if (oldway):
   h0    = 8974.87    # altitude (m)
   vac   = 133.73     # speed (m/s)

   vshk  = 336.0      # from LAMB
   
   tfall = 44.4       # fall time of weapon (agrees with runge kutta (RK) sim)
   
   theta = 158.0      # angle to get a/c on correct departure path
   
   er    = 6.2        # approx dist from impact point bomb released (km, from RK sim)
   
   shkarv = 55.0      # approx time after detonation shock arrives (from maneuv.py)
   
#  what to assume is the maximum speed attained during maneuver
   
   vmax  = vac
   vmax  = 350.0 * 0.44704   # mph->m/s
   vmax  = 150.00
else:
   vac     = p['vac']
   vmax    = p['vmax']
   shkarv  = p['tshk']
   h0      = p['hAircraft']
   hdet    = p['hDetonation']
   R       = p['R']
   tfall   = p['tfall']
   theta   = p['bestang']
   appd    = p['approx']

hd0   = 262.0      # heading (deg)

er    = 5.5        # approx dist from impact point bomb released (km, from RK sim)
   
bank  = 60.0       # angle a/c banked during right hand turn

cmtr2mi = m2ft / 5280.0

#  hiroshima

lat_h =  34.3853
lon_h = 132.4533

#  take flight log azimuth and get azimuth from hiroshima to the plane
#  to compute the bomb away point

head = (hd0 + 180) % 360

(balat,balon) = oblateDisplaceLatitudeLongitude(lat_h,lon_h,head,er)

fo = open('manu.out','w')

print('vac %10.3f vmax %10.3f vmax/vac %10.3f' % (vac,vmax,vmax/vac))

#  compute turn radius (from youtube video)

th = theta * D2R
R  =  vmax * vmax /(9.8 * math.tan(bank * D2R))

tab = R * th / vmax

print('Radius of maneuver: %10.3f (m) time avail: %10.3f (sec)' % (R,tab))

#  determine the center lat,lon of maneuver

er  = R / 1000.0

(clat,clon) = oblateDisplaceLatitudeLongitude(balat,balon,0.0,er)

fo.write('%15.5f %15.5f %2d %2d\n' % (lon_h,lat_h,1,0))
fo.write('%15.5f %15.5f %2d %2d\n' % (balon,balat,1,1))

fo.write('\n')
fo.write('%15.5f %15.5f %2d %2d\n' % (clon,clat,2,2))
fo.write('\n')

#  generate points along the maneuver

az0 =  0.0 + 180.0
az1 = theta + 180.0
daz =  5.0

az = az0

while (az <= az1):
   (nlat,nlon) = oblateDisplaceLatitudeLongitude(clat,clon,az,er)

   fo.write('%15.5f %15.5f %2d %2d\n' % (nlon,nlat,3,0))

   if (az > az0):
      (faz,xer) = oblateEarthAngleAzimuth(olat,olon,nlat,nlon)
      faz = faz * R2D

   olat = nlat
   olon = nlon

   az = az + daz

fo.write('\n')

#  how much of the bomb fall time did we use up

dt = tfall - tab

print('time left after maneuver %10.3f' % (dt))

#  fly away at this azimuth for time remaining

er = vac * dt / 1000.0

(btlat,btlon) = oblateDisplaceLatitudeLongitude(nlat,nlon,faz,er)

fo.write('%15.5f %15.5f %2d %2d\n' % (nlon,nlat,4,0))
fo.write('%15.5f %15.5f %2d %2d\n' % (btlon,btlat,4,1))
fo.write('\n')

#  how much time is left before shock arrives

trem = (tfall + shkarv) - (tab + tfall)

print('trem %10.3f' % (trem))

er   = trem * vmax / 1000.0

fo.write('%15.5f %15.5f %2d %2d\n' % (btlon,btlat,5,0))
(salat,salon) = oblateDisplaceLatitudeLongitude(btlat,btlon,faz,er)
fo.write('%15.5f %15.5f %2d %2d\n' % (salon,salat,5,1))
fo.write('\n')

#  add safe range circle

az0 =   0.0
az1 = 360.0
daz =   5.0

er = rsafe
az = az0

while (az <= az1):
   (slat,slon) = oblateDisplaceLatitudeLongitude(lat_h,lon_h,az,er)
   fo.write('%15.5f %15.5f %2d %2d\n' % (slon,slat,6,1))
   az = az + daz

fo.close()

#  final distance

(sa2braz,sa2brer) = oblateEarthAngleAzimuth(lat_h,lon_h,salat,salon)

print('burst to shock arrival %10.3f (m) ' % (sa2brer * earthRadiusMeters),end='')
print('%10.3f (miles)' % (sa2brer * earthRadiusMeters * cmtr2mi))
