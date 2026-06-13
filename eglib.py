import math
import numpy as np

PI                = math.pi
HALFPI            = 0.50 * PI
TWOPI             = 2.0 * PI
D2R               = PI / 180.0
R2D               = 180.0 / PI
m2ft              = 100.0 / 30.48

earthFlat         = 1.0 / 298.257223563
earthRotationRate =         7.2921151467e-5
earthRadiusMeters =   6378137.0
earthRadiusKm     = earthRadiusMeters / 1000.0
earthGravAcc      = 9.80665

def oblateDisplaceLatitudeLongitude(olat,olon,az,er):

#  solution of the geodetic direct problem after t vincenty
#  modified rainsford method with helmert's elliptical terms
#  effective in any az and at any distance short of antipodal
#  computes baz 'backwards' azimuth (deg) but not returned
#
#  inputs:
#    olat   starting latitude point (deg)
#    olon   starting longitude point (deg)
#    az     angle WRT North to displace (deg)
#    er     angle along surface of earth to displace (km)
#
#  outputs:
#    dlat   displaced latitude (deg)
#    dlon   displaced lonigitude (deg)

   epsilon = 5.0e-12

   olat  = olat * D2R
   olon  = olon * D2R
   az    = az * D2R
   er    = er / earthRadiusKm

   r   = 1.0 - earthFlat
   tu  = r * math.sin(olat) / math.cos(olat)
   sf  = math.sin(az)
   cf  = math.cos(az)
   baz = 0.0

   if (cf != 0.0):
      baz = 2.0 * math.atan2(tu,cf)
   
   cu  = 1.0 / np.sqrt(1.0 + tu * tu)
   su  = tu * cu
   sa  = cu * sf
   c2a = 1.0 -sa * sa
   x   = np.sqrt((1.0 / r / r - 1.0) * c2a + 1.0) + 1.0
   x   = (x - 2.0) / x
   c   = 1.0 - x
   c   = (x * x / 4.0 + 1.0) / c
   d   = (0.375 * x * x - 1.0) * x
   tu  =  er / r / c
   y   = tu
   
   sy  = math.sin(y)
   cy  = math.cos(y)
   cz  = math.cos(baz + y)
   e   = 2.0 * cz * cz - 1.0
   c   = y
   x   = e * cy
   y   = e + e - 1.0
   y   = (((4.0 * sy * sy - 3.0) * y * cz * d / 6.0 + x) * d / 4.0 - cz) * sy * d + tu

   while (abs(y - c) > epsilon):
      sy  = math.sin(y)
      cy  = math.cos(y)
      cz  = math.cos(baz + y)
      e   = 2.0 * cz * cz - 1.0
      c   = y
      x   = e * cy
      y   = e * e - 1.0
      y   = (((4.0 * sy * sy - 3.0) * y * cz * d / 6.0 + x) * d / 4.0 - cz) * sy * d + tu

   baz = cu * cy * cf - su * sy
   c   = r * np.sqrt(sa * sa + baz * baz)
   d   = su * cy + cu * sy * cf

   dlat = math.atan2(d,c)

   c   = cu * cy - su * sy * cf
   d   = ((e * cy * c + cz) * sy * c + y) * sa
   x   = math.atan2(sy * sf,c)
   c   = ((-3.0 * c2a + 4.0) * earthFlat + 4.0) * c2a * earthFlat / 16.0
   d   = ((e * cy * c + cz) * sy * c + y) * sa

   dlon = olon + x - (1.0 - c) * d * earthFlat

   baz = math.atan2(sa,baz) + PI

   dlat = dlat * R2D
   dlon = dlon * R2D
   baz  = baz * R2D

#   if (dlon < 0.0):
#      dlon = dlon + 360.0

   return (dlat,dlon)

def oblateEarthAngleAzimuth(lat1,lon1,lat2,lon2):

#   solution of the geodetic inverse problem after t vincenty 
#   modified rainsfords method with helmerts elliptical terms
#   effective in any azimuth and at any distance
#
#  inputs:
#    lat1     latitude of point 1 (deg)
#    lon1     longitude of point 1 (deg)
#    lat2     latitude of point 2 (deg)
#    lon2     longitude of point 2 (deg)
#
#  outputs
#    az       azimuth from point 1 to point 2 (radians)
#    er       great circle path length (radians) (er * 6378.137 = km)

   er    = 0.0
   az    = 0.0

   eps   = 1.0e-8

   maximumTries = 128

   lat1 = lat1 * D2R
   lon1 = lon1 * D2R
   lat2 = lat2 * D2R
   lon2 = lon2 * D2R

   diflat = lat1 - lat2
   diflon = lon1 - lon2

   r       = 1.0 - earthFlat
   tu1     = r * math.sin(lat1) / math.cos(lat1)
   tu2     = r * math.sin(lat2) / math.cos(lat2)
   cu1     = 1.0 / np.sqrt(tu1 * tu1 + 1.0)
   su1     = cu1 * tu1
   cu2     = 1.0 / np.sqrt(tu2 * tu2 + 1.0)
   s       = cu1 * cu2
   baz     = s * tu2

   azimuth = baz * tu1

   x   = lon2 - lon1

   d   = 1.0e10

   counter = 0

   c2a = 0.0
   cy  = 0.0
   cz  = 0.0
   e   = 0.0
   sy  = 0.0
   y   = 0.0
   
   the_diff = abs(d-x)

   while (the_diff > eps):
      counter = counter + 1

      if (counter > maximumTries):
         print("too many tries")
         return 0.0,-1.0

      sx  = math.sin(x)
      cx  = math.cos(x)
      tu1 = cu2 * sx
      tu2 = baz - su1 * cu2 * cx
      sy  = np.sqrt(tu1 * tu1 + tu2 * tu2)
      cy  = s * cx + azimuth
      y   = math.atan2(sy,cy)

#  max(1.0e-10,sy) helps with case when point 1 = point 2
#  from dividing by zero

      sa  = s * sx / max(1.0e-10,sy)
      c2a = 1.0 - sa * sa
      cz  = azimuth + azimuth

      if (c2a > 0.0):
         cz = -cz / c2a + cy

      e   = 2.0 * cz * cz  - 1.0
      c   = ((-3.0 * c2a + 4.0) * earthFlat + 4.0) * c2a * earthFlat / 16.0
      d   = x
      x   = ((e * cy * c + cz) * sy * c + y) * sa
      x   = (1.0 - c) * x * earthFlat + lon2 - lon1

      the_diff = abs(d-x)

#   make forward azimuth (0-360)

   az = TWOPI + math.atan2(tu1,tu2)
   az = az % TWOPI

   x  = np.sqrt((1.0 / r  / r - 1.0)  * c2a + 1.0) + 1.0
   x  = (x - 2.0)  / x
   c  = 1.0 - x
   c  = (x  * x / 4.0 + 1.0) / c
   d  = (0.375 * x  * x - 1.0) * x
   x  = e * cy
   s  = 1.0 - e - e
   er = ((((sy * sy * 4.0 - 3.0) * s * cz * d / 6.0 - x) * d / 4.0 + cz ) * sy * d + y) * c * r

   return (az,er)

def sphericalDisplaceLatitudeLongitude(lat,lon,az,er):

#  displace a latitude, longitude on a spherical earth

#  inputs:
#    olat   starting latitude point (deg)
#    olon   starting longitude point (deg)
#    az     angle WRT North to displace (deg)
#    er     angle along surface of earth to displace (km)
#
#  outputs:
#    dlat   displaced latitude (deg)
#    dlon   displaced lonigitude (deg)

   azr  = az * D2R
   err  = er / earthRadiusKm
   latr = lat * D2R
   lonr = lon * D2R
   
   cosaz = math.cos(azr)
   sinaz = math.sin(azr)

   coser = math.cos(err)
   siner = math.sin(err)

   colat = math.cos(latr)
   silat = math.sin(latr)

   dlat = math.asin(cosaz * colat * siner + silat * coser)

   dlon = math.atan2(sinaz * siner * colat, coser - silat * math.sin(dlat))
   dlon = lonr + dlon

#   if (dlon < 0.0):
#      dlon = dlon + TWOPI

#   if (dlon > TWOPI):
#      dlon = dlon - TWOPI

   dlon = dlon * R2D
   dlat = dlat * R2D

   return (dlat,dlon)

def sphericalEarthAngleAzimuth(lat1,lon1,lat2,lon2):

#  determine angle az azimuth between two points on spherical earth
#  inputs:
#    lat1     latitude of point 1 (deg)
#    lon1     longitude of point 1 (deg)
#    lat2     latitude of point 2 (deg)
#    lon2     longitude of point 2 (deg)
#
#  outputs
#    az       azimuth from point 1 to point 2 (radians)
#    er       great circle path length (radians) (er * 6378.137 = km)

   lat1r = lat1 * D2R
   lon1r = lon1 * D2R
   lat2r = lat2 * D2R
   lon2r = lon2 * D2R

   theta1 = HALFPI - lat1r
   theta2 = HALFPI - lat2r

   ctheta1 = math.cos(theta1)
   stheta1 = math.sin(theta1)
   ctheta2 = math.cos(theta2)
   stheta2 = math.sin(theta2)

   cang = stheta1 * stheta2 * math.cos(lon1r-lon2r)  + ctheta1 * ctheta2
   cang = max(-1.0,min(1.0,cang))

   er   = math.acos(cang)
   az   = 0.0

   if (er > 0.0):
      sang = math.sqrt(1.0 - cang * cang)
      caz  = (ctheta2 - ctheta1 * cang) / (sang * stheta1)
      saz  = -stheta2 * math.sin(lon1r - lon2r) / sang
      az   = math.atan2(saz,caz)
      if (az < 0.0):
         az = az + TWOPI

   return (az,er)

