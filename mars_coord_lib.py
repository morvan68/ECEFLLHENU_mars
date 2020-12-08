# -*- coding: utf-8 -*-
# @Author: wuzida
# @Date:   2018-12-28 15:20:54
# @Last Modified by:   wuzida
# @Last Modified time: 2018-12-29 16:49:16
import math

#set to mars2000, from
# Report of the IAU/IAG Working Group on cartographic coordinates and 
# rotational elements: 2006 P. Kenneth Seidelmann

a = 3396190.000
b = 3376200.000
# inverse flattening	169.8944472

f = (a - b) / a
e_sq = f * (2-f)

OMEGA_MARS = 7.088235959e-5  # Mars' rotation rate (rad/s)
# 24.6229 hours = 3600 * 24.6229
def geodetic_to_ecef(lat, lon, h):
    # (lat, lon) in degrees
    # h in meters
    lamb = math.radians(lat)
    phi = math.radians(lon)
    s = math.sin(lamb)
    N = a / math.sqrt(1 - e_sq * s * s)

    sin_lambda = math.sin(lamb)
    cos_lambda = math.cos(lamb)
    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)

    x = (h + N) * cos_lambda * cos_phi
    y = (h + N) * cos_lambda * sin_phi
    z = (h + (1 - e_sq) * N) * sin_lambda

    return x, y, z

def ecef_to_enu(x, y, z, lat0, lon0, h0):
    
    lamb = math.radians(lat0)
    phi = math.radians(lon0)
    s = math.sin(lamb)
    N = a / math.sqrt(1 - e_sq * s * s)

    sin_lambda = math.sin(lamb)
    cos_lambda = math.cos(lamb)
    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)

    x0 = (h0 + N) * cos_lambda * cos_phi
    y0 = (h0 + N) * cos_lambda * sin_phi
    z0 = (h0 + (1 - e_sq) * N) * sin_lambda

    xd = x - x0
    yd = y - y0
    zd = z - z0

    t = -cos_phi * xd -  sin_phi * yd

    xEast = -sin_phi * xd + cos_phi * yd
    yNorth = t * sin_lambda  + cos_lambda * zd
    zUp = cos_lambda * cos_phi * xd + cos_lambda * sin_phi * yd + sin_lambda * zd

    return xEast, yNorth, zUp

def enu_to_ecef(xEast, yNorth, zUp, lat0, lon0, h0):
    
    lamb = math.radians(lat0)
    phi = math.radians(lon0)
    s = math.sin(lamb)
    N = a / math.sqrt(1 - e_sq * s * s)

    sin_lambda = math.sin(lamb)
    cos_lambda = math.cos(lamb)
    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)

    x0 = (h0 + N) * cos_lambda * cos_phi
    y0 = (h0 + N) * cos_lambda * sin_phi
    z0 = (h0 + (1 - e_sq) * N) * sin_lambda

    t = cos_lambda * zUp - sin_lambda * yNorth

    zd = sin_lambda * zUp + cos_lambda * yNorth
    xd = cos_phi * t - sin_phi * xEast 
    yd = sin_phi * t + cos_phi * xEast

    x = xd + x0 
    y = yd + y0 
    z = zd + z0 

    return x, y, z

def ecef_to_geodetic( x, y, z, a = 3396190.0, b = 3376200.000 ):
   # Convert from ECEF cartesian coordinates to 
   # latitude, longitude and height.
   # lat lon in degrees, height in m
   # Mars2000 is default
   
    x2 = x ** 2
    y2 = y ** 2
    z2 = z ** 2

#    a = 6378137.0000    # earth radius in meters
#    b = 6356752.3142    # earth semiminor in meters 
    
    e = math.sqrt (1-(b/a)**2)
    b2 = b*b 
    e2 = e ** 2
    ep = e*(a/b)
    r = math.sqrt(x2+y2)
    r2 = r*r 
    E2 = a ** 2 - b ** 2
    F = 54*b2*z2
    G = r2 + (1-e2)*z2 - e2*E2
    c = (e2*e2*F*r2)/(G*G*G)
    s = ( 1 + c + math.sqrt(c*c + 2*c) )**(1/3) 
    P = F / (3 * (s+1/s+1)**2 * G*G)
    Q = math.sqrt(1+2*e2*e2*P) 
    ro = -(P*e2*r)/(1+Q) + math.sqrt((a*a/2)*(1+1/Q) - (P*(1-e2)*z2)/(Q*(1+Q)) - P*r2/2) 
    tmp = (r - e2*ro) ** 2 
    U = math.sqrt( tmp + z2 ) 
    V = math.sqrt( tmp + (1-e2)*z2 )
    zo = (b2*z)/(a*V) 

    height = U*( 1 - b2/(a*V) )
    
    lat = math.atan( (z + ep*ep*zo)/r )

    temp = math.atan(y/x)
    if x >=0 :
        longit = temp
    elif (x < 0) & (y >= 0):
        longit = math.pi + temp
    else :
        longit = temp - math.pi

    lat0 = math.degrees(lat)
    lon0 = math.degrees(longit)
    h0 = height 

    return lat0, lon0, h0

def geodetic_to_enu(lat, lon, h, lat_ref, lon_ref, h_ref):

    x, y, z = geodetic_to_ecef( lat, lon, h)
    
    return ecef_to_enu( x, y, z, lat_ref, lon_ref, h_ref)

def enu_to_geodetic(xEast, yNorth, zUp, lat_ref, lon_ref, h_ref):

    x,y,z = enu_to_ecef( xEast, yNorth, zUp, lat_ref, lon_ref, h_ref)

    return ecef_to_geodetic( x,y,z)
def mci_to_mcmf(Pos_ECI, T):
    """T in jd"""
    #https://space.stackexchange.com/questions/38807/transform-eci-to-ecef
    """
    In the IAU language, ECEF and ECI do not exist per se. But to the spirit of your question, yes, ITRS is the equivalent of ECEF, and GCRS is the equivalent of ECI. Note that the IAU frames move differently than the simplified models, so check with your customer what frames they need.  
    """
    ...
    # Reference Date : 12:00 UT 1 Jan 2000 (JD 2451545.0)
    
    #precession
    
    #nutation
    
    #rotation
    mcmf_pos = (Pos_ECI.x.value, Pos_ECI.y.value, Pos_ECI.z.value )
    
    C_ENU2ECEF = np.array([[ np.cos(lon), -np.sin(lon) , 0],
                           [ np.sin(lon), np.cos(lat),   0],
                           [     0      ,    0       ,   1]])
    
    #polar motion
    
    return mcmf_pos

def mcmf_to_mci( Pos_ECEF, T):
    ...
    mci_pos = (Pos_ECEF.x.value, Pos_ECEF.y.value, Pos_ECEF.z.value )
    
    return mci_pos

def xyz2llh(x,y,z):
    '''
    alternative from web, currently not used.
    currently still earth not mars
    https://gis.stackexchange.com/questions/265909/converting-from-ecef-to-geodetic-coordinates
    Function to convert xyz ECEF to llh
    convert cartesian coordinate into geographic coordinate
    ellipsoid definition: WGS84
      a= 6,378,137m
      f= 1/298.257

    Input
      x: coordinate X meters
      y: coordinate y meters
      z: coordinate z meters
    Output
      lat: latitude rad
      lon: longitude rad
      h: height meters
    '''
    # --- WGS84 constants
    a = 6378137.0
    f = 1.0 / 298.257223563
    # --- derived constants
    b = a - f*a
    e = math.sqrt(math.pow(a,2.0)-math.pow(b,2.0))/a
    clambda = math.atan2(y,x)
    p = math.sqrt(pow(x,2.0)+pow(y,2))
    h_old = 0.0
    # first guess with h=0 meters
    theta = math.atan2(z,p*(1.0-math.pow(e,2.0)))
    cs = math.cos(theta)
    sn = math.sin(theta)
    N = math.pow(a,2.0)/math.sqrt(math.pow(a*cs,2.0)+math.pow(b*sn,2.0))
    h = p/cs - N
    while abs(h-h_old) > 1.0e-6:
        h_old = h
        theta = math.atan2(z,p*(1.0-math.pow(e,2.0)*N/(N+h)))
        cs = math.cos(theta)
        sn = math.sin(theta)
        N = math.pow(a,2.0)/math.sqrt(math.pow(a*cs,2.0)+math.pow(b*sn,2.0))
        h = p/cs - N
    llh = {'lon':clambda, 'lat':theta, 'height': h}
    return llh
