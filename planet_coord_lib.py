import math

from astropy.time import Time

#mars2000, from
# Report of the IAU/IAG Working Group on cartographic coordinates and 
# rotational elements: 2006 P. Kenneth Seidelmann

planet = 'Mars'
print('planet lib, ', planet)

if planet == 'Mars':
    OMEGA_PLANET = 7.088235959e-5  # rotation rate (rad/s)
    # 24.6229 hours
elif planet == 'Earth':
    OMEGA_PLANET = 7.2921158553e-5
else:
    print('planet not found')
    exit(0)
def geodetic_to_ecef(lat, lon, h, a = 3396190.000, b = 3376200.000 ):
    # (lat, lon) in radians
    # h in meters
    
    f = (a - b) / a #flattening
    e_sq = f * (2-f)

    s = math.sin(lat)
    N = a / math.sqrt(1 - e_sq * s * s)

    x = (h + N) * math.cos(lat) * math.cos(lon)
    y = (h + N) * math.cos(lat) * math.sin(lon)
    z = (h + (1 - e_sq) * N) * math.sin(lat)

    return x, y, z

def ecef_to_geodetic( x, y, z, a = 3396190.0, b = 3376200.000 ):
   # Convert from ECEF cartesian coordinates to 
   # return latitude, longitude and height.
   # lat lon in radians, height in m
   # Mars2000 is default a, b
   
    x2 = x ** 2
    y2 = y ** 2
    z2 = z ** 2

    E2 = a ** 2 - b ** 2

    e = math.sqrt (1-(b/a)**2)
    b2 = b*b 
    e2 = e ** 2
    ep = e*(a/b)
    r = math.sqrt(x2+y2)
    r2 = r*r 
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

    return lat, longit, height

def ecef_to_enu(x, y, z, lat0, lon0, h0, a = 3396190.0, b = 3376200.0 ):
    """lat lon of reference in radians"""
    lamb = lat0
    phi = lon0
    
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

def enu_to_ecef(xEast, yNorth, zUp, lat0, lon0, h0, a = 3396190.0, b = 3376200.000 ):
    """lat lon of reference in radians"""

    lamb = lat0
    phi = lon0
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

def geodetic_to_enu(lat, lon, h, lat_ref, lon_ref, h_ref):
    """lat lon of reference in radians"""

    x, y, z = geodetic_to_ecef( lat, lon, h)
    
    return ecef_to_enu( x, y, z, lat_ref, lon_ref, h_ref)

def enu_to_geodetic(xEast, yNorth, zUp, lat_ref, lon_ref, h_ref):
    """lat lon of reference in radians"""

    x,y,z = enu_to_ecef( xEast, yNorth, zUp, lat_ref, lon_ref, h_ref)

    return ecef_to_geodetic( x,y,z)
def pcpf_to_pci( Pos_ECEF, T):
    """
    T is astropy time
    PosECEF is np array
    
    Reference Date : 12:00 UT 1 Jan 2000 (JD 2451545.0)
    """
    T0 = Time(2451545.0, format='jd', scale='utc')
    dt = T - T0
    dt.format = 'sec'

    xx = Pos_ECEF[0]
    yy = Pos_ECEF[1]
    zz = Pos_ECEF[2]

    #precession

    #nutation

    #rotation
    agl = dt.value  * OMEGA_PLANET #radians
    agl = agl % (2 * math.pi)
#    print('pcpf to pci angle, ', agl)
#    C_ENU2ECEF = np.array([[ np.cos(lon), -np.sin(lon) , 0],
#                          [ np.sin(lon), np.cos(lat),   0],
#                           [     0      ,    0       ,   1]])
    mci_pos = ( xx * math.cos(-agl) - yy * math.sin(-agl),
                xx * math.sin(-agl) + yy * math.cos(-agl),
                zz )

    #polar motion

    return mci_pos

def pci_to_pcpf(Pos_ECI, T):
    """
    T is astropy time
    Pos_ECI is np array

    https://space.stackexchange.com/questions/38807/transform-eci-to-ecef
    In the IAU language, ECEF and ECI do not exist per se. But to the spirit of your question, yes,
    ITRS is the equivalent of ECEF, and GCRS is the equivalent of ECI. Note that the IAU frames move
    differently than the simplified models, so check with your customer what frames they need.
    """
    # Reference Date : 12:00 UT 1 Jan 2000 (JD 2451545.0)
    T0 = Time(2451545.0, format='jd', scale='utc')
#    print('t, ',T, T0)
    dt = T - T0
    dt.format = 'sec'
#    print('dtif astropy T, ',dt)
#    print('dtif astropy Tv, ',dt.value)

    xx = Pos_ECI[0]
    yy = Pos_ECI[1]
    zz = Pos_ECI[2]

    #precession
    
    #nutation
    
    #rotation
    agl = dt.value * OMEGA_PLANET #radians
    agl = agl % (2 * math.pi)
#    print('pci to pcpf angle rad, deg, ', agl, math.degrees(agl) )

#    C_ENU2ECEF = np.array([[ np.cos(lon), -np.sin(lon) , 0],
#                          [ np.sin(lon), np.cos(lat),   0],
#                           [     0      ,    0       ,   1]])
    mcmf_pos = (xx * math.cos(agl) - yy * math.sin(agl),
                xx * math.sin(agl) + yy * math.cos(agl),
                zz )
    
    #polar motion
    
    return mcmf_pos

def xyz2llh( x,y,z, a = 3396190.0, b = 3376200.000 ):
    '''
    alternative from web, currently not used.
    https://gis.stackexchange.com/questions/265909/converting-from-ecef-to-geodetic-coordinates
    Function to convert xyz ECEF to lat lon h
    convert cartesian coordinate into geographic coordinate
    Input
      x: coordinate X meters
      y: coordinate y meters
      z: coordinate z meters
    Output
      lat: latitude rad
      lon: longitude rad
      h: height meters
    '''
    f = (a - b) / a
    
    # --- derived constants
    e = math.sqrt(math.pow(a,2.0) - math.pow(b,2.0)) / a
    clambda = math.atan2(y,x)
    p = math.sqrt(pow(x,2.0) + pow(y,2))
    h_old = 0.0
    # first guess with h=0 meters
    theta = math.atan2(z,p*(1.0 - math.pow(e,2.0)))
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
    llh = [theta, clambda, h]
    return llh
if __name__ == '__main__':
    
    x = 1000
    y = 1000
    z = 1000
    
    scale = 3
    x, y, z = scale*2566618.9, scale*1452025.2, scale*1679729.6
    a = 3396190.0
    b = 3376200.0

#    #wgs84
#    a = 6378137.0  # Semi-major axis
#    f = 1 / 298.257223563  # Flattening
#    b = a * (1 - f)  # Semi-minor axis
    
    print( 'a,b, ', a, b)
    
    aa = ecef_to_geodetic( x, y, z, a=a, b=b)
    print( 'ecef, ', math.degrees(aa[0]), math.degrees(aa[1]), aa[2] )
    
    bb = xyz2llh( x,y,z, a=a, b=b)
    print( 'xyz, ', math.degrees(bb[0]), math.degrees(bb[1]), bb[2] )