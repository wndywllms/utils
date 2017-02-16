import numpy as np
import math

def angular_separation(ra1,dec1,RA,DEC):
    separation = np.sqrt((np.cos(DEC*np.pi/180.))**2.*(ra1-RA)**2.  + (dec1-DEC)**2.)
    return separation

def angsep(ra1deg, dec1deg, ra2deg, dec2deg):
    """Returns angular separation between two coordinates (all in degrees)"""
    #import math
    import astropy.coordinates as ac
    
    
    #ra1rad=ra1deg*math.pi/180.0
    #dec1rad=dec1deg*math.pi/180.0
    #ra2rad=ra2deg*math.pi/180.0
    #dec2rad=dec2deg*math.pi/180.0
    
    ra1deg = ac.Angle('{deg}d'.format(deg=str(ra1deg)))
    dec1deg = ac.Angle('{deg}d'.format(deg=str(dec1deg)))
    ra2deg = ac.Angle('{deg}d'.format(deg=str(ra2deg)))
    dec2deg = ac.Angle('{deg}d'.format(deg=str(dec2deg)))
    sep = ac.angle_utilities.angular_separation(ra1deg, dec1deg, ra2deg, dec2deg)
    sep = ac.Angle(sep)
    
    ## calculate scalar product for determination
    ## of angular separation
    #x=math.cos(ra1rad)*math.cos(dec1rad)*math.cos(ra2rad)*math.cos(dec2rad)
    #y=math.sin(ra1rad)*math.cos(dec1rad)*math.sin(ra2rad)*math.cos(dec2rad)
    #z=math.sin(dec1rad)*math.sin(dec2rad)
    
    #rad=math.acos(x+y+z)
        
    ## use Pythargoras approximation if rad < 1 arcsec
    #if rad<0.000004848:
        #rad=math.sqrt((math.cos(dec1rad)*(ra1rad-ra2rad))**2+(dec1rad-dec2rad)**2)
        
    ## Angular separation
    #deg=rad*180/math.pi
    #return deg
    return sep.degree

def ra_to_degrees(ra_str, delim=' '):
    '''
    converts array of strings or single string ra values to decimal degrees
    '''
    # string or array
    if isinstance(ra_str, str):
        t = ra_str.split(delim)
        ra_deg = 15.*(float(t[0]) + float(t[1])/60. + float(t[2])/3600.)
        return ra_deg
    else:
        ra_deg = np.zeros(len(ra_str))
        for i,ra_s in enumerate(ra_str):
            t = ra_s.split(delim)
            ra_deg[i] = 15.*(float(t[0]) + float(t[1])/60. + float(t[2])/3600.)
        return ra_deg

def dec_to_degrees(dec_str, delim=' '):
    '''
    converts array of strings or single string dec values to decimal degrees
    '''
    # string or array
    if isinstance(dec_str,str):
        t = dec_str.split(delim)
        dec_deg = (float(t[0]) + float(t[1])/60. + float(t[2])/3600.)
        return dec_deg
    else:
        dec_deg = np.zeros(len(dec_str))
        for i,dec_s in enumerate(dec_str):
            t = dec_s.split(delim)
            dec_deg[i] = (float(t[0]) + float(t[1])/60. + float(t[2])/3600.)
        return dec_deg        

def ra_to_str(dra, ndec=2,delim=':'):
    '''
    converts a single decimal degrees ra to hh:mm:ss.s
    '''
    dra = dra/15.
    dd = math.floor(dra)
    dfrac = dra - dd
    dmins = dfrac*60.
    dm = math.floor(dmins)
    dsec = (dmins-dm)*60.
    if round(dsec, ndec) == 60.00:
        dsec = 0.
        dm += 1
    if dm == 60.:
        dm = 0.
        dd += 1
    sra = '%02d%s%02d%s%05.2f' %(dd,delim,dm,delim,dsec)  
    return sra

def dec_to_str(ddec,ndec=1,delim=':'):
    '''
    converts a single decimal degrees dec to dd:mm:ss.s
    '''
    if ddec  >=0 :
        dd = math.floor(ddec)
        dfrac = ddec - dd
        dmins = dfrac*60.
        dm = math.floor(dmins)
        dsec = (dmins-dm)*60.
        if round(dsec, ndec) == 60.0:
            dsec = 0.
            dm += 1
        if dm == 60.:
            dm = 0.
            dd += 1
        sdec = '%02d%s%02d%s%04.1f' %(dd,delim,dm,delim,dsec)
    else:
        ddec = abs(ddec)
        dd = math.floor(ddec)
        dfrac = ddec - dd
        dmins = dfrac*60.
        dm = math.floor(dmins)
        dsec = (dmins-dm)*60.
        if round(dsec, ndec) == 60.0:
            dsec = 0.
            dm += 1
        if dm == 60.:
            dm = 0.
            dd += 1
        sdec = '-%02d%s%02d%s%04.1f' %(dd,delim,dm,delim,dsec)
    return sdec
    
def match_celestial (ra1, dec1, ra2, dec2, radius, verbose=1):
    '''
    input: ra1,dec1 - sky coordinates (in degrees) of array1
           ra2,dec2 - sky coordinates (in degrees) of array2
           radius - radius (in degrees) within which to provide matches
    returns: ind1 - matching indicies of array 1
             ind2 - matching indicies of array 2
             dist - matching sky distances in arcsec
    '''
    if verbose > 1:
        print 'sky matching with radius %.3f arcsec' %(radius*3600.)
    # radius in degrees
    # coordinates in degrees
    # make everything arcsec
    #ra1 *= 3600.
    #ra2 *= 3600.
    #dec1 *= 3600.
    #dec2 *= 3600.
    #radius *= 3600.
    ind1 = np.array([],dtype=int)
    ind2 = np.array([],dtype=int)
    dist = np.array([])
    # for each source in array 1, find all matching sources in array 2
    for i in range(len(ra1)):
        cd = np.cos(dec1[i]*np.pi/180.)**2.  # for dec in degrees
        sep = np.sqrt( (cd)*(ra1[i] - ra2)**2. + (dec1[i] - dec2)**2. )
        temp_ind2 = np.where(sep<radius)[0]
        #print temp_ind2
        if len(temp_ind2) > 0:
            ind2 = np.concatenate((ind2,temp_ind2))
            # store the index 1 time the number of matches
            ind1 = np.concatenate((ind1,i*np.ones(len(temp_ind2))))
            dist = np.concatenate((dist,sep[temp_ind2]/3600.))
    dist *= 3600.
    return ind1, ind2, dist


def coordinates_to_name(ra,dec,mode='truncate', nsigra=3, nsigdec=2, prefix='J'):
    
    drh = ra/15.  # decimal hours
    rh = int(math.floor(drh))   # integer decimal hrs
    dfrh = drh - rh  # fractional decimal hrs
    drm = dfrh*60.   # decimal mins
    rm = int(math.floor(drm))  # integer decimal mins
    dfrm = drm-rm    # fractional decimal mins
    drs = dfrm*60.   # decimal secs
    rs = int(math.floor(drs))  # integer decimal secs
    dfrs = drs-rs    # fractional decimal secs
    frs = dfrs*10**nsigra
    if mode == 'truncate':
        frs = int(math.floor(frs))
    else:
        frs = int(np.round(frs))
    
    
    if dec/abs(dec) > 0:
        decsign = '+'
    else:
        decsign = '-'
        
    dec = abs(dec)
        
    ddd = dec  # decimal degrees
    dd = int(math.floor(ddd))   # integer decimal hrs
    dfdd = ddd - dd  # fractional decimal hrs
    ddm = dfdd*60.   # decimal mins
    dm = int(math.floor(ddm))  # integer decimal mins
    dfdm = ddm-dm    # fractional decimal mins
    dds = dfdm*60.   # decimal secs
    ds = int(math.floor(dds))  # integer decimal secs
    dfds = dds-ds    # fractional decimal secs
    fds = dfds*10**nsigdec
    if mode == 'truncate':
        fds = int(math.floor(fds))
    else:
        fds = int(np.round(fds))
    
    
    
    
        
    if  nsigdec>0:
        sname = '{prefix}{rh:02d}{rm:02d}{rs:02d}.{frs:0{nsigra}d}{dsign:1}{dd:02d}{dm:02d}{ds:02d}.{fds:0{nsigdec}d}'.format(prefix=prefix, rh=rh, rm=rm, rs=rs, frs=frs, nsigra=nsigra, dsign=decsign, dd=dd, dm=dm, ds=ds, fds=fds, nsigdec=nsigdec)
    else:
        sname = '{prefix}{rh:02d}{rm:02d}{rs:02d}.{frs:0{nsigra}d}{dsign:1}{dd:02d}{dm:02d}{ds:02d}'.format(prefix=prefix, rh=rh, rm=rm, rs=rs, frs=frs, nsigra=nsigra, dsign=decsign, dd=dd, dm=dm, ds=ds, fds=fds, nsigdec=nsigdec)

    #if mode == truncate:
    
    return sname
    