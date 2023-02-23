import numpy as np
import casacore.tables as pt
import casacore.quanta as qa
import casacore.measures as pm


from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
import astropy.units as u


def get_optimum_size(size):
    """
    Gets the nearest optimum image size

    Taken from the casa source code (cleanhelper.py)

    Parameters
    ----------
    size : int
        Target image size in pixels

    Returns
    -------
    optimum_size : int
        Optimum image size nearest to target size

    """

    def prime_factors(n, douniq=True):
        """ Return the prime factors of the given number. """
        factors = []
        lastresult = n
        sqlast=int(np.sqrt(n))+1
        if n == 1:
            return [1]
        c=2
        while 1:
            if (lastresult == 1) or (c > sqlast):
                break
            sqlast=int(np.sqrt(lastresult))+1
            while 1:
                if(c > sqlast):
                    c=lastresult
                    break
                if lastresult % c == 0:
                    break
                c += 1
            factors.append(c)
            lastresult /= c
        if (factors==[]): factors=[n]
        return  np.unique(factors).tolist() if douniq else factors

    n = int(size)
    if (n%2 != 0):
        n+=1
    fac=prime_factors(n, False)
    for k in range(len(fac)):
        if (fac[k] > 7):
            val=fac[k]
            while (np.max(prime_factors(val)) > 7):
                val +=1
            fac[k]=val
    newlarge=np.product(fac)
    for k in range(n, newlarge, 2):
        if ((np.max(prime_factors(k)) < 8)):
            return k
    return newlarge




def time_smearing(resolution,delta_Theta):
    #http://www.cv.nrao.edu/course/astr534/Interferometers1.html
    # output and input units are equal
    # Condition that delta_Theta * delta_T << 1.37E4 * Resolution
    # where delta_Theta is the angular radius of the image and delta_T is the time smoothing.

    # Given resolution and field of view (delta_theta) this returns a condition on the time smearing at the correlator
    delta_T = 1.37E4*resolution/delta_Theta

    print('Time averaging should be less than %s'%delta_T)


    # Time average smearing amplitude loss in http://adsabs.harvard.edu/full/1989ASPC....6..247B
    Reduction = 1-1.22E-9*(delta_Theta/resolution)**2.0 * delta_T**2.0

    print('At radius %s a source will only have %s percent of its flux if data smoothed in time to %s'%(delta_Theta,Reduction,delta_T))
    
    return delta_T

def time_smearing2(delta_T,resolution,delta_Theta):
    #http://www.cv.nrao.edu/course/astr534/Interferometers1.html
    # output and input units are equal
    # Condition that delta_Theta * delta_T << 1.37E4 * Resolution
    # where delta_Theta is the angular radius of the image and delta_T is the time smoothing.
    # Time average smearing amplitude loss in http://adsabs.harvard.edu/full/1989ASPC....6..247B

    # Given resolution and field of view (delta_theta) this returns a condition on the time smearing at the correlator
    Reduction = 1-1.22E-9*(delta_Theta/resolution)**2.0 * delta_T**2.0

    return Reduction
    
def bandwidth_smearing(freq,resolution,delta_Theta):

    #http://www.cv.nrao.edu/course/astr534/Interferometers1.html
    # Output and input units are equal
    # Condition that delta_Theta * delta_freq << Freq * resolution
    # where delta_Theta is the offset from the pointing centre and delta_freq is the bandwidth smoothing.

    # Bandwidth smearing amplitude loss in http://adsabs.harvard.edu/full/1989ASPC....6..247B
    
    # Given resolution, freq and offset this gives the condition for the delta_freq

    import scipy.special

    delta_freq = freq*resolution/delta_Theta

    print('Bandwidth averaging should be much less than %s'%delta_freq)

    
    beta = (delta_freq/freq) * (delta_Theta/resolution)
    gamma = 2*(np.log(2)**0.5)
    Reduction = ((np.pi**0.5)/(gamma * beta)) * (scipy.special.erf(beta*gamma/2.0))

    print('At radius %s a source will only have %s percent of its flux if data smoothed in freq to %s'%(delta_Theta,Reduction,delta_freq))

    return delta_freq

def bandwidth_smearing2(delta_freq,freq,resolution,delta_Theta):
    '''
    #http://www.cv.nrao.edu/course/astr534/Interferometers1.html
    # Output and input units are equal
    # Condition that delta_Theta * delta_freq << Freq * resolution
    # where delta_Theta is the offset from the pointing centre and delta_freq is the bandwidth smoothing.

    # Given resolution, freq and offset this gives the condition for the delta_freq
    '''
    import scipy.special


    beta = (delta_freq/freq) * (delta_Theta/resolution)
    gamma = 2*(np.log(2)**0.5)
    Reduction = ((np.pi**0.5)/(gamma * beta)) * (scipy.special.erf(beta*gamma/2.0))

    return Reduction


def check_flagged(ms):
    ''' check_flagged(ms) check the flag fraction of a measurement set
    from ddf-pipeline
    '''
    
    import pyrap.tables as pt
    t = pt.table(ms, readonly=True, ack=False)
    tc = t.getcol('FLAG').flatten()
    return float(np.sum(tc))/len(tc)


def get_ms_time_resolution(ms):
    '''
    return timestep (in s) and ntimes
    '''
    import pyrap.tables as pt
    t = pt.table(ms, readonly=True, ack=False)
    
    times = t.getcol('TIME')
    t.close()
    
    unq_times = np.unique(times)
    ntimes = len(unq_times)
    unq_times = np.sort(unq_times)
    
    dt = np.diff(unq_times)
    
    
    print ('MS {} has {} timesteps of {:.2f}s'.format(ms, ntimes, dt[0]))
    
    # check that all timesteps are the same
    if not np.all(dt==dt[0]):
        print ('Not all timesteps in ms {} are the same, returning full array of timesteps'.format(ms))
        return dt, ntimes
    
    return dt[0], ntimes


def get_ms_freq_resolution(ms):
    '''
    return freq (in kHz) and nfreq
    '''
    t = pt.table(ms+'::SPECTRAL_WINDOW', readonly=True, ack=False)
    freq = t.getcol('CHAN_FREQ')[0]
    
    nfreqs = len(freq)
    freq = np.sort(freq)
    dfreq = np.diff(freq)
    
    print ('MS {} has {} freqsteps of {:.2f} kHz'.format(ms, nfreqs, dfreq[0]*1e3))
    
    # check that all timesteps are the same
    if not np.all(dfreq==dfreq[0]):
        print ('Not all freqsteps in ms {} are the same, returning full array of freqsteps'.format(ms))
        return dfreq, nfreqs
    
    t.close()
    return dfreq[0]*1e3, nfreqs


def get_timerange(ms):
    t = pt.table(ms +'/OBSERVATION', readonly=True, ack=False)
    t.close()
    return t.getcell('TIME_RANGE',0)


def get_uvw(ms,unit='m'):
    with pt.table(ms, ack = False) as t:
        uvw = t.getcol('UVW')
        wavelength = 2.99e8 / np.mean(t.SPECTRAL_WINDOW[0]['CHAN_FREQ'])
        #uvw = np.linalg.norm(uvw, axis=1)
        if unit == 'lambda':
            uvw = uvw /wavelength
        return uvw
        
        
def get_pntg(ms):
    """
    Get phase centre
    """
    field_no = 0
    ant_no   = 0
    with pt.table(ms + "/FIELD", ack = False) as field_table:
        direction = field_table.getcol("PHASE_DIR")
    RA        = direction[ant_no, field_no, 0]
    Dec       = direction[ant_no, field_no, 1]

    if (RA < 0):
        RA += 2 * np.pi
    
    return (RA,Dec)
    
def get_antpos(ms):
    me = pm.measures()
    with pt.table(ms + "/ANTENNA", ack = False) as ant_table:
        pos = ant_table.getcol('POSITION')
        x = qa.quantity( pos[:,0], 'm' )
        y = qa.quantity( pos[:,1], 'm' )
        z = qa.quantity( pos[:,2], 'm' )
        position =  me.position( 'wgs84', x, y, z )
    return position
    
def get_elev1(ms):
    me = pm.measures()
    pos = get_antpos(ms)
    me.doframe(pos)
    
    # Get the first pointing of the first antenna
    with  pt.table(ms + '/FIELD', ack=False) as field_table:
        field_no = 0
        phase_dir = field_table.getcol('PHASE_DIR')
        ra = phase_dir[ 0, field_no, 0 ]
        if ra<0: ra += 2*numpy.pi
        dec = phase_dir[ 0, field_no, 1 ]
    # Get a ordered list of unique time stamps from the measurement set
    #time_table = pt.taql('select TIME from $1 orderby distinct TIME', tables=[ms])
    with pt.table(ms, readonly=True, ack=False) as t:
        time = t.getcol('TIME')
    
    time1 = time/3600.0
    time1 = time1 - np.floor(time1[0]/24)*24
    ra_qa  = qa.quantity( ra, 'rad' )
    dec_qa = qa.quantity( dec, 'rad' )
    pointing =  me.direction('j2000', ra_qa, dec_qa)
    t = qa.quantity(time[0], 's')
    t1 = me.epoch('utc', t)
    me.doframe(t1)
    # Loop through all time stamps and calculate the elevation of the pointing
    elevation = []
    for t in time:
        t_qa = qa.quantity(t, 's')
        t1 = me.epoch('utc', t_qa)
        me.doframe(t1)
        a = me.measure(pointing, 'azel')
        el = a['m1']
        elevation.append(el['value']/np.pi*180)
    elevation = np.array(elevation)
    return time1, elevation


def get_elevation(ms):
    ''' get time and elevation the astropy way
    input:
      ms - measurement set name
    returns:
      time - astropy.time.core.Time
      elevation - astropy.coordinates.angles.Latitude'''
    with pt.table(ms, readonly=True, ack=False) as t:
        time = t.getcol('TIME')
        time = np.unique(time)
        atime=Time(time*u.s,format='mjd',scale='utc')
    with pt.table(ms + '/ANTENNA', ack=False) as ant_table:
        pos = ant_table.getcol('POSITION')
        ant_id = 0
        x =  pos[ant_id,0]
        y =  pos[ant_id,1]
        z =  pos[ant_id,2]
        array_position = EarthLocation(x=x*u.m, y=y*u.m, z=z*u.m) # ellipsoid='wgs84')
    with  pt.table(ms + '/FIELD', ack=False) as field_table:
        field_id = 0
        phase_dir = field_table.getcol('PHASE_DIR')
        ra = phase_dir[ 0, field_id, 0 ]
        if ra<0: ra += 2*numpy.pi
        dec = phase_dir[ 0, field_id, 1 ]
        aphase_dir = SkyCoord(ra,dec,unit=u.rad)
        
    
    aphase_dir_altaz = aphase_dir.transform_to(AltAz(obstime=atime,location=array_position))
    return atime, aphase_dir_altaz.alt

