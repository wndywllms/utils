import numpy as np

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

    print 'Time averaging should be less than %s'%delta_T


    # Time average smearing amplitude loss in http://adsabs.harvard.edu/full/1989ASPC....6..247B
    Reduction = 1-1.22E-9*(delta_Theta/resolution)**2.0 * delta_T**2.0

    print 'At radius %s a source will only have %s percent of its flux if data smoothed in time to %s'%(delta_Theta,Reduction,delta_T)
    
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

    print 'Bandwidth averaging should be much less than %s'%delta_freq

    
    beta = (delta_freq/freq) * (delta_Theta/resolution)
    gamma = 2*(np.log(2)**0.5)
    Reduction = ((np.pi**0.5)/(gamma * beta)) * (scipy.special.erf(beta*gamma/2.0))

    print 'At radius %s a source will only have %s percent of its flux if data smoothed in freq to %s'%(delta_Theta,Reduction,delta_freq)

    return delta_freq

def bandwidth_smearing2(delta_freq,freq,resolution,delta_Theta):.
    #http://www.cv.nrao.edu/course/astr534/Interferometers1.html
    # Output and input units are equal 
    # Condition that delta_Theta * delta_freq << Freq * resolution
    # where delta_Theta is the offset from the pointing centre and delta_freq is the bandwidth smoothing.

    # Given resolution, freq and offset this gives the condition for the delta_freq
    import scipy.special


    beta = (delta_freq/freq) * (delta_Theta/resolution)
    gamma = 2*(np.log(2)**0.5)
    Reduction = ((np.pi**0.5)/(gamma * beta)) * (scipy.special.erf(beta*gamma/2.0))

    return Reduction

