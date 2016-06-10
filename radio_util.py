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
    import numpy

    def prime_factors(n, douniq=True):
        """ Return the prime factors of the given number. """
        factors = []
        lastresult = n
        sqlast=int(numpy.sqrt(n))+1
        if n == 1:
            return [1]
        c=2
        while 1:
                if (lastresult == 1) or (c > sqlast):
                    break
                sqlast=int(numpy.sqrt(lastresult))+1
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
        return  numpy.unique(factors).tolist() if douniq else factors

    n = int(size)
    if (n%2 != 0):
        n+=1
    fac=prime_factors(n, False)
    for k in range(len(fac)):
        if (fac[k] > 7):
            val=fac[k]
            while (numpy.max(prime_factors(val)) > 7):
                val +=1
            fac[k]=val
    newlarge=numpy.product(fac)
    for k in range(n, newlarge, 2):
        if ((numpy.max(prime_factors(k)) < 8)):
            return k
    return newlarge



