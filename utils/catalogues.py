from astropy.coordinates import SkyCoord
import numpy as np


def name_from_coords1(ra,dec,prefix):
    sc = SkyCoord(ra,dec,frame='icrs', unit='deg')
    sc = sc.to_string(style='hmsdms',sep='',precision=2)
    name = prefix+sc.replace(' ','')[:-1]
    return name
    

def name_from_coords(ra,dec,prefix=''):
    '''
    wrapper for name_from_coords1 to handle arays of ra,dec
    '''
    
    if isinstance(ra, float):
        sc = SkyCoord(ra,dec,frame='icrs', unit='deg')
        sc = sc.to_string(style='hmsdms',sep='',precision=2)
        name = prefix+sc.replace(' ','')[:-1]
        return name_from_coords1(ra,dec,prefix=prefix)
    else:
        names = []
        for rai, deci in zip(ra,dec):
            names.append(name_from_coords1(rai,deci,prefix=prefix)) 
        return np.array(names)
