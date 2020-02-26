from astropy.coordinates import SkyCoord
import numpy as np

import shapely
from shapely.geometry import Polygon
from shapely.ops import cascaded_union

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


def ellipse(x0,y0,a,b,pa,n=200):
    theta=np.linspace(0,2*np.pi,n,endpoint=False)
    st = np.sin(theta)
    ct = np.cos(theta)
    pa = np.deg2rad(pa+90)
    sa = np.sin(pa)
    ca = np.cos(pa)
    p = np.empty((n, 2))
    p[:, 0] = x0 + a * ca * ct - b * sa * st
    p[:, 1] = y0 + a * sa * ct + b * ca * st
    return Polygon(p)


class Make_Shape(object):
    '''From lotss-catalogue/process_lgz.py -- Basic idea taken from remove_lgz_sources.py -- maybe should be merged with this one day
    but the FITS keywords are different.
    '''
    def __init__(self,clist):
        '''
        clist: a list of components that form part of the source, with RA, DEC, DC_Maj...
        '''
        ra=np.mean(clist['RA'])
        dec=np.mean(clist['DEC'])

        ellist=[]
        for r in clist:
            n_ra=r['RA']
            n_dec=r['DEC']
            x=3600*np.cos(dec*np.pi/180.0)*(ra-n_ra)
            y=3600*(n_dec-dec)
            newp=ellipse(x,y,r['DC_Maj']+0.1,r['DC_Min']+0.1,r['DC_PA'])
            ellist.append(newp)
        self.cp=cascaded_union(ellist)
        self.ra=ra
        self.dec=dec
        self.h=self.cp.convex_hull
        a=np.asarray(self.h.exterior.coords)
        #for i,e in enumerate(ellist):
        #    if i==0:
        #        a=np.asarray(e.exterior.coords)
        #    else:
        #        a=np.append(a,e.exterior.coords,axis=0)
        mdist2=0
        bestcoords=None
        for r in a:
            dist2=(a[:,0]-r[0])**2.0+(a[:,1]-r[1])**2.0
            idist=np.argmax(dist2)
            mdist=dist2[idist]
            if mdist>mdist2:
                mdist2=mdist
                bestcoords=(r,a[idist])
        self.mdist2=mdist2
        self.bestcoords=bestcoords
        self.a=a

    def length(self):
        return np.sqrt(self.mdist2)

    def pa(self):
        p1,p2=self.bestcoords
        dp=p2-p1
        angle=(180*np.arctan2(dp[1],dp[0])/np.pi)-90
        if angle<-180:
            angle+=360
        if angle<0:
            angle+=180
        return angle

    def width(self):
        p1,p2=self.bestcoords
        d = np.cross(p2-p1, self.a-p1)/self.length()
        return 2*np.max(d)
