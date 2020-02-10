# Various routines for makeing fits cutouts
# using just a fits file
# or varius cutout servers

from astroquery.skyview import SkyView
import astropy.coordinates as coord
import astropy.units as u

from lxml import html
import requests
import shutil

import subprocess as sub
import numpy as np
import astropy.io.fits as pf
#import pywcs as pw
import astropy.wcs as pw
import aplpy as ap
import os
from .sky_util import *

VERBOSE = 10
clobber=True

def download_file(url,outname):
    if os.path.isfile(outname):
        print('File',outname,'already exists, skipping')
    else:
        print('Downloading',outname)
        while True:
            try:
                response = requests.get(url, stream=True,verify=False,timeout=120)
                if response.status_code!=200:
                    print('Warning, HTML status code',response.status_code)
                    if response.status_code>=500 and response.status_code<600:
                        raise RuntimeError('Retry!')
            except requests.exceptions.ConnectionError:
                print('Connection error! sleeping 60 seconds before retry...')
                sleep(60)
            except RuntimeError:
                print('Transient error reported, retrying after sleep')
                sleep(60)
            except requests.exceptions.Timeout:
                print('Timeout: retrying download')
            else:
                break
        with open(outname, 'wb') as out_file:
            shutil.copyfileobj(response.raw, out_file)
        del response
        

def plot_image(fits,ax, rms=np.nan, F=np.nan, cont=None, contcol='r', stretch='sqrt'):
  ax.set_frame_color('k')
  ax.set_tick_color('k')
  ax.tick_labels.set_font(size='xx-small')
  ax.tick_labels.set_xformat('ddd.dd')
  ax.tick_labels.set_yformat('ddd.dd')
  
  
  head = pf.getheader(fits)
  try:
    image = pf.getdata(fits)
  except:
    print('problem with image: ',fits)
    return
  ximgsize = head.get('NAXIS1')
  yimgsize = head.get('NAXIS2')
  pixscale = head.get('CDELT1')   
   
    
  if np.isnan(rms):
    rms = np.std(image)
    av = np.median(image)
    n_one = np.ones(image.shape)
    for i in range(5):
      #print rms, av
      within1sigma = np.where((image - av*n_one) < 3.*rms*n_one)
      av = np.median(image[within1sigma])
      rms = np.std(image[within1sigma])
  if np.isnan(F):
    F = image.max()
  
  
  dat = pf.getdata(fits)
  mean = np.median(dat)
  sig = np.std(dat)
  for i in range(10):
    dat = np.ma.masked_where(abs(dat - mean) > 5.*sig,dat).compressed()
    mean = np.median(dat)
    sig = np.std(dat)
  #, vmid = mean+20*sig
  ax.show_grayscale(vmin = mean-0.5*sig, vmax = mean+15*sig, stretch=stretch, invert=True)
  

  
  if cont!=None:
    image = pf.getdata(cont)
    rms = np.std(image)
    av = np.median(image)
    n_one = np.ones(image.shape)
    for i in range(5):
      within1sigma = np.where((image - av*n_one) < 3.*rms*n_one)
      av = np.median(image[within1sigma])
      rms = np.std(image[within1sigma])
    cont_levels = 3.*rms*np.array([-2.**(1 * j ) for j in range(0, 2 )] + [ 2.**(1. * j ) for j in range(0, 10 ) ] )
    ax.show_contour(cont, hdu=0, layer='contour', levels=cont_levels, filled=False, cmap=None, colors=contcol, returnlevels=False, convention=None, slices=[0,1], smooth=None, kernel='gauss')
  #else:
    #cont_levels = 3.*rms*np.array([-2.**(1 * j ) for j in range(0, 2 )] + [ 2.**(1. * j ) for j in range(0, 10 ) ] )
    #ax.show_contour(fits, hdu=0, layer='contour', levels=cont_levels, filled=False, cmap=None, colors=contcol, returnlevels=False, convention=None, slices=[0,1], smooth=None, kernel='gauss')
    
  return


def postage(fitsim,postfits,ra,dec,s=2./60, verbose=0):
  '''
  input:
    fitsim   - name of fits file from which to make cutout
    postfits - name of cutout
    ra       - degrees
    dec      - degrees
    s        - size of cutout in degrees (default 2')
  '''
  
  head = pf.getheader(fitsim)
  
  hdulist = pf.open(fitsim)
  # Parse the WCS keywords in the primary HDU
  wcs = pw.WCS(hdulist[0].header)

  # Some pixel coordinates of interest.
  #skycrd = np.array([ra,dec])
  skycrd = np.array([[ra,dec,0,0]], np.float_)

  # Convert pixel coordinates to world coordinates
  pixel = wcs.wcs_sky2pix(skycrd, 1)

  x = pixel[0][0]
  y = pixel[0][1]
  pixsize = abs(wcs.wcs.cdelt[0])
  if np.isnan(s):
    s = 25.
  N = s/pixsize
  if verbose > 2:
      print("making cutout %s from fits %s at ra,dec=%.5f, %.5f and size %.3f'" %(fitsim,postfits,ra,dec,s*60))
      print('x=%.5f, y=%.5f, N=%i' %(x,y,N))

  ximgsize = head.get('NAXIS1')
  yimgsize = head.get('NAXIS2')

  if x ==0:
    x = ximgsize/2
  if y ==0:
    y = yimgsize/2

  offcentre = False
  # subimage limits: check if runs over edges
  xlim1 =  x - (N/2)
  if(xlim1<1):
    xlim1=1
    offcentre=True
  xlim2 =  x + (N/2)
  if(xlim2>ximgsize):
    xlim2=ximgsize
    offcentre=True
  ylim1 =  y - (N/2)
  if(ylim1<1):
    ylim1=1
    offcentre=True
  ylim2 =  y + (N/2)
  if(ylim2>yimgsize):
    offcentre=True
    ylim2=yimgsize

  xl = int(xlim1)
  yl = int(ylim1)
  xu = int(xlim2)
  yu = int(ylim2)
  if verbose > 2:
    print('postage stamp is %i x %i pixels' %(xu-xl,yu-yl))

  # make fits cutout
  inps = fitsim + '[%0.0f:%0.0f,%0.0f:%0.0f]' %(xl,xu,yl,yu)

  # using fitscopy make the cutout
  if os.path.isfile(postfits): os.system('rm '+postfits)
  os.system( 'fitscopy %s %s' %(inps,postfits) )
  return



def cutout_from_local_file(fitscut, ra, dec, imsize, local_file="", clobber=False):
    """
args:
    fitscut - name of cutout
    ra      - degrees
    dec     - degrees
    imsize  - size of cutout in degrees
kwargs:
    local_file     - name of fits file from which to make cutout
    """
    import pyfits as pf
    import pywcs as pw
    
    if os.path.exists(fitscut):
        if clobber:
            os.system('rm '+fitscut)
        else:
            print('file exists and clobber is False: ',fitscut)
            return
        
    
    
    head = pf.getheader(local_file)
    
    hdulist = pf.open(local_file)
    # Parse the WCS keywords in the primary HDU
    wcs = pw.WCS(hdulist[0].header)

    # Some pixel coordinates of interest.
    #skycrd = np.array([ra,dec])
    skycrd = np.array([[ra,dec,0,0]], np.float_)

    # Convert pixel coordinates to world coordinates
    pixel = wcs.wcs_sky2pix(skycrd, 1)

    x = np.floor(pixel[0][0])
    y = np.floor(pixel[0][1])
    pixsize = abs(wcs.wcs.cdelt[0])
    N = imsize/pixsize
    #if VERBOSE > 2:
    print("making cutout %s from fits %s at ra,dec=%.5f, %.5f and size %.3f'" %(local_file,fitscut,ra,dec,imsize*60))
    #if VERBOSE > 5:
    print('x=%.5f, y=%.5f, N=%i' %(x,y,N))

    ximgsize = head.get('NAXIS1')
    yimgsize = head.get('NAXIS2')

    if x ==0:
        x = ximgsize/2
    if y ==0:
        y = yimgsize/2

    offcentre = False
    # subimage limits: check if runs over edges
    xlim1 =    x - (N/2)
    if(xlim1<1):
        xlim1=1
        offcentre=True
    xlim2 =    x + (N/2)
    if(xlim2>ximgsize):
        xlim2=ximgsize
        offcentre=True
    ylim1 =    y - (N/2)
    if(ylim1<1):
        ylim1=1
        offcentre=True
    ylim2 =    y + (N/2)
    if(ylim2>yimgsize):
        offcentre=True
        ylim2=yimgsize

    xl = int(xlim1)
    yl = int(ylim1)
    xu = int(xlim2)
    yu = int(ylim2)
    if VERBOSE > 5:
        print('postage stamp is %i x %i pixels' %(xu-xl,yu-yl))

    # make fits cutout
    inps = local_file + "[%0.0f:%0.0f,%0.0f:%0.0f]" %(xl,xu,yl,yu)

    # using fitscopy make the cutout
    if os.path.isfile(fitscut): os.system('rm '+fitscut)
    os.system( 'fitscopy %s %s' %(inps,fitscut) )
    return

def cutout_from_local_files2(fitscut, ra, dec, imsize, local_file_path="", local_file_list="", clobber=False):
    """
args:
    fitscut - name of cutout
    ra      - degrees
    dec     - degrees
    imsize  - size of cutout in degrees
kwargs:
    local_file     - name of fits file from which to make cutout
    """

    if os.path.exists(fitscut):
        if clobber:
            os.system('rm '+fitscut)
        else:
            print('file exists and clobber is False: ',fitscut)
            return

    import astropy.coordinates as ac
    from astropy.table import Table
    
    t = Table.read(local_file_list, format='ascii')
    C = ac.SkyCoord(ra,dec,unit='deg')
    idx, sep, dist =  ac.match_coordinates_sky(C, ac.SkyCoord(t['ra'],t['dec'],unit='deg'))
    fitsname = t['name'][idx]
    
    print('nearest fitsimage is',fitsname)
    
    cutout_from_local_file(fitscut, ra, dec, imsize, local_file=local_file_path+'/'+fitsname)
    
    return


def cutout_from_local_files(fitscut, ra, dec, imsize, local_file_path="", local_file_list="", clobber=False):
    """
args:
    fitscut - name of cutout
    ra      - degrees
    dec     - degrees
    imsize  - size of cutout in degrees
kwargs:
    local_file     - name of fits file from which to make cutout
    """
    import pyfits as pf
    import pywcs as pw
    
    
    if os.path.exists(fitscut):
        if clobber:
            os.system('rm '+fitscut)
        else:
            print('file exists and clobber is False: ',fitscut)
            return
    
    imdat = np.genfromtxt(local_file_list, dtype='S15,f,f,f,f', names=['images','ramin','ramax','decmin','decmax'])

    ind = np.where((ra>=imdat['ramin']) & (ra<imdat['ramax']) & (dec>=imdat['decmin']) & (dec<imdat['decmax']))[0][0]
    local_file = local_file_path+'/'+imdat['images'][ind]

    head = pf.getheader(local_file)
    
    hdulist = pf.open(local_file)
    # Parse the WCS keywords in the primary HDU
    wcs = pw.WCS(hdulist[0].header)

    # Some pixel coordinates of interest.
    #skycrd = np.array([ra,dec])
    skycrd = np.array([[ra,dec,0,0]], np.float_)

    # Convert pixel coordinates to world coordinates
    pixel = wcs.wcs_sky2pix(skycrd, 1)

    x = pixel[0][0]
    y = pixel[0][1]
    #pixsize = abs(wcs.wcs.cdelt[0])
    pixsize = abs(wcs.wcs.cd[0][0])
    N = imsize/pixsize
    #if VERBOSE > 2:
    print("making cutout %s from fits %s at ra,dec=%.5f, %.5f and size %.3f'" %(local_file,fitscut,ra,dec,imsize*60))
    #if VERBOSE > 5:
    print('x=%.5f, y=%.5f, N=%i' %(x,y,N))

    ximgsize = head.get('NAXIS1')
    yimgsize = head.get('NAXIS2')

    if ximgsize is None:
        dat = pf.getdata(local_file)
        ximgsize, yimgsize = dat.shape

    if x ==0:
        x = ximgsize/2
    if y ==0:
        y = yimgsize/2

    offcentre = False
    # subimage limits: check if runs over edges
    xlim1 =    x - (N/2)
    if(xlim1<1):
        xlim1=1
        offcentre=True
    xlim2 =    x + (N/2)
    if(xlim2>ximgsize):
        xlim2=ximgsize
        offcentre=True
    ylim1 =    y - (N/2)
    if(ylim1<1):
        ylim1=1
        offcentre=True
    ylim2 =    y + (N/2)
    if(ylim2>yimgsize):
        offcentre=True
        ylim2=yimgsize

    xl = int(xlim1)
    yl = int(ylim1)
    xu = int(xlim2)
    yu = int(ylim2)
    if VERBOSE > 5:
        print('postage stamp is %i x %i pixels' %(xu-xl,yu-yl))

    # make fits cutout
    inps = local_file + "[%0.0f:%0.0f,%0.0f:%0.0f]" %(xl,xu,yl,yu)

    # using fitscopy make the cutout
    if os.path.isfile(fitscut): os.system('rm '+fitscut)
    os.system( 'fitscopy %s %s' %(inps,fitscut) )
    return


def download_panstarrs(fitsname,ra,dec,f='i',imsize=0.08, clobber=False):
    '''
    imsize in deg
    '''
    
    if os.path.exists(fitsname):
        if clobber:
            os.system('rm '+fitsname)
        else:
            print('file exists and clobber is False: ',fitsname)
            return
        
    s = int(imsize*3600.*4)  #arcsec to pix seems to be 4pix=1arcsec
    
    wgeturl = "http://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos={ra:f}+{dec:f}&filter=color&filter={f:s}&filetypes=stack&auxiliary=data&size={s:d}&output_size=0&verbose=0&autoscale=99.500000&catlist=".format(ra=ra,dec=dec,f=f,s=s)
    cmd = ['wget', wgeturl, '-O', 'ttt']
    print(' '.join(cmd))
    p = sub.Popen(cmd)
    p.wait()
    with open('ttt','r') as ff:
        lines = ff.readlines()
    for line in lines:
        if 'Download' in line:
            break
    i = line.find('Download FITS cutout')
    
    line1 = line[i:]
    i = line1.find('http:')
    if i == -1:
        i = line1.find('href="//')
        line2 = line1[i+8:]
        i = line2.find('">')
        line3 = line2[:i]
        fits = line3
        print(fits)
    else:
        line2 = line1[i:]
        i = line2.find('">')
        line3 = line2[:i]
        fits = line3
        print(fits)
    #fitsname = 'panstars_{f:s}_{ra:f}+{dec:f}.fits'.format(name=fits,f=f,ra=ra,dec=dec)
    cmd = ['wget',fits, '-O', fitsname ]
    print(' '.join(cmd))
    p=sub.Popen(cmd)
    p.wait()
    try:
        pf.open(fitsname)
    except:
        print('error downloading panstarrs cutout')
        return -1
    os.system('rm -rf ttt')
    return fitsname

def get_first(ra,dec):
    url="http://archive.stsci.edu/"
    page=requests.get(url+"vlafirst/search.php?RA=%.7f&DEC=%.6f&Radius=30.0&action=Search" % (ra,dec),verify=False)
    print(page.status_code)

    tree=html.fromstring(page.text)
    table=tree.xpath('//tbody')
    links=[]
    dists=[]
    for row in table[0].getchildren():
        td=row.getchildren()
        links.append(td[0].getchildren()[0].attrib['href'])
        dists.append(float(td[8].text))

    index=np.argmin(dists)
    path=links[index]
   
    outname=path.split('/')[-1]
    download_file(url+path,outname)
    return outname


def get_nvss(ra,dec,size=1000):
    # uses skyview
    coords=coord.SkyCoord(ra, dec, unit=(u.deg, u.deg))
    paths = SkyView.get_images(position=coords.to_string('hmsdms'),survey='NVSS',width=size*u.arcsec)
    hdu = paths[0]
    hdu[0].header['BMAJ']=45.0/3600.0
    hdu[0].header['BMIN']=45.0/3600.0
    hdu[0].header['BPA']=0
    hdu[0].header['RESTFREQ']=1.4e9
    filename='NVSS-'+coords.to_string('hmsdms').replace(' ','')+'.fits'
    hdu.writeto(filename)
    return filename


def cutout_from_server(fitscut, url="", clobber=False):
    """ get a fits cutout from a server
args:
    fitscut - name of cutout
kwargs:
    url     - server url
    """
    import pyfits as pf
    
    
    if os.path.exists(fitscut):
        if clobber:
            os.system('rm '+fitscut)
        else:
            print('file exists and clobber is False: ',fitscut)
            return
    
    if VERBOSE > 2:
        print("cutout from server %s" %(fitscut))
  
    if (not os.path.isfile(fitscut)) or clobber:
        if VERBOSE > 5: print('...downloading')
        N = np.random.uniform(1,100000)
        #print url
        base = ["wget", "-q" ,"-O", 'tempout%05i' %(N)]
        # with timeout
        base = ["timeout", "60", "wget", "-q" ,"-O", 'tempout%05i' %(N)]
        base.extend(url.split())
        if VERBOSE > 5:
            print(' '.join(base))
        sub.Popen( base ).wait()
        size = os.path.getsize('tempout%05i' %(N))

        
        try:
            print('try open')
            pf.open('tempout%05i' %(N))
        except IOError:
            #sub.Popen(['rm', 'tempout%05i' %(N)]).wait()
            print("Failed")
            return False
        except:
            #sub.Popen(['rm', 'tempout%05i' %(N)]).wait()
            return False
        
        #sub.Popen(['mv', 'tempout%05i' %(N), fitscut]).wait()
        os.system('mv tempout%05i' %(N)+' '+fitscut)
        #print fitscut
        return True


def get_WSRT_cutout(fitscut, ra, dec, imsize, clobber=False):
    """ get a fits cutout from WSRT server
args:
    fitscut - name of cutout
    ra      - degrees
    dec     - degrees
    imsize  - size of cutout in degrees
returns
    result  - success?
    """
    
    if os.path.exists(fitscut):
        if clobber:
            os.system('rm '+fitscut)
        else:
            print('file exists and clobber is False: ',fitscut)
            return
    
    sra = ra_to_str( ra ).replace(':','+')
    sdec = dec_to_str( dec ).replace(':','+')
    imsize = imsize*60.  # in arcmin
    
    url = "http://www.astron.nl/wow/testcode.php?POS={ra}+%s2B{dec}&Equinox=J2000&SIZE={imsize:.1f}%2C+{imsize:.1f}&cutout=1&surv_id=5".format(ra=sra, dec=sdec, imsize=imsize)
    
    result = cutout_from_server(fitscut, url=url)
    return result
    
def get_NDWFS_cutout(fitscut, ra, dec, imsize, band="I", clobber=False):
    """ get a fits cutout from NDWFS server
args:
    fitscut - name of cutout
    ra      - degrees
    dec     - degrees
    imsize  - size of cutout in degrees
kwargs
    band    - NDWFS band
returns
    result  - success?
    """
    
    if os.path.exists(fitscut):
        if clobber:
            os.system('rm '+fitscut)
        else:
            print('file exists and clobber is False: ',fitscut)
            return
    
    print(ra,dec)
    sra = ra_to_str( ra )
    sdec = dec_to_str( dec )
    imsize = imsize*60.  # in arcmin
    
    #url = "http://archive.noao.edu/ndwfs/cutout.php?ra={ra}&dec={dec}&rawidth={imsize:.1f}&decwidth={imsize:.1f}&filter={band}".format(ra=sra, dec=sdec, imsize=imsize, band=band)
    
    url = "http://r2.sdm.noao.edu/ndwfs/cutout.php?ra={ra}&dec={dec}&rawidth={imsize:.1f}&decwidth={imsize:.1f}&filter={band}".format(ra=sra, dec=sdec, imsize=imsize, band=band)
    
    result = cutout_from_server(fitscut, url=url)
    return result

def get_SDWFS_cutout(fitscut, ra, dec, imsize, band="I1", clobber=False):
    """ get a fits cutout from SDWFS server
args:
    fitscut - name of cutout
    ra      - degrees
    dec     - degrees
    imsize  - size of cutout in degrees
kwargs
    band    - SDWFS band
returns
    result  - success?
    """
    
    if os.path.exists(fitscut):
        if clobber:
            os.system('rm '+fitscut)
        else:
            print('file exists and clobber is False: ',fitscut)
            return
        
    url = "http://irsa.ipac.caltech.edu/cgi-bin/Subimage/nph-subimage?origfile=/irsadata/SPITZER/SDWFS//images/combined_epochs/{band}_bootes.v32.fits&ra={ra:f}&dec={dec:f}&xsize={imsize:f}".format(band=band, ra=ra, dec=dec, imsize=imsize)
    
    result = cutout_from_server(fitscut, url=url)
    return result

def get_first_cutout(fitscut, ra, dec, imsize, clobber=False):
    """ get a fits cutout from FIRST server
args:
    fitscut - name of cutout
    ra      - degrees
    dec     - degrees
    imsize  - size of cutout in degrees
returns
    result  - success?
    """
    
    if os.path.exists(fitscut):
        if clobber:
            os.system('rm '+fitscut)
        else:
            print('file exists and clobber is False: ',fitscut)
            return
    
    sra = ra_to_str( ra )
    sdec = dec_to_str( dec )
    imsize = imsize * 60.
    
    url = "third.ucllnl.org/cgi-bin/firstimage?RA={ra} {dec}&Dec=&Equinox=J2000&ImageSize={imsize:.1f}&MaxInt=10&FITS=1&Download=1".format(ra=sra, dec=sdec, imsize=imsize)
    
    print(url)
    result = cutout_from_server(fitscut, url=url)
    return result


def get_nvss_cutout(fitscut, ra, dec, imsize):
    """ get a fits cutout from NVSS server
args:
    fitscut - name of cutout
    ra      - degrees
    dec     - degrees
    imsize  - size of cutout in degrees
returns
    result  - success?
    """
    sra = ra_to_str( ra )
    sdec = dec_to_str( dec )
    imsize = imsize * 60.
    
    url = "third.ucllnl.org/cgi-bin/firstimage?RA={ra} {dec}&Dec=&Equinox=J2000&ImageSize={imsize:.1f}&MaxInt=10&FITS=1&Download=1".format(ra=sra, dec=sdec, imsize=imsize)
    
    print(url)
    result = cutout_from_server(fitscut, url=url)
    return result

def get_NDWFS_cutout_MAGES(fitsname, ra,dec, user, password, imsize = 2., band='I', verbose=0, clobber=False):
    
  if os.path.exists(fitsname):
    if clobber:
        os.system('rm '+fitsname)
    else:
        print('file exists and clobber is False: ',fitsname)
        return
    
  print('cutout %s' %(fitsname))
  #imsize in arcmin
  sra = ra_to_str( ra )
  sdec = dec_to_str( dec )
  
  
  bands = ['FUV', 'NUV', 'Bw', 'R', 'I', 'K', 'J', 'Ks', 'z', '3.6', '4.5', '5.8', '8.0', '24']
  ndwfs_cutout_bands = np.array((np.nan, np.nan, 1, 2, 3, 4, 28, 30, np.nan, 5, 6, 7, 8, 25))
  bandid = ndwfs_cutout_bands[bands.index(band)]
  if verbose > 2:
    print(band, bandid)
  '''
1 'Bw'
2 'R'
3 'I'
4 'K'
5 'IRAC_3.6 all epochs'
6 'IRAC_4.5 all epochs'
7 'IRAC_5.8 all epochs'
8 'IRAC_8.0 all epochs'
9 'IRAC_3.6 epoch 1'
10 'IRAC_4.5 epoch 1'
11 'IRAC_5.8 epoch 1'
12 'IRAC_8.0 epoch 1'
13 'IRAC_3.6 epoch 2'
14 'IRAC_4.5 epoch 2'
15 'IRAC_5.8 epoch 2'
16 'IRAC_8.0 epoch 2'
17 'IRAC_3.6 epoch 3'
18 'IRAC_4.5 epoch 3'
19 'IRAC_5.8 epoch 3'
20 'IRAC_8.0 epoch 3'
21 'IRAC_3.6 epoch 4'
22 'IRAC_4.5 epoch 4'
23 'IRAC_5.8 epoch 4'
24 'IRAC_8.0 epoch 4'
25 'MIPS 24'
26 'MIPS 70'
27 'MIPS 160'
28 'J'
29 'H'
30 'Ks'
  '''
  if np.isnan(bandid):
      print('not in mages')
      return
  
  url=" --user=%s --password=%s  http://www.noao.edu/ndwfs/cutout.dh5.php?ra=%s&dec=%s&rawidth=%.1f&decwidth=%.1f&s_filter_id=%i" %(user, password, sra, sdec, imsize, imsize, bandid)
  
  result = cutout_from_server(fitsname, url=url)
  return result
  
