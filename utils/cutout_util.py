# Various routines for makeing fits cutouts
# using just a fits file
# or varius cutout servers

import subprocess as sub
import numpy as np
import pyfits as pf
import pywcs as pw
import aplpy as ap
import os
from sky_util import *

VERBOSE = 10
clobber=True

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
    print 'problem with image: ',fits
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
  #dat = np.ma.masked_where(dat<1.,dat).compressed()
  #dat = np.ma.masked_where(dat>10000.,dat).compressed()
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
  
  ##s in degrees
  #tra = ra/15.
  #dra = math.floor(tra)
  #dramin = math.floor((tra - dra)*60.)
  #drasec = ((tra - dra)*60. - math.floor((tra - dra)*60.))*60.
  #sra = '%02i%02i%02i' %(dra,dramin,drasec)
  #ddec = math.floor(dec)
  #ddecmin = math.floor((dec - ddec)*60.)
  #ddecsec = ((dec - ddec)*60. - math.floor((dec - ddec)*60.))*60.
  #sdec = '%02i%02i' %(ddec,ddecmin)
  #Src_name = 'J%s+%s' %(sra,sdec)

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
      print "making cutout %s from fits %s at ra,dec=%.5f, %.5f and size %.3f'" %(fitsim,postfits,ra,dec,s*60)
      print 'x=%.5f, y=%.5f, N=%i' %(x,y,N)

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
    print 'postage stamp is %i x %i pixels' %(xu-xl,yu-yl)

  # make fits cutout
  inps = fitsim + '[%0.0f:%0.0f,%0.0f:%0.0f]' %(xl,xu,yl,yu)

  # using fitscopy make the cutout
  if os.path.isfile(postfits): os.system('rm '+postfits)
  os.system( 'fitscopy %s %s' %(inps,postfits) )
  return



def cutout_from_local_file(fitscut, ra, dec, imsize, local_file=""):
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
    print "making cutout %s from fits %s at ra,dec=%.5f, %.5f and size %.3f'" %(local_file,fitscut,ra,dec,imsize*60)
    #if VERBOSE > 5:
    print 'x=%.5f, y=%.5f, N=%i' %(x,y,N)

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
        print 'postage stamp is %i x %i pixels' %(xu-xl,yu-yl)

    # make fits cutout
    inps = local_file + "[%0.0f:%0.0f,%0.0f:%0.0f]" %(xl,xu,yl,yu)

    # using fitscopy make the cutout
    if os.path.isfile(fitscut): os.system('rm '+fitscut)
    os.system( 'fitscopy %s %s' %(inps,fitscut) )
    return

def cutout_from_local_files2(fitscut, ra, dec, imsize, local_file_path="", local_file_list=""):
    """
args:
    fitscut - name of cutout
    ra      - degrees
    dec     - degrees
    imsize  - size of cutout in degrees
kwargs:
    local_file     - name of fits file from which to make cutout
    """


    import astropy.coordinates as ac
    from astropy.table import Table
    
    t = Table.read(local_file_list, format='ascii')
    C = ac.SkyCoord(ra,dec,unit='deg')
    idx, sep, dist =  ac.match_coordinates_sky(C, ac.SkyCoord(t['ra'],t['dec'],unit='deg'))
    fitsname = t['name'][idx]
    
    print 'nearest fitsimage is',fitsname
    
    cutout_from_local_file(fitscut, ra, dec, imsize, local_file=local_file_path+'/'+fitsname)
    
    return


def cutout_from_local_files(fitscut, ra, dec, imsize, local_file_path="", local_file_list=""):
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
    print "making cutout %s from fits %s at ra,dec=%.5f, %.5f and size %.3f'" %(local_file,fitscut,ra,dec,imsize*60)
    #if VERBOSE > 5:
    print 'x=%.5f, y=%.5f, N=%i' %(x,y,N)

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
        print 'postage stamp is %i x %i pixels' %(xu-xl,yu-yl)

    # make fits cutout
    inps = local_file + "[%0.0f:%0.0f,%0.0f:%0.0f]" %(xl,xu,yl,yu)

    # using fitscopy make the cutout
    if os.path.isfile(fitscut): os.system('rm '+fitscut)
    os.system( 'fitscopy %s %s' %(inps,fitscut) )
    return


def download_panstarrs(fitsname,ra,dec,f='i',imsize=0.08):
    '''
    imsize in deg
    '''
    s = int(imsize*3600.*4)  #arcsec to pix seems to be 4pix=1arcsec
    
    wgeturl = "http://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos={ra:f}+{dec:f}&filter=color&filter={f:s}&filetypes=stack&auxiliary=data&size={s:d}&output_size=0&verbose=0&autoscale=99.500000&catlist=".format(ra=ra,dec=dec,f=f,s=s)
    cmd = ['wget', wgeturl, '-O', 'ttt']
    print ' '.join(cmd)
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
        print fits
    else:
        line2 = line1[i:]
        i = line2.find('">')
        line3 = line2[:i]
        fits = line3
        print fits
    #fitsname = 'panstars_{f:s}_{ra:f}+{dec:f}.fits'.format(name=fits,f=f,ra=ra,dec=dec)
    cmd = ['wget',fits, '-O', fitsname ]
    print ' '.join(cmd)
    p=sub.Popen(cmd)
    p.wait()
    try:
        pf.open(fitsname)
    except:
        print 'error downloading panstarrs cutout'
        return -1
    os.system('rm -rf ttt')
    return fitsname

def cutout_from_server(fitscut, url=""):
    """ get a fits cutout from a server
args:
    fitscut - name of cutout
kwargs:
    url     - server url
    """
    import pyfits as pf
    
    if VERBOSE > 2:
        print "cutout from server %s" %(fitscut)
  
    if (not os.path.isfile(fitscut)) or clobber:
        if VERBOSE > 5: print '...downloading'
        N = np.random.uniform(1,100000)
        #print url
        base = ["wget", "-q" ,"-O", 'tempout%05i' %(N)]
        # with timeout
        base = ["timeout", "60", "wget", "-q" ,"-O", 'tempout%05i' %(N)]
        base.extend(url.split())
        if VERBOSE > 5:
            print ' '.join(base)
        sub.Popen( base ).wait()
        size = os.path.getsize('tempout%05i' %(N))
        #os.system('fitscheck tempout%05i > tempoutfits%05i' %(N,N))
        
        #with open('tempoutfits%05i' %(N)) as t:
            #fline = t.readline().strip()
            #if "Header missing END card" in fline:
                #print "Failed"
                #return False
            #elif "Empty or corrupt FITS file" in fline:
                #print "Failed"
                #return False
            #elif "Checksum not found" in fline:
                #print 'checksum error'
                #sub.Popen(['mv', 'tempout%05i' %(N), fitscut]).wait()
                ##print fitscut
                #return True
        
        try:
            print 'try open'
            pf.open('tempout%05i' %(N))
        except IOError:
            #sub.Popen(['rm', 'tempout%05i' %(N)]).wait()
            print "Failed"
            return False
        except:
            #sub.Popen(['rm', 'tempout%05i' %(N)]).wait()
            return False
        
        #sub.Popen(['mv', 'tempout%05i' %(N), fitscut]).wait()
        os.system('mv tempout%05i' %(N)+' '+fitscut)
        #print fitscut
        return True


def get_WSRT_cutout(fitscut, ra, dec, imsize):
    """ get a fits cutout from WSRT server
args:
    fitscut - name of cutout
    ra      - degrees
    dec     - degrees
    imsize  - size of cutout in degrees
returns
    result  - success?
    """
    
    sra = ra_to_str( ra ).replace(':','+')
    sdec = dec_to_str( dec ).replace(':','+')
    imsize = imsize*60.  # in arcmin
    
    url = "http://www.astron.nl/wow/testcode.php?POS={ra}+%s2B{dec}&Equinox=J2000&SIZE={imsize:.1f}%2C+{imsize:.1f}&cutout=1&surv_id=5".format(ra=sra, dec=sdec, imsize=imsize)
    
    result = cutout_from_server(fitscut, url=url)
    return result
    
def get_NDWFS_cutout(fitscut, ra, dec, imsize, band="I"):
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
    
    print ra,dec
    sra = ra_to_str( ra )
    sdec = dec_to_str( dec )
    imsize = imsize*60.  # in arcmin
    
    url = "http://archive.noao.edu/ndwfs/cutout.php?ra={ra}&dec={dec}&rawidth={imsize:.1f}&decwidth={imsize:.1f}&filter={band}".format(ra=sra, dec=sdec, imsize=imsize, band=band)
    
    result = cutout_from_server(fitscut, url=url)
    return result

def get_SDWFS_cutout(fitscut, ra, dec, imsize, band="I1"):
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
        
    url = "http://irsa.ipac.caltech.edu/cgi-bin/Subimage/nph-subimage?origfile=/irsadata/SPITZER/SDWFS//images/combined_epochs/{band}_bootes.v32.fits&ra={ra:f}&dec={dec:f}&xsize={imsize:f}".format(band=band, ra=ra, dec=dec, imsize=imsize)
    
    result = cutout_from_server(fitscut, url=url)
    return result

def get_first_cutout(fitscut, ra, dec, imsize):
    """ get a fits cutout from FIRST server
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
    
    print url
    result = cutout_from_server(fitscut, url=url)
    return result


#def get_nvss_cutout(fitscut, ra, dec, imsize):
    #""" get a fits cutout from NVSS server
#args:
    #fitscut - name of cutout
    #ra      - degrees
    #dec     - degrees
    #imsize  - size of cutout in degrees
#returns
    #result  - success?
    #"""
    #sra = ra_to_str( ra )
    #sdec = dec_to_str( dec )
    #imsize = imsize * 60.
    
    #url = "third.ucllnl.org/cgi-bin/firstimage?RA={ra} {dec}&Dec=&Equinox=J2000&ImageSize={imsize:.1f}&MaxInt=10&FITS=1&Download=1".format(ra=sra, dec=sdec, imsize=imsize)
    
    #print url
    #result = cutout_from_server(fitscut, url=url)
    #return result

def get_NDWFS_cutout_MAGES(fitsname, ra,dec, imsize = 2., band='I', verbose=0):
    
  print 'cutout %s' %(fitsname)
  #imsize in arcmin
  sra = ra_to_str( ra )
  sdec = dec_to_str( dec )
  
  
  bands = ['FUV', 'NUV', 'Bw', 'R', 'I', 'K', 'J', 'Ks', 'z', '3.6', '4.5', '5.8', '8.0', '24']
  ndwfs_cutout_bands = np.array((np.nan, np.nan, 1, 2, 3, 4, 28, 30, np.nan, 5, 6, 7, 8, 25))
  bandid = ndwfs_cutout_bands[bands.index(band)]
  if verbose > 2:
    print band, bandid
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
      print 'not in mages'
      return
  
  url=" --user=NDWFScutout --password=givemedata!  http://www.noao.edu/ndwfs/cutout.dh5.php?ra=%s&dec=%s&rawidth=%.1f&decwidth=%.1f&s_filter_id=%i" %(sra, sdec, imsize, imsize, bandid)
  
  result = cutout_from_server(fitsname, url=url)
  return result
  
  
  #if (not os.path.isfile(fitsname)) or clobber:
    #print '...downloading'
    #N = np.random.uniform(1,100000)
    #sub.Popen(["wget", "-q" ,"-O", 'tempout%05i' %(N) ,"--user=NDWFScutout", "--password=givemedata!" , url]).wait()
    #if verbose > 2:
        #print ' '.join(["wget", "-q" ,"-O", 'tempout%05i' %(N),"--user=NDWFScutout", "--password=givemedata!" , url])
    #size = os.path.getsize('tempout%05i' %(N))

    #if size > 100.:
      #sub.Popen(['mv', 'tempout%05i' %(N), fitsname]).wait()
      #print fitsname
    #else:
      #sub.Popen(['rm', 'tempout%05i' %(N)]).wait()
  #return 

    
#def get_NDWFS_cutout(ra,dec, fitsname, band='I', imsize = 2., verbose=0):
    
  #print 'cutout %s' %(fitsname)
  ##imsize in arcmin
  #sra = ra_to_str( ra )
  #sdec = dec_to_str( dec )
  #url = "http://archive.noao.edu/ndwfs/cutout.php?ra=%s&dec=%s&rawidth=%.1f&decwidth=%.1f&filter=%s"    %(sra, sdec, imsize, imsize, band)
  
  ##url="http://www.noao.edu/ndwfs/cutout-form.dh5.php?ra=%s&dec=%s&rawidth=%.1f&decwidth=%.1f&filter=%s"    %(sra, sdec, imsize, imsize, band)
  
  ##url="http://www.noao.edu/ndwfs/cutout.dh5.php?ra=14:31:18.2&dec=35:15:50&rawidth=2.0&decwidth=2.0&s_filter_id=25"
  
  #if (not os.path.isfile(fitsname)) or clobber:
    #print '...downloading'
    #N = np.random.uniform(1,100000)
    #sub.Popen(["wget", "-q" ,"-O", 'tempout%05i' %(N) , url]).wait()
    #if verbose > 2:
        #print ' '.join(["wget", "-q" ,"-O", 'tempout%05i' %(N) , url])
    #size = os.path.getsize('tempout%05i' %(N))

    #if size > 91.:
      #sub.Popen(['mv', 'tempout%05i' %(N), fitsname]).wait()
      #print fitsname
    #else:
      #sub.Popen(['rm', 'tempout%05i' %(N)]).wait()
  #return 

#def get_SDWFS_cutout(ra,dec, fitsname, band='I1', imsize = 1./60., verbose=0):
  #print 'cutout %s' %(fitsname)
  ##imsize in degrees
  #url = "http://irsa.ipac.caltech.edu/cgi-bin/Subimage/nph-subimage?origfile=/irsadata/SPITZER/SDWFS//images/combined_epochs/%s_bootes.v32.fits&ra=%f&dec=%f&xsize=%f"    %(band, ra, dec, imsize)
  
  #if (not os.path.isfile(fitsname)) or clobber:
    #print '...downloading'
    #N = np.random.uniform(1,100000)
    #sub.Popen(["wget", "-q" ,"-O", 'tempout%05i' %(N) , url]).wait()
    #if verbose > 2:
        #print ' '.join(["wget", "-q" ,"-O", 'tempout%05i' %(N) , url])
    #size = os.path.getsize('tempout%05i' %(N))

    #if size > 91.:
      #sub.Popen(['mv', 'tempout%05i' %(N), fitsname]).wait()
      #print fitsname
    #else:
      #sub.Popen(['rm', 'tempout%05i' %(N)]).wait()
  #return 

  
  
#def get_MAGES_cutout_from_file(ra,dec, fitsname, band='M1', imsize = 1./60., verbose=0):
    
  #if (not os.path.isfile(fitsname)) or clobber:
    #imsize = imsize*3600.  #in arcsec
    #if verbose > 2:
        #print 'cutout %s' %(fitsname)
    ##imsize in degrees

    #if band == 'M1':
        #bigfits = 'mages/mages_24_Emos0.5_BF4_tile2.fits'

    #head = pf.getheader(bigfits)

    #hdulist = pf.open(bigfits)
    ## Parse the WCS keywords in the primary HDU
    #wcs = pw.WCS(hdulist[0].header)

    ## Some pixel coordinates of interest.
    #skycrd = np.array([ra,dec])
    #skycrd = np.array([[ra,dec,0,0]], np.float_)
    #pixel = wcs.wcs_sky2pix(skycrd, 1)
    ## Some pixel coordinates of interest.

    #x = pixel[0][0]
    #y = pixel[0][1]
    ## calculate number of pixels to copy
    #pixsize = abs(wcs.wcs.cdelt[0])
    #if np.isnan(imsize):
        #imsize = 25.
    #N = 10.*(imsize/pixsize)
    #if verbose > 2:
        #print 'x=%.5f, y=%.5f, N=%i' %(x,y,N)

    #ximgsize = head.get('NAXIS1')
    #yimgsize = head.get('NAXIS2')

    #if x ==0:
        #x = ximgsize/2
    #if y ==0:
        #y = yimgsize/2

    #offcentre = False
    ## subimage limits: check if runs over edges
    #xlim1 =  x - (N/2)
    #if(xlim1<1):
        #xlim1=1
        #offcentre=True
        #xlim2 =  x + (N/2)
    #if(xlim2>ximgsize):
        #xlim2=ximgsize
        #offcentre=True
        #ylim1 =  y - (N/2)
    #if(ylim1<1):
        #ylim1=1
        #offcentre=True
        #ylim2 =  y + (N/2)
    #if(ylim2>yimgsize):
        #offcentre=True
        #ylim2=yimgsize

    #xl = int(xlim1)
    #yl = int(ylim1)
    #xu = int(xlim2)
    #yu = int(ylim2)
    #if verbose > 2:
        #print 'postage stamp is %i x %i pixels' %(xu-xl,yu-yl)

    ## make fits cutout
    #inps = bigfits + '[%0.0f:%0.0f,%0.0f:%0.0f]' %(xl,xu,yl,yu)
  
    #sub.Popen(["fitscopy", inps, fitsname]).wait()

    #return 

  
#def get_first_cutout(fitsname,ra,dec,imsize=3.0, verbose=0):
    #sra = ra_to_str( ra )
    #sdec = dec_to_str( dec )
    #url = "third.ucllnl.org/cgi-bin/firstimage?RA=%s %s&Dec=&Equinox=J2000&ImageSize=%.1f&MaxInt=10&FITS=1&Download=1" %(sra,sdec,imsize)
    ## get the first image
    #if not os.path.isfile(fitsname) or clobber:
        #print 'getting first image: ',url
        #sub.Popen(["wget","-q","-O",fitsname,url]).wait()
        #if verbose > 2:
            #print ' '.join(["wget","-q","-O",fitsname,url])
    #return
    
#def get_WSRT_cutout(ra,dec, fitsname, band='I', imsize=2., verbose=0, clobber=1):
    
    ##fitsname  = 'cut3050-5.fits.gz'
    ##'http://www.astron.nl/wow/zip/%s' %(fitsname)
  #print 'cutout %s' %(fitsname)
  ##imsize in arcmin
  #sra = ra_to_str( ra ).replace(':','+')
  #sdec = dec_to_str( dec ).replace(':','+')
  #url = "http://www.astron.nl/wow/testcode.php?POS=%s+%s%s&Equinox=J2000&SIZE=%.1f%s+%.1f&cutout=1&surv_id=5"    %(sra, r'%2B', sdec, imsize, r'%2C', imsize)
  
  ##url="http://www.noao.edu/ndwfs/cutout-form.dh5.php?ra=%s&dec=%s&rawidth=%.1f&decwidth=%.1f&filter=%s"    %(sra, sdec, imsize, imsize, band)
  
  ##url="http://www.noao.edu/ndwfs/cutout.dh5.php?ra=14:31:18.2&dec=35:15:50&rawidth=2.0&decwidth=2.0&s_filter_id=25"
  
  #if (not os.path.isfile(fitsname)) or clobber:
    #print '...downloading'
    #N = np.random.uniform(1,100000)
    #sub.Popen(["wget", "-q" ,"-O", 'tempout%05i' %(N) , url]).wait()
    #if verbose > 2:
        #print ' '.join(["wget", "-q" ,"-O", 'tempout%05i' %(N) , url])
    #size = os.path.getsize('tempout%05i' %(N))

    
    
    #t = open('tempout%05i' %(N),'r')
    #tt = t.readlines()
    #ttt = ' '.join(tt)
    
    ##print 'tempout%05i' %(N)
    ##cmd = ["cat", "tempout%05i"%(N), "|","grep" ,"'No data available'" ,"|" ,"wc", "-l" ]
    ##wc = sub.Popen(["cat", "tempout%05i"%(N), "|","grep" ,"'No data available'" ,"|" ,"wc", "-l" ], stdout=sub.PIPE, stderr=sub.PIPE)
    ##print cmd
    ##print wc
    ##if wc == 0:
        ##print 'no data'
        ##return 0
        
    
    ##os.system('cat tempout%05i |grep fits > templine%05i' %(N,N))
    
    ##t = open('templine%05i' %(N),'r')
    ##tt = t.readlines()
    ##ttt = tt[0]
    #i2 = ttt.find('No data available')
    #if i2 > 0:
        #print 'no coverage'
        #sub.Popen(['rm', 'tempout%05i' %(N)]).wait()
        #return 0
        
    #i2 = ttt.find('fits')
    #fitspath = ttt[i2-14:i2+7]
    #url = 'http://www.astron.nl/wow/%s' %(fitspath)
    #sub.Popen(["wget", "-q" ,"-O", 'tempout%05i' %(N) , url]).wait()
    #if verbose > 2:
        #print ' '.join(["wget", "-q" ,"-O", 'tempout%05i' %(N) , url])
    
    #sub.Popen(['mv', 'tempout%05i' %(N), fitsname+'.gz']).wait()
    #sub.Popen(['gunzip', fitsname+'.gz']).wait()
    #print fitsname
    #return 1
