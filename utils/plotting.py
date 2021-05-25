import numpy as np
import matplotlib as mpl
#mpl.use('Agg')
mpl.rc_file('~/.config/matplotlib/matplotlibrc')  # <-- the file containing your settings
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import os
import datetime

def fig_save_many(f, name, types=[".png"], dpi=200,metadata=None):
    if not isinstance(types, list):
        types = [types]
    for ext in types:
        if ext[0] != '.':
            ext = '.'+ext
        if 'png' in ext:
            metadata = {'Description': 'Created by '+ os.path.abspath((sys.argv[0])) +' at '+datetime.datetime.now().strftime("%c")}
            f.savefig(name+ext, dpi=dpi, metadata=metadata)
        else:
            f.savefig(name+ext, dpi=dpi)
    return

def paper_single(TW=6.64, AR=0.74, FF=1., fontsize=16.0, fontf='serif', fonts=["Times New Roman", "Computer Modern Roman", "STIXGeneral"], fontss=['Tahoma', 'DejaVu Sans', 'Lucida Grande', 'Verdana']):
    '''paper_single(TW = 6.64, AR = 0.74, FF = 1.)
    TW = 3.32
    AR = 0.74
    FF = 1.
    '''
    mpl.rcParams['font.family'] = fontf
    if fontf == 'sans-serif':
        mpl.rcParams['mathtext.fontset'] = 'dejavusans'
        mpl.rc('text', usetex=False) 
    else:
        mpl.rcParams['mathtext.fontset'] = 'cm'
        mpl.rc('text', usetex=True) 
    mpl.rcParams['font.sans-serif'] = fontss
    mpl.rcParams['font.serif'] = fonts
    mpl.rcParams['font.size'] = fontsize
    mpl.rc('figure', figsize=(FF*TW, FF*TW*AR), dpi=100)
    mpl.rc('figure.subplot', left=0.15, right=0.95, bottom=0.15, top=0.92)
    mpl.rc('lines', linewidth=1.75, markersize=8.0, markeredgewidth=0.75)
    #mpl.rc('font', size=fontsize, family=fontf, serif=fontst, sans-serif=fontsts)
    mpl.rc('xtick', labelsize='small')
    mpl.rc('ytick', labelsize='small')
    mpl.rc('xtick.major', width=1.0, size=8)
    mpl.rc('ytick.major', width=1.0, size=8)
    mpl.rc('xtick.minor', width=1.0, size=4)
    mpl.rc('ytick.minor', width=1.0, size=4)
    mpl.rc('axes', linewidth=1.5)
    mpl.rc('legend', fontsize='small', numpoints=1, labelspacing=0.4, frameon=False) 
    mpl.rc('savefig', dpi=300)
    
    
#def paper_single_mjh(TW = 8.0, AR = 0.75, FF = 1.):
    #'''paper_single(TW = 6.64, AR = 0.74, FF = 1.)
    #TW = 3.32
    #AR = 0.74
    #FF = 1.
    ##mpl.rc('figure', figsize=(4.5,3.34), dpi=200)
    #mpl.rc('figure', figsize=(FF*TW, FF*TW*AR), dpi=200)
    #mpl.rc('figure.subplot', left=0.18, right=0.97, bottom=0.18, top=0.9)
    #mpl.rc('lines', linewidth=1.0, markersize=4.0)
    #mpl.rc('font', size=9.0, family="serif", serif="CM")
    #mpl.rc('xtick', labelsize='small')
    #mpl.rc('ytick', labelsize='small')
    #mpl.rc('axes', linewidth=0.75)
    #mpl.rc('legend', fontsize='small', numpoints=1, labelspacing=0.4, frameon=False) 
    #mpl.rc('text', usetex=True) 
    #mpl.rc('savefig', dpi=300)
    #'''
    ##import matplotlib as mpl
    ## textwidht = 42pc = 42 * 12 pt = 42 * 12 * 1/72.27 inches
    ## columnsep = 2pc
    ## ... colwidth = 20pc
    
    ##mpl.rc('figure', figsize=(4.5,3.34), dpi=200)
    #mpl.rcdefaults()
    
    ##mpl.rc('figure', figsize=(FF*TW, FF*TW*AR), dpi=100)
    ##mpl.rc('figure.subplot', left=0.15, right=0.95, bottom=0.15, top=0.92)
    ##mpl.rc('lines', linewidth=1.75, markersize=8.0, markeredgewidth=0.75)
    #mpl.rc('font', size=20.0, family="serif", serif="Times")
    ##mpl.rc('xtick', labelsize='small')
    ##mpl.rc('ytick', labelsize='small')
    ##mpl.rc('xtick.major', width=1.0, size=8)
    ##mpl.rc('ytick.major', width=1.0, size=8)
    ##mpl.rc('xtick.minor', width=1.0, size=4)
    ##mpl.rc('ytick.minor', width=1.0, size=4)
    ##mpl.rc('axes', linewidth=1.5)
    ##mpl.rc('legend', fontsize='small', numpoints=1, labelspacing=0.4, frameon=False) 
    #mpl.rc('text', usetex=True) 
    #mpl.rc('savefig', dpi=300)
    
    
def paper_single_mult_ax(nrows=1, ncols=1, **kwargs):
    #import matplotlib as mpl
    paper_single(FF=max(nrows,ncols))
    f, ax = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
    plt.minorticks_on()
    ylocator6 = plt.MaxNLocator(5)
    xlocator6 = plt.MaxNLocator(6)
    if len(ax.shape) > 1:
        for axrow in ax:
            for axcol in axrow:
                axcol.xaxis.set_major_locator(xlocator6)
                axcol.yaxis.set_major_locator(ylocator6)
    else:
        for axcol in ax:
            axcol.xaxis.set_major_locator(xlocator6)
            axcol.yaxis.set_major_locator(ylocator6)
    return f, ax
    
    
def paper_single_ax(TW=6.64, AR=0.74, FF=1., fontsize=16.0,  fonts=["Times New Roman", "Computer Modern Roman", "STIXGeneral"], fontss=['Tahoma', 'DejaVu Sans', 'Lucida Grande', 'Verdana'], fontf='serif', projection=None):
    paper_single(TW=TW, AR=AR, FF=FF, fontsize=fontsize, fonts=fonts, fontss=fontss, fontf=fontf)
    f = plt.figure()
    ax = plt.subplot(111, projection=projection)
    plt.minorticks_on()
    ylocator6 = plt.MaxNLocator(5)
    xlocator6 = plt.MaxNLocator(6)
    ax.xaxis.set_major_locator(xlocator6)
    ax.yaxis.set_major_locator(ylocator6)
    return f, ax

def paper_double_ax():
    paper_single(TW = 12)
    f = plt.figure()
    ax = plt.subplot(111)
    plt.minorticks_on()
    ylocator6 = plt.MaxNLocator(5)
    xlocator6 = plt.MaxNLocator(6)
    ax.xaxis.set_major_locator(xlocator6)
    ax.yaxis.set_major_locator(ylocator6)
    return f, ax

def paper_double_mult_ax(nrows=1, ncols=1, setticks=True, FF=1, TW=6.97*2, AR=0.74, **kwargs):
    paper_single()
    
    mpl.rc('figure', figsize=(FF*TW, FF*TW*AR), dpi=200)
    mpl.rc('figure.subplot', left=0.1, right=0.97, bottom=0.1, top=0.97)
    mpl.rc('font', size=24.0, family="serif", serif="CM")
    
    f, ax = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
    plt.minorticks_on()
    if setticks:
        ylocator6 = plt.MaxNLocator(5)
        xlocator6 = plt.MaxNLocator(6)
        if len(ax.shape) > 1:
            for axrow in ax:
                for axcol in axrow:
                    axcol.xaxis.set_major_locator(xlocator6)
                    axcol.yaxis.set_major_locator(ylocator6)
        else:
            for axcol in ax:
                axcol.xaxis.set_major_locator(xlocator6)
                axcol.yaxis.set_major_locator(ylocator6)
    return f, ax


def plot_np_hist(ax, xbins, ydata, **kwargs):
    '''
    plot output of np.histogram with matplotlib
    ax - axes to plot into
    xbins - bin edges (dimension N+1)
    ydata - bin counts (dimension N)
    **kwargs - matplotlib.pyplot.plot arguments
    '''
    left,right = xbins[:-1],xbins[1:]
    X = np.array([left,right]).T.flatten()
    Y = np.array([ydata,ydata]).T.flatten()
    ax.plot(X,Y,**kwargs)
    return 

def add_colorbar(f,ax,c,label='',side='right',size='5%',pad=0.05):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(side, size=size, pad=pad)
    cbar=f.colorbar(c, cax=cax, orientation='vertical')
    if label != '':
        cbar.set_label(label)
    return


    

def invert_lim(ax,axis):
    assert axis.lower() in ['x','y'] , "Axis should be 'x' or 'y'"
    if axis.lower() =='x':
        invert_xlim(ax)
    else:
        invert_ylim(ax)
    return

def invert_xlim(ax):
    x1,x2 = ax.get_xlim()
    ax.set_xlim(x2,x1)
    return


def invert_ylim(ax):
    y1,y2 = ax.get_ylim()
    ax.set_ylim(y2,y1)
    return

def plot_equal(ax, col='k', **kwargs):
    '''
    plot one-to-one line
    '''
    x1,x2 = ax.get_xlim()
    y1,y2 = ax.get_ylim()
    
    l1 = min(x1,y1)
    l2 = max(x2,y2)
    
    ax.plot([l1,l2], [l1,l2],c=col, **kwargs)
    
    ax.set_xlim(x1,x2)
    ax.set_ylim(y1,y2)
    
    return

def format_log_axis10(ax, axis='both'):
    '''
    format a log axis as ... 0.1, 0, 1, 10 ...
    axis should be 'x','y',or 'both'
    default axis='both'
    '''
    assert axis in ['x','y','both'], 'invalid axis value'
    import matplotlib.ticker as ticker
    if axis in ['y','both']:
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
    if axis in ['x','both']:
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
    return

def set_attrib(ax, xlabel=None, ylabel=None, xlim=None, ylim=None, xtick_spacing=None, ytick_spacing=None, xtick_min_spacing=None, ytick_min_spacing=None, title=None):
    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    if xtick_spacing:
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(xtick_spacing))
    if xtick_min_spacing:
        ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(xtick_min_spacing))
    if ytick_spacing:
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(ytick_spacing))
    if ytick_min_spacing:
        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(ytick_min_spacing))
    return


def donley_mask(f_ch1, f_ch2, f_ch3, f_ch4, mags=True):
    '''Select sources using the Donley+ 2012 criteria
    returns  - array 1 where source in stern wedge, 0 where outside and -1 where any mag is nan
    '''
    
    if mags:
        f_ch1 = -2.5*np.log(f_ch1)
        f_ch2 = -2.5*np.log(f_ch2)
        f_ch3 = -2.5*np.log(f_ch3)
        f_ch4 = -2.5*np.log(f_ch4)
    
    x = np.log10(f_ch3/f_ch1)
    y = np.log10(f_ch4/f_ch2)
    
    s1 = (x >= 0.08)
    s2 = (y >= 0.15)
    s3 = (y >= 1.21*x - 0.27)
    s4 = (y <= 1.21*x + 0.27)
    s5 = (f_ch2 > f_ch1)
    s6 = (f_ch3 > f_ch2)
    s7 = (f_ch4 > f_ch3)
    nan1 = np.isnan(f_ch1)
    nan2 = np.isnan(f_ch2)
    nan3 = np.isnan(f_ch3)
    nan4 = np.isnan(f_ch4)
    nanmask = nan1 & nan2 & nan3 & nan4
    donleymask = s1 & s2 & s3 & s4 & s5 & s6 & s7 
    #maskind = np.where(mask)[0]
    donleymask[nanmask] = np.nan
    return donleymask

def stern_mask(m_ch1, m_ch2, m_ch3, m_ch4):
    '''Select sources in the Stern wedge
    returns  - array 1 where source in stern wedge, 0 where outside and -1 where any mag is nan
    '''
    s1 = ( (m_ch3 - m_ch4) > 0.6 )
    s2 = ( (m_ch1 - m_ch2) > 0.2*(m_ch3 - m_ch4) +0.18 )
    s3 = ( (m_ch1 - m_ch2) > 2.5*(m_ch3 - m_ch4) - 3.5 )
    nan1 = np.isnan(m_ch1)
    nan2 = np.isnan(m_ch2)
    nan3 = np.isnan(m_ch3)
    nan4 = np.isnan(m_ch4)
    nanmask = nan1 & nan2 & nan3 & nan4
    sternmask = s1 & s2 & s3
    #maskind = np.where(mask)[0]
    sternmask[nanmask] = np.nan
    return sternmask

def plot_stern_wedge(ax):
    # Stern+ 2005
    # [5.8]-[8.0] > 0.6 U [3.6]-[4.6] > 0.2([5.8]-[8.0]) +0.18 U ([3.6]-[4.5]) > 2.5([5.8]-[8.0]) - 3.5
    # intercepts (0.6, ymax) (0.6, 0.3) (1.6,0.5) (xmax, 2.5*xmax-3.5)  [0.6, -2.0]

    xmin,xmax = ax.get_xlim()
    ymin,ymax = ax.get_ylim()
    yd = 2.5*xmax -3.5
    if yd > ymax:
        xd = xmax
    else:
        xd = (ymax + 3.5)/2.5
        yd = ymax 
    ax.plot([0.6, 0.6, 1.6, xd], [ymax, 0.3, 0.5, yd], 'k')
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    return

def MIR_flux_to_mag(f1,f2,f3,f4,zp=25.):
    ## Vega mags  : add the extra (to go from AB to vega)
    ## 
    #m1 = -2.5*np.log10(f1) + 23.9 -2.788
    #m2 = -2.5*np.log10(f2) + 23.9 -3.255
    #m3 = -2.5*np.log10(f3) + 23.9 -3.743
    #m4 = -2.5*np.log10(f4) + 23.9 -4.372
    
    m1 = -2.5*np.log10(f1) + zp -2.788
    m2 = -2.5*np.log10(f2) + zp -3.255
    m3 = -2.5*np.log10(f3) + zp -3.743
    m4 = -2.5*np.log10(f4) + zp -4.372
    return [m1,m2,m3,m4]


def MIR_AGN(m1,m2,m3,m4, mode='flux'):
    if mode == 'flux':
        ## Vega mags  : add the extra
        m1,m2,m3,m4 = MIR_flux_to_mag(m1,m2,m3,m4)
    agn_mask = stern_mask(m1, m2, m3, m4)
    return agn_mask

def plot_mir_colours(m1,m2,m3,m4, ax=None, spec_ind=None, col=None, vmin=None, vmax=None):
    
    if ax is None:
        f, ax = paper_single_ax()
    
    ax.minorticks_on()
    
    A = 0.5
    if len(m1) > 1000:
        A = 0.1
    elif len(m1) > 10000:
        A = 0.01
        
    if col is None:
        ax.plot(m3 -m4,  m1-m2, 'k.', alpha=A)
    else:
        if vmin is None:
            vmin = min(col)
        if vmax is None:
            vmax = max(col)
        c=ax.scatter(m3 -m4,  m1-m2, c=col, marker='o', edgecolor='none', vmin=vmin, vmax=vmax)
        try: plt.colorbar(c)
        except: print("error")
    
    if spec_ind is not None:
        if col is None:
            ax.plot((m3-m4)[spec_ind],  (m1-m2)[spec_ind], 'rx', alpha=A)
        else:
            ax.plot((m3-m4)[spec_ind],  (m1-m2)[spec_ind], 'rx', alpha=A)
        
    
    plot_stern_wedge(ax)
    
    ax.set_xlabel('[5.8] - [8.0]')
    ax.set_ylabel('[3.6] - [4.5]')
    
    #ax.set_xlim(-5,5)
    #ax.set_ylim(-0.4, 1.2)
    
    #plt.savefig(savename)
    
    return f,ax


def hist2d(ax, xdat, ydat, xyrange, bins, thresh=2, cmap=plt.cm.Greys, log=False, scatterother=False):
    import scipy

    tt = ax.get_aspect()

    # histogram the data
    hh, locx, locy = scipy.histogram2d(xdat, ydat, range=xyrange, bins=bins)
    mhh = np.mean(hh)
    shh = np.std(hh)
    if log:
        lhh = np.log10(hh)
    else:
        lhh = hh
    posx = np.digitize(xdat, locx)
    posy = np.digitize(ydat, locy)


    #select points within the histogram
    ind = (posx > 0) & (posx <= bins[0]) & (posy > 0) & (posy <= bins[1])
    hhsub = hh[posx[ind] - 1, posy[ind] - 1] # values of the histogram where the points are
    xdat1 = xdat[ind][hhsub < thresh] # low density points
    ydat1 = ydat[ind][hhsub < thresh]
    lhh[hh  < thresh] = np.nan # fill the areas with low density by NaNs

    ar = (0.6/0.65)*(np.diff(xyrange[0])/np.diff(xyrange[1]))[0]
    c = ax.imshow(np.flipud(lhh.T),extent=np.array(xyrange).flatten(), interpolation='none', cmap=cmap, aspect=ar)  
    
    ax.set_aspect(tt)
    
    if scatterother:
        ax.plot(xdat1, ydat1, 'k,')    
    
    
    return c


def make_ax3():
    paper_single(TW=8, AR=0.9)
    f = plt.figure()
    
    from matplotlib.ticker import NullFormatter, MaxNLocator

    nullfmt   = NullFormatter()         # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.6
    bottom_h = bottom+height+0.02
    left_h = left+width+0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    ax = plt.axes(rect_scatter)
    plt.minorticks_on()
    axx = plt.axes(rect_histx)
    plt.minorticks_on()
    axy = plt.axes(rect_histy)
    plt.minorticks_on()

    # no labels
    axx.xaxis.set_major_formatter(nullfmt)
    axy.yaxis.set_major_formatter(nullfmt)
    
    
    axy.xaxis.set_major_locator(MaxNLocator(3))
    axx.yaxis.set_major_locator(MaxNLocator(3))
    
    return f,ax,axx,axy



def plot_opt_radio_overlay(opt_cutout, radio_cutout, outname, markers=None):
    
    import astropy.io.fits as fits
    from astropy.stats import sigma_clipped_stats
    import aplpy as ap
    
    def flatten(fitsin):
        fh = fits.open(fitsin)
        data = fh[0].data.squeeze() # drops the size-1 axes
        header = fh[0].header
        mywcs = wcs.WCS(header).celestial
        new_header = mywcs.to_header()
        new_fh = fits.PrimaryHDU(data=data, header=new_header)
        new_fh.writeto(fitsin, clobber=True)
        return

    def  get_noise(fitsfile, vm=False):
        dat = fits.getdata(fitsfile)
        if np.sum(np.isfinite(dat))==0:
            return np.nan, np.nan, np.nan
        mean, median, stddev = sigma_clipped_stats(dat, sigma=5, iters=10)
        if not vm:
            return mean, median, stddev
        else:
            vmax=np.median(dat[dat>(mean+5*stddev)])
            return mean, median, stddev, vmax

    
    
    f = plt.figure()
    #ax1 = f.add_subplot(121)
    ##ax2 = f.add_subplot(122)
    
    
    m1, m2, rms, vm = get_noise(opt_cutout, vm=True)
    
    ax0 = ap.FITSFigure(opt_cutout,figure=f,north=True)
    ax0.show_colorscale(vmin=m1+1*rms, vmax=vm, stretch='log')


    
    m1, m2, rms, vm = get_noise(radio_cutout, vm=True)
    radiomax=np.nanmax(fits.getdata(radio_cutout))
    drlimit=2000
    #print radiomax/drlimit,rms*2.0
    minlevel=max([radiomax/drlimit,rms*2.0])
    levels=minlevel*2.0**np.linspace(0,14,30)
    ax0.show_contour(radio_cutout, levels=levels, colors='y')
    

    
    ax0.hide_axis_labels()
    ax0.hide_tick_labels()
    
    ax0.add_scalebar(1/60.)  # size in deg
    ax0.scalebar.set_color('w')
    ax0.scalebar.set_corner('top right')
    ax0.scalebar.set_label("$1'$")  # length in degrees
    ax0.scalebar.set_font(size='small')
    ax0.ticks.hide()
    
    if markers is not None:
        ax0.show_markers(markers['ra'],markers['dec'])
    
    f.savefig(outname)
    
    plt.close()
    


    return


def nanhist(x,**kwargs):
    
    n,b,p = plt.hist(x[np.isfinite(x)],**kwargs)
    return n,b,p

#xmin,xmax = ax.get_xlim()
#ymin,ymax = ax.get_ylim()
#ax.plot([0.6, 0.6, 1.6, xmax], [ymax, 0.3, 0.5, 2.5*xmax -3.5], 'k')


markers=['o','s','D','p','H','v','^','<','>','1','2','3','4','h']
mylinestyles = ['-', '--','-.', ':' ]
mylinestyles = ['-', '-','-', '-' ,'-','-','-']

#sample_colours = ['#ff0000', '#ff0e00', '#ff1b00', '#ff2900', '#ff3700', '#ff4500', '#ff5200', '#ff6000', '#ff6e00', '#ff7c00', '#ff8a00', '#ff9700' ', ''#ffa500']
#sample_colours = ['#ff0000', '#ff5200', '#ffa500']
erg_colours = ['c','m']
sample_colours = ['#fecc5c','#fd8d3c', '#e31a1c','darkred']
sample_markers = ['s', '^', 'o']
sample_colours = ['#fecc5c','#fd8d3c', '#e31a1c','darkred']

def sample_color_range(Ncol):
    return plt.cm.YlOrRd(np.linspace(0, 1, Ncol))
sample_markers = ['s', '^', 'o','D','d','h']


pcut_colours = ['#e31a1c', '#fd8d3c', '#fecc5c']

#massbin_colours = ['#eaad15', '#bfbb40', '#95ca6a', '#6ada95', '#40e8bf','#15f7ea']
#massbin_colours = [ '#6ada95','#95ca6a', '#bfbb40', '#eaad15']
#massbin_colours = [ '#6ada95','#95ca6a', '#eaad15']
#massbin_colours = [ '#6ada95','#eaad15', '#eaad15']
#massbin_colours = ['#5ab4ac','#d8b365']
# http://colorbrewer2.org/
massbin_colours = ['#018571','#80cdc1','#dfc27d','#a6611a']
massbin_markers = ['s', '^','o','D']
#massbin_colours = [ '#6ada95', '#eaad15']
 
#massbin_colours = ['#005d61','#006155','#005d61','#004d61','#004561']
#massbin_colours = ['#005d61','#00614d','#006155','#00615d','#005d61','#005561','#004d61','#004561']


