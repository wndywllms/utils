import numpy as np

def stern_mask(m_ch1, m_ch2, m_ch3, m_ch4):
    '''Select sources in the Stern wedge
    returns  - array 1 where source in stern wedge, 0 where outside and -1 where any mag is nan
    '''
    s1 = ( (m_ch3 - m_ch4) > 0.6 )
    s2 = ( (m_ch1 - m_ch2) > 0.2*(m_ch3 - m_ch4) +0.18 )
    s3 = ( (m_ch1 - m_ch2) > 2.5*(m_ch3 - m_ch4) - 3.5 )
    #nan1 = np.isnan(m_ch1)
    #nan2 = np.isnan(m_ch2)
    #nan3 = np.isnan(m_ch3)
    #nan4 = np.isnan(m_ch4)
    #nanmask = nan1 & nan2 & nan3 & nan4
    nan1 = np.isfinite(m_ch1)
    nan2 = np.isfinite(m_ch2)
    nan3 = np.isfinite(m_ch3)
    nan4 = np.isfinite(m_ch4)
    finmask = nan1 & nan2 & nan3 & nan4
    sternmask = s1 & s2 & s3 & finmask
    #maskind = np.where(mask)[0]
    #m=1*sternmask -1*nanmask
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

def MIR_flux_to_mag(f1,f2,f3,f4):
    ## Vega mags  : add the extra
    m1 = -2.5*np.log10(f1) + 23.9 -2.788
    m2 = -2.5*np.log10(f2) + 23.9 -3.255
    m3 = -2.5*np.log10(f3) + 23.9 -3.743
    m4 = -2.5*np.log10(f4) + 23.9 -4.372
    return [m1,m2,m3,m4]


def MIR_AGN(m1,m2,m3,m4, mode='flux'):
    if mode == 'flux':
        ## Vega mags  : add the extra
        m1,m2,m3,m4 = MIR_flux_to_mag(m1,m2,m3,m4)
    agn_mask = stern_mask(m1, m2, m3, m4)
    return agn_mask
