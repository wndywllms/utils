#!/usr/bin/python

import matplotlib as mpl

import os
import sys
from pylab import *
from scipy import *
from numpy import *
import argparse



# define column names and formats and default values
DEFAULTS = ['', '', '', '', '', 0.0, 0.0, 0.0, 0.0, 0.0, np.zeros((4,)), 0.0, 0.0, 0.0]
NAMES = ['Name', 'Type', 'Patch', 'Ra', 'Dec', 'I', 'Q', 'U', 'V', 'ReferenceFrequency', 'SpectralIndex', 'MajorAxis', 'MinorAxis', 'Orientation']
FORMATS = ['S100', 'S10', 'S100', 'S15', 'S15', 'd', 'd', 'd', 'd', 'd', '(4,)d', 'd', 'd', 'd']

############################################################################

def ra_to_degrees(ra_str, delim=':'):
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

def dec_to_degrees(dec_str, delim='.'):
    '''
    converts array of strings or single string dec values to decimal degrees
    '''
    # string or array
    if isinstance(dec_str,str):
        if delim == '.':
            #dec_str[dec_str.find('.')] = ':'  # replace first
            #dec_str[dec_str.find('.')] = ':'  # replace second
            dec_str = dec_str.replace('.',':',2)
            t = dec_str.split(':')
        else:
            t = dec_str.split(delim)
        dec_deg = (float(t[0]) + float(t[1])/60. + float(t[2])/3600.)
        return dec_deg
    else:
        dec_deg = np.zeros(len(dec_str))
        for i,dec_s in enumerate(dec_str):
            if delim == '.':
                #dec_str[dec_str.find('.')] = ':'  # replace first
                #dec_str[dec_str.find('.')] = ':'  # replace second
                dec_s = dec_s.replace('.',':',2)
                t = dec_s.split(':')
            else:
                t = dec_s.split(delim)
            dec_deg[i] = (float(t[0]) + float(t[1])/60. + float(t[2])/3600.)
        return dec_deg        

def ra_to_str(dra, delim=':'):
    '''
    converts a single decimal degrees ra to hh:mm:ss.s
    '''
    if isinstance(dra,float):
        dra = dra/15.
        dd = math.floor(dra)
        dfrac = dra - dd
        dmins = dfrac*60.
        dm = math.floor(dmins)
        dsec = (dmins-dm)*60.
        sra = '{dd:02.0f}{delim}{dm:02.0f}{delim}{dsec:04.1f}'.format(dd=dd,dm=dm,dsec=dsec, delim=delim)  
    else:
        sra = np.array(len(dra),dtype='string')
        for i,ra in enumerate(dra):
            ddra = ra/15.
            dd = math.floor(ra)
            dfrac = ddra - dd
            dmins = dfrac*60.
            dm = math.floor(dmins)
            dsec = (dmins-dm)*60.
            sra[i] = '{dd:02.0f}{delim}{dm:02.0f}{delim}{dsec:04.1f}'.format(dd=dd,dm=dm,dsec=dsec, delim=delim)  
        
    return sra
def dec_to_str(ddec, delim='.'):
    '''
    converts a single decimal degrees dec to dd:mm:ss.s
    '''
    if isinstance(ddec,float):
        dd = math.floor(ddec)
        dfrac = ddec - dd
        dmins = dfrac*60.
        dm = math.floor(dmins)
        dsec = (dmins-dm)*60.
        sdec = '{dd:02.0f}{delim}{dm:02.0f}{delim}{dsec:04.1f}'.format(dd=dd,dm=dm,dsec=dsec, delim=delim)  
    else:
        sdec = np.array(len(ddec),dtype='string')
        for i,dec in enumerate(ddec):
            dd = math.floor(dec)
            dfrac = dec - dd
            dmins = dfrac*60.
            dm = math.floor(dmins)
            dsec = (dmins-dm)*60.
            sdec[i] = '{dd:02.0f}{delim}{dm:02.0f}{delim}{dsec:04.1f}'.format(dd=dd,dm=dm,dsec=dsec, delim=delim)  
    return sdec

     

def read_bbs(infile, verbose=False):
    print '* reading input bbs skymodel: %s' %infile
    
    with open(infile,'r') as t:
        # read the first line as format line
        line1 = t.readline()
        t2 = line1.strip()
        i1 = t2.find('(')
        i2 = t2.find(')')
        if i1 > 0:
            t2 = t2[i1+1:i2]
        elif 'FORMAT = ' in t2:       
            t2 = t2.replace('FORMAT = ','')
        elif 'format = ' in t2:       
            t2 = t2.replace('format = ','')
        t2 = t2.replace(' ','')  
        # deal with the specind possibility
        i1 = t2.find('[')
        i2 = t2.find(']')
        t2seg = t2[0:i1].replace(',',';') + t2[i1:i2] + t2[i2:].replace(',',';')
        innames = t2seg.split(';')
        Ncol = len(innames)    
        #writefile.write("# "+line1)
        
        # store what formats and defaults are given
        informats = []
        indefaults = []
        for i in range(len(innames)):
            # check if default value
            if '=' in innames[i]:
                innamesplit = innames[i].split('=')
                innames[i] = innamesplit[0]
                defval = innamesplit[1].replace("'","")
                # case of spec ind
                if '[' in defval:
                    defvals = defval.replace('[','').replace(']','').split(',')
                    for idv, defval in enumerate(defvals):
                        if defval != '':
                            DEFAULTS[NAMES.index(innamesplit[0])][idv] = float(defval)
                else:
                    # update the default, check for the case of '[]'
                    if defval != '[]':
                        if isinstance(DEFAULTS[NAMES.index(innamesplit[0])], float) : DEFAULTS[NAMES.index(innamesplit[0])] =float(defval)
                        else: DEFAULTS[NAMES.index(innamesplit[0])] = defval
                #informats.append('S1')
            # add to list
            informats.append(FORMATS[NAMES.index(innames[i])])    
            indefaults.append(DEFAULTS[NAMES.index(innames[i])])
            
        # Now get the values    
        textarray = []  
        #textarray_patches = []
        t = open(infile,'r')
        # read rest of lines
        t2 = t.readlines()[1:]
        t.close()
        for line in t2:
            #blank lines
            if line=='\n': continue
            # comments are '#'
            if '#' in line: continue
            line = line.strip()
            C = line.strip().split(', ')
            # patch lines
            #if len(C) == 5:
                #textarray_patches.append(C)
            #else:
            if 1:
                # deal with the specind
                i1 = line.find('[')
                i2 = line.find(']')
                lineseg = line[0:i1].replace(',',';') + line[i1:i2] + line[i2:].replace(',',';')
                C = lineseg.split('; ')
                textarray.append(C)
                
        textarray = np.array(textarray)
        nsources = textarray.shape[0]
        # build the bbs record arrya
        dt = np.dtype(dict(names=NAMES,formats=FORMATS))
        bbsdata = np.recarray((nsources,), dtype=dt)
            
        # store all defualts
        for col in bbsdata.dtype.names:
            colind = NAMES.index(col)
            for i in range(nsources):
                bbsdata[col][i] = DEFAULTS[colind]
        
        # move the real values in
        for name in innames:
            if name == 'SpectralIndex':
                for n in range(nsources):      
                    ni = innames.index(name)
                    if ni < len(textarray[n]):
                        val = textarray[n][ni]
                        vals = val.replace('[','').replace(']','').split(',')
                        for iv, val in enumerate(vals):
                            bbsdata[name][n][iv] = float(val)
            else:
                for n in range(nsources):
                    i = innames.index(name)
                    if i >= len(textarray[n]): continue
                    newval = textarray[n][i]
                    if len(newval) == 0: continue
                    if FORMATS[i] == 'd': bbsdata[name][n] = float(newval)
                    else: bbsdata[name][n] = newval
                    
    #bbsdata
    data = bbsdata[np.where(bbsdata.Name != '')]
    patchdata = bbsdata[np.where(bbsdata.Name == '')]
    
    # if there are no patches, put all sources in their own patch
    if len(patchdata) == 0:
        
        patchdata = np.recarray((nsources,), dtype=dt)
        # store all defualts
        for col in patchdata.dtype.names:
            colind = NAMES.index(col)
            for i in range(nsources):
                patchdata[col][i] = DEFAULTS[colind]
        #patchdata = bbsdata[np.where(bbsdata.Name != '')]
        # add the names of the patches 
        for iPatch, ThisPatch in enumerate(patchdata):
            patchdata.Patch[iPatch] = data.Name[iPatch]
            data.Patch[iPatch] = data.Name[iPatch]
        
    # now we definitely have patches - update their positions and fluxes
        
    for iPatch, ThisPatch in enumerate(patchdata):
        iPatchSources = np.where(data.Patch == ThisPatch.Patch)
        patchdata.I[iPatch] = np.sum(data.I[iPatchSources])
        #print bbsdata['Ra'][iPatchSources]
        tra = ra_to_degrees(data.Ra[iPatchSources])
        tdec = dec_to_degrees(data.Dec[iPatchSources])
        patchdata.Ra[iPatch] = ra_to_str(np.average(tra, weights=data.I[iPatchSources]))
        patchdata.Dec[iPatch] = dec_to_str(np.average(tdec, weights=data.I[iPatchSources]))
   
    if verbose:
        print 'BBS recarray: ', data
        print 'BBS patches recarray: ', patchdata
    
        
    return data, patchdata


def sort_by_flux(bbs_data, bbsdata_patches):
    ### NOTE should always sort by patch flux, because even if no patches are specified in the input bbs file, we automatically put each source in its own patch
    print '* sorting by I flux'
    if len(bbsdata_patches) > 0:
        bbsdata_patches.sort(order='I')
        return bbs_data, bbsdata_patches[::-1]
    else:
        bbs_data.sort(order='I')
        return bbs_data[::-1], bbsdata_patches

def write_bbsdata_to_file(bbsdata, bbsdata_patches, outbbs):
    print '* writing output bbs skymodel: %s ' %outbbs
    with open(outbbs, 'w') as f:
        write_columns = ['Name', 'Type', 'Patch', 'Ra', 'Dec', 'I', 'Q', 'U', 'V', 'ReferenceFrequency', 'SpectralIndex', 'MajorAxis', 'MinorAxis', 'Orientation']
        fmtline = "# (%s) = format" %(', '.join(write_columns))
        f.write(fmtline+'\n')
        # if there are patches
        if len(bbsdata_patches) > 0:
            for ipatch in range(len(bbsdata_patches)):
                patchname = bbsdata_patches.Patch[ipatch]
                patchline = " , , %s, %s, %s, %f" %(patchname, bbsdata_patches.Ra[ipatch], bbsdata_patches.Dec[ipatch], bbsdata_patches.I[ipatch])
                f.write(patchline+'\n')
                isources = np.where(bbsdata.Patch == patchname)[0]
                for i in isources:
                    source_columns = []
                    for col in write_columns:
                        t = bbsdata[col][i]
                        if isinstance(t, np.ndarray):
                            strarray = '[{s}]'.format(s=', '.join(['{tt:3f}'.format(tt=tt) for tt in t]))
                            source_columns.append(strarray)
                        else: source_columns.append(str(t))
                    f.write(', '.join(source_columns)+'\n')
            # are there any left out sources (i.e without patch specified... this should only happen if there is one patch in the file...
            for i in range(len(bbsdata)):
                if bbsdata['Patch'][i] not in bbsdata_patches.Patch:
                    source_columns = []
                    for col in write_columns:
                        t = bbsdata[col][i]
                        if isinstance(t, np.ndarray): source_columns.append(str(list(t)))
                        else: source_columns.append(str(t))
                    f.write(', '.join(source_columns)+'\n')
                    
        # no patches, just write sources
        # this should never happen!!
        else:
            print "why do you have no patches?! there is a problem here..."
            print " writing only the sources"
            for i in range(len(bbsdata)):
                source_columns = []
                for col in write_columns:
                    t = bbsdata[col][i]
                    if isinstance(t, np.ndarray): source_columns.append(str(list(t)))
                    else: source_columns.append(str(t))
                f.write(', '.join(source_columns)+'\n')
        
    return


def write_region_file(bbsdata, bbsdata_patches, reg, gcol='blue', pcol='green', col='yellow'):
    '''write a region file for ds9 based on sources/patches'''
    print '* writing output file: %s' %reg
    #region file
    if os.path.isfile(reg):
        inp = raw_input('region file %s exists, overwrite (y)? ' %(reg))
        if inp in ['n','N','no','No']: return
    os.system('rm -f '+reg)

    #Create ds9 region file
    regfile=open(reg,'w')
    regfile.write('# Region file format: DS9 version 4.1 \n')
    regfile.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n')
    regfile.write('fk5 \n')


    #Create ds9 region file
    regfile=open(reg,'w')
    #write the sources
    for i in range(len(bbsdata)):
        ra   = bbsdata.Ra[i]
        dec  = bbsdata.Dec[i]
        src  = bbsdata.Name[i]
        dec  = dec.replace(".",":",2)
        typ = bbsdata.Type[i]

        if typ == 'GAUSSIAN':
            a = bbsdata.MajorAxis[i]
            b = bbsdata.MinorAxis[i]
            pa = bbsdata.Orientation[i]
            color = gcol
            regfile.write('ellipse(%s,%s,%f",%f",%f) # color=%s text={%s} \n' %(ra, dec, a, b, pa, color, src))

        elif typ == 'POINT':
            color = pcol
            regfile.write('point(%s,%s) # point=circle color=%s text={%s} \n' %(ra, dec, color, src))

        else:
            color = col
            regfile.write('point(%s,%s) # point=circle color=%s text={%s} \n' %(ra, dec, color, src))
    # write the patches
    for i in range(len(bbsdata_patches)):
        ra   = bbsdata_patches.Ra[i]
        dec  = bbsdata_patches.Dec[i]
        src  = bbsdata_patches.Patch[i]
        dec  = dec.replace(".",":",2)

        color = col
        regfile.write('point('+ra+','+dec+') # point=circle color='+color+' text={'+src+'} \n')

    regfile.close()
    
    return

"""
def put_sources_in_patches(bbsdata, bbsdata_patches, patch_root='Field'):
    '''puts each source in its own patch
    '''
    print '* putting sources in patches'
    nsources = bbsdata.shape[0]
    textarray_patches = []
    for i in range(nsources):
        # it has no patch, make a patch for it
        if bbsdata.Patch[i] == '':
            bbsdata.Patch[i] = patch_root + str(i)
            textarray_patches.append([bbsdata.Patch[i], bbsdata.Ra[i], bbsdata.Dec[i]  ])
        # is it an existing patch, then do not rename it but store the patch
        else:
            ipatch = np.where(bbsdata_patches.Patch == bbsdata.Patch[i])[0]
            textarray_patches.append([bbsdata.Patch[i], bbsdata_patches.Ra[ipatch], bbsdata_patches.Dec[ipatch]  ])
            
    textarray_patches = np.array(textarray_patches)
    npatches = textarray_patches.shape[0]
    
    # build the patch recarray from scratch
    names_patches = ['Patch','Ra','Dec']
    formats_patches = ['S30','S30','S30']
    dt = dtype(dict(names=names_patches, formats=formats_patches))
    bbsdata_patches = np.recarray((npatches,), dtype=dt)
    if npatches > 0:
        patch0 = textarray_patches[0]
        for ip,item in enumerate(patch0):
            if len(item) > 1 :
                if ':' in item: raind = ip
                elif '.' in item: decind = ip
                else: patchind = ip
        for i in range(npatches):
            bbsdata_patches['Ra'][i] = textarray_patches[i][raind]
            bbsdata_patches['Dec'][i] = textarray_patches[i][decind]
            bbsdata_patches['Patch'][i] = textarray_patches[i][patchind]
                
    return bbsdata, bbsdata_patches
"""

def angular_separation(ra1,dec1,RA,DEC):
    separation = np.sqrt((np.cos(DEC*np.pi/180.))**2.*(ra1-RA)**2.  + (dec1-DEC)**2.)
    return separation

def rename_patch(bbsdata, bbsdata_patches, ra_c, dec_c, rad = 10./3600., patchroot='NEWPATCH'):
    '''renames all patches within rad of position ra_c, dec_c'''
    
    if isinstance(ra_c, str):
        ra_c = ra_to_degrees(ra_c)
    if isinstance(dec_c, str):
        dec_c = dec_to_degrees(dec_c)
    
    #Ra = ra_to_degrees(bbsdata.Ra)
    #Dec = dec_to_degrees(bbsdata.Dec)
    #sep = angular_separation(Ra, Dec, ra_c, dec_c)
    #ind = np.where(sep<rad)[0]
    
    patchRa = ra_to_degrees(bbsdata_patches.Ra)
    patchDec = dec_to_degrees(bbsdata_patches.Dec)
    patchsep = angular_separation(patchRa, patchDec, ra_c, dec_c)
    patchind = np.where(patchsep<rad)[0]
    otherpatchind = np.where(patchsep>=rad)[0]
    print ra_c, dec_c, patchind
    if patchroot == '':
        patchroot = bbsdata_patches.Patch[patchind[0]]
    
    if len(patchind) > 1:
        print 'more than 1 patch being combined'
    else:
        print "only 1 patch"
    for pi in patchind:
        oldpatch = bbsdata_patches.Patch[pi]
        bbsdata_patches.Patch[pi] = patchroot
        # get the sources and rename their patches
        spind = np.where(bbsdata.Patch == oldpatch)[0]
        for spi in spind:
            bbsdata.Patch[spind] = patchroot
    
    # keep only one of the patches that now have the same name
    print patchind
    selectpatches = np.append(otherpatchind, patchind[0])
    bbsdata_patches = bbsdata_patches[selectpatches]
    
    # update the patch info
    iPatch = np.where(bbsdata_patches.Patch == patchroot)[0]
    if len(iPatch) > 1:
        print "we have some problem: too many patches with the new name"
    ThisPatch = bbsdata_patches[iPatch]
    iPatchSources = np.where(bbsdata.Patch == ThisPatch.Patch)
    bbsdata_patches.I[iPatch] = np.sum(bbsdata.I[iPatchSources])
    #print bbsdata['Ra'][iPatchSources]
    tra = ra_to_degrees(bbsdata.Ra[iPatchSources])
    tdec = dec_to_degrees(bbsdata.Dec[iPatchSources])
    bbsdata_patches.Ra[iPatch] = ra_to_str(np.average(tra, weights=bbsdata.I[iPatchSources]))
    bbsdata_patches.Dec[iPatch] = dec_to_str(np.average(tdec, weights=bbsdata.I[iPatchSources]))
        
    ## it assumes all patches are unique: TODO
    #for i in ind:
        #ipatch = np.where(bbsdata_patches.Patch==bbsdata.Patch[i])
        #if verbose > 2: print bbsdata.Patch[i]
        #bbsdata.Patch[i] = patchroot
        #bbsdata_patches.Patch[ipatch] = bbsdata.Patch[i]
        #if verbose > 2: print bbsdata.Patch[i], ipatch, bbsdata_patches[ipatch]
        
    ##update the patchlist
    #keep_patches = []
    #patches = bbsdata_patches.Patch
    #for patch in bbsdata.Patch:
        #if patch in patches:
            #patches_ind = np.where(patch==patches)[0]
            ##if verbose > 2: print patches_ind
            #keep_patches.append(patches_ind[0])
    #bbsdata_patches = bbsdata_patches[keep_patches]
    
    return bbsdata, bbsdata_patches

def remove_sources(bbsdata, bbsdata_patches, ra_c, dec_c, rad = 10./3600.):
    '''removes all sources (and their patches) within rad of position ra_c, dec_c
    [ this can be used for combining two gsm.py outputs... to include brighter sources further out]
    e.g. > gsm.py 4C28.37.gsmsky5deg.model 221.68041667 27.95027778 5. 3.
    > gsm.py 4C28.37.gsmsky3deg.model 221.68041667 27.95027778 3. 1.
    then remove_sources < 3deg in gsmsky5deg.model and put in the gsmsky3deg.model'''
    Ra = ra_to_degrees(bbsdata.Ra)
    Dec = dec_to_degrees(bbsdata.Dec)
    sep = angular_separation(Ra, Dec, ra_c, dec_c)
    ind = np.where(sep>=rad)[0]
    # keep the sources outside the radius
    bbsdata = bbsdata[ind]
    
    # keep only the patches that still exist
    patches = np.unique(bbsdata.Patch)
    ind_incl = []
    # store the indices of patches still in the source list
    for i, patch in enumerate(bbsdata_patches.Patch):
        if patch in patches: ind_incl.append(i)
    bbsdata_patches = bbsdata_patches[ind_incl]
    
    return bbsdata, bbsdata_patches

def plot_sky(bbsdata, bbsdata_patches, ra_c, dec_c, figname='bbs.skymodel.png', scale='log'):
    print '* plotting image: %s' %(figname)
    import pylab as pl
    fig = pl.figure()
    ax1 = pl.subplot(111)
    sources = bbsdata.Name
    Ra = ra_to_degrees(bbsdata.Ra)
    Dec = dec_to_degrees(bbsdata.Dec)
    I = bbsdata.I
    Ra_p = ra_to_degrees(bbsdata_patches.Ra)
    Dec_p = dec_to_degrees(bbsdata_patches.Dec)
    
    if len(Ra_p) > 0:
        if np.sum(Ra_p) > 0:
            ax1.scatter(Ra_p, Dec_p, s=40, c='gray', marker='x')
    
    if scale == 'log':
        c = ax1.scatter(Ra, Dec, c=np.log10(I), edgecolor='none')
    else:
        c = ax1.scatter(Ra, Dec, c=I, edgecolor='none')
    x1,x2 = ax1.get_xlim()
    y1,y2 = ax1.get_ylim()
    cbar = fig.colorbar(c)
    cbar.set_label('log I')
        
    t = np.arange(0, 2.*np.pi+0.001, np.pi/360.)
    circs = [1., 2., 5.] ## degrees
    if not np.isnan(ra_c): RA = ra_c
    else: RA = np.median(Ra)
    if not np.isnan(dec_c): DEC = dec_c
    else: DEC = np.median(Dec)
    for R in circs:
        ax1.plot(RA+R*np.cos(t), DEC+R*np.sin(t), linestyle='dotted', color='gray' )
    
    ax1.set_xlim(x2,x1)
    ax1.set_ylim(y1,y2)    
    ax1.set_xlabel('Ra [deg]')
    ax1.set_ylabel('Dec [deg]')
    pl.savefig(figname)
    return

def append_sky_models(bbsdata1, bbsdata_patches1, bbsdata2, bbsdata_patches2):
    nsources1 = bbsdata1.shape[0]
    nsources2 = bbsdata2.shape[0]
    
    # build the bbs record arrya
    dt = dtype(dict(names=NAMES,formats=FORMATS))
    bbsdata = np.recarray((nsources1+nsources2,), dtype=dt)
        
    # store all defualts
    for col in bbsdata.dtype.names:
        colind = NAMES.index(col)
        for i in range(nsources1+nsources2):
            bbsdata[col][i] = DEFAULTS[colind]
    
    # store all values
    for i in range(nsources1):
        bbsdata[i] = bbsdata1[i]
    for i in range(nsources2):
        bbsdata[nsources1+i] = bbsdata2[i]
        
    
    npatches1 = bbsdata_patches1.shape[0]
    npatches2 = bbsdata_patches2.shape[0]
    # build the patch recarray
    names_patches = ['Patch','Ra','Dec']
    formats_patches = ['S30','S30','S30']
    dt = dtype(dict(names=names_patches, formats=formats_patches))
    bbsdata_patches = np.recarray((npatches1+npatches2,), dtype=dt)
    # store all values
    for i in range(npatches1):
        bbsdata_patches[i] = bbsdata_patches1[i]
    for i in range(npatches2):
        bbsdata_patches[npatches1+i] = bbsdata_patches2[i]

        
    return bbsdata, bbsdata_patches

def apply_fluxcut(bbsdat, bbsdat_patches, fcut):
    selind = np.where(bbsdat_patches.I > fcut)[0]
    bbsdat_patches = bbsdat_patches[selind]
    selsourceind = []
    for i in range(len(bbsdat)):
        if bbsdat.Patch[i] in bbsdat_patches.Patch:
            selsourceind.append(i)
    bbsdat = bbsdat[selsourceind]
    return bbsdat, bbsdat_patches

########################################################################

#dat1, dat_patches1 = read_bbs('4C28.37.gsmsky6deg.model', verbose=verbose)
#dat1, dat_patches1 = remove_sources(dat1, dat_patches1, ra_centre, dec_centre, rad = 3.)
#dat2, dat_patches2 = read_bbs('4C28.37.gsmsky3deg.model', verbose=verbose)
#dat, dat_patches = append_sky_models(dat1, dat_patches1, dat2, dat_patches2)
#dat = sort_by_flux(dat)
#dat, dat_patches = put_sources_in_patches(dat, dat_patches)  
#dat, dat_patches = rename_patch(dat, dat_patches, ra_centre, dec_centre, rad = 10./3600., patchroot='4C28p37') 
#dat, dat_patches = rename_patch(dat, dat_patches, '15:04:59.40000000', '+26.00.46.90800000', rad = 10./3600., patchroot='3C310')
#plot_sky(dat, dat_patches, ra_centre, dec_centre)
#write_bbsdata_to_file(dat, dat_patches, '4C28.37.gsmsky6deg.comb.model')


