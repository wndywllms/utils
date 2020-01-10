# written by : Wendy Williams
# created : June 2013
# Utilities Library
# fits_util.py : routines for reading and handling FITS tables into python
#                record arrays

#import pyfits as pf
import astropy.io.fits as pf
import numpy as np

def  writeto(*args, **kwargs):
    '''
shorthand to astropy.io.fits.writeto    
    '''
    pf.writeto(*args, **kwargs)
    return
  
def readfits(*args, **kwargs):
    '''
    convenience - use astropy
    '''
    dat = pf.getdata(*args, **kwargs)
    return dat

def load_fits(fitsname, fields=[], ext=0, verbose=0):
    '''load_fits : load table data from fits file
Parameters
----------
required
  fitsname : str
     name of fitsfile to be read 
optional
  fields : list of str
    list of fields to be returned in array 
    if field ids are specified in fields list then load 
    only those fields (saves a bit on memory for large datasets)
  verbose : int
    print extra information
Returns
-------
  a : python record array
    data columns can be accessed by name
    '''
    if verbose > 0:
        print('loading data: %s' %(fitsname))
        
    # load fits and get data cast to recarray
    t = pf.open(fitsname, ext=ext)
    data = t[1].data
    t.close()
    #t = data.view(np.recarray)
    t = data
    
    if verbose > 2:
        print(data.dtype.names)
        
    # select only some fields
    if len(fields) > 0:
        # copy first column
        fi = t.dtype.names.index(fields[0])
        rec = np.recarray((len(t),),formats=t.dtype[fi],names=fields[0])
        rec[fields[0]] = t[fields[0]]
        
        a = rec
        
        # copy the other columns
        for fi, name in enumerate(fields[1:]):
            rec2 = t[name]
        
            d0=rec.dtype.descr
            d1=rec2.dtype.descr

            if rec.shape[0]==rec2.shape[0]:
                if len(d1)>1:
                    d2=d0+[(name,d1,rec2.shape[1:])]
                else:
                    d2=d0+[(name,d1[0][1])]
                a=np.zeros(rec.shape[0],dtype=d2)

            else:
                d2=d0+[(name,d1,(rec.shape[0],rec2.shape))]
                a=np.zeros(rec.shape[0],dtype=d2)

            for field in rec.dtype.fields:
                
                a[field] = rec[field]
                

            a[name]=rec2
            
            # update rec, to add next column            
            rec = a
        # print out the columns
        if verbose > 2:
            print(a.dtype.names)
    # choosing all columns
    else:
        a = t
            
    #if verbose > 2:
        #f = 1024.*1024 # for MB
        #print '              memory status: mem %.2f, resident %.2f, stacksize %.2f' %(memory()/f, resident()/f, stacksize()/f)
    
    #return a
    return a.view(np.recarray)
    
def append_field(rec, name, rec2):
  '''Append a field (column) to a record array

Parameters
----------
rec : recarray
   Input record arrya array
name  : str
   Name of field to be appended
rec2 : array
  Array to be appended as field "name"

Returns
-------
a : recarray
  Record Array with new appended field
  '''
  d0=rec.dtype.descr
  d1=rec2.dtype.descr

  if rec.shape[0]==rec2.shape[0]:
      if len(d1)>1:
          d2=d0+[(name,d1,rec2.shape[1:])]
      else:
          d2=d0+[(name,d1[0][1])]
      a=np.zeros(rec.shape[0],dtype=d2)

  else:
      d2=d0+[(name,d1,(rec.shape[0],rec2.shape))]
      a=np.zeros(rec.shape[0],dtype=d2)

  for field in rec.dtype.fields:
      
      a[field] = rec[field]

  a[name]=rec2

  #return a
  return a.view(np.recarray)
  
def append_record(rec, name, rec2):
  '''Append a record (row) to a record array. Assuming they have the same dtypes.

Parameters
----------
rec : recarray
   Input record arrya array
rec2 : recarray
  Array to be appended

Returns
-------
a : recarray
  Record Array with new appended field
  '''
  d0=rec.dtype.descr
  d1=rec2.dtype.descr

  if rec.shape[0]==rec2.shape[0]:
      if len(d1)>1:
          d2=d0+[(name,d1,rec2.shape[1:])]
      else:
          d2=d0+[(name,d1[0][1])]
      a=np.zeros(rec.shape[0],dtype=d2)

  else:
      d2=d0+[(name,d1,(rec.shape[0],rec2.shape))]
      a=np.zeros(rec.shape[0],dtype=d2)

  for field in rec.dtype.fields:
      
      a[field] = rec[field]

  a[name]=rec2

  return a.view(np.recarray)
