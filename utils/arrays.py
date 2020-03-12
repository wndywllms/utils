import numpy as np

def find_x_in_y(x, y, ydata=None):
    '''
    find the corresponding elements (rows or subarrays) in ydata indexed by y that match to x
    also handles if values in x are missing from y
    
    length of y should correspond to the first dimension of xdata
    returns ydata_x
    
    Example:
    
        x = np.array([2,1,5,10,100,6])
        y = np.array([3,5,7,100,9,8,6,6])
        ydata = np.array(['a','b','c','d','e','f','g','h'])
    
    find_x_in_y(x, y) returns
        masked_array(data=[--, --, 5, --, 100, 6],
              mask=[ True,  True, False,  True, False, False], 
              fill_value=999999)
        array([ True,  True, False,  True, False, False])

    
    find_x_in_y(x, y, ydata) returns
        masked_array(data=[--, --, 'b', --, 'd', 'g'],
              mask=[ True,  True, False,  True, False, False],
              fill_value='N/A', dtype='|S1')
        array([ True,  True, False,  True, False, False])
     
    '''
    
    if ydata is None:
        ydata = y.copy()
    
    assert len(y) == ydata.shape[0]
    

    index_ysort = np.argsort(y)
    sorted_y = y[index_ysort]
    sorted_ydata = ydata[index_ysort]
    
    sorted_index = np.searchsorted(sorted_y, x)

    yindex_of_x = np.take(index_ysort, sorted_index, mode="clip")
    #ydata_of_x = np.take(sorted_ydata, sorted_index, mode="clip", axis=0)
    
    # where the names do not match
    mask_no_value = y[yindex_of_x] != x 

    y_x_index = np.ma.array(yindex_of_x, mask=mask_no_value)
    ydata_x = ydata[y_x_index]
    
    return ydata_x, mask_no_value
