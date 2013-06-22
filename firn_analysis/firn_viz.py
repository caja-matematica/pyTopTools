from pyTopTools import pyImage as P
#from pyTopTools import convert2vtk as C
import numpy as np


def bmp_block( bottom, top, height, fpath, fprefix ):
    """
    Returns list of frames in a block
    """
    ims = []
    # grab each image from base --> base+height
    for x in range( height ):
        image_num = bottom + x
        im = P.PyImage( fpath + fprefix + str( image_num ) + '.bmp' )
        im.bmp2array()
        ims.append( im.data )
    return ims

def arr2vtk( data, vtkname ):
    """
    Convert numpy 3D array to VTK format

    cim : ChompImage object for cubical corners
    """
    # convert ims --> array!
    #cim.get_cubical_corners()
    
    pts = C.data2vtk( data, vtkname, dtype='arr' )
    return pts


if __name__ == "__main__":

    bmp_path = '/data/CT_Firn_Sample/output23-10-3/'
    prefix = 'K09b-23-10-'

    h = 30
    # input from shell script
    low = 3200
    high = low + h
    
    print "height = ", h
    print low, high
    
    # do the homology computations
    B = bmp_block( low, high, h, bmp_path, prefix )
    #cim = P.ChompImage( np.array( B ) )
    
    #vpts = arr2vtk( np.asarray( B ), 'vtktest' )
    
    # Requires Canopy or at least tvtk installed
    from pyTopTools import structured_3d_VTK as svtk
    svtk.numpy2vtk( np.array( B ), 'firn3d' )
