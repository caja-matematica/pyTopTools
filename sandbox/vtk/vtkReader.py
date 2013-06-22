import vtk
import numpy as np

def vtk2numpy( coord ):
    sc = coord.strip().split(' ')
    pt = [ float(x) for x in sc]
    return np.array( pt )

def find_volume( voxel ):
    """
    voxel holds numpy array with the eight points defining a voxel in
    VTK format.

    Return volume of voxel.
    """
    p = voxel 
    length = np.absolute( p[1][0]-p[0][0] )
    width = np.absolute( p[2][1]-p[0][1] )
    height = np.absolute( p[4][2]-p[0][2] )
    return  length * width * height

def construct_voxel( base_point, mx, my ,mz ):
    """
    Construct a VTK voxel given a base point in the lower left corner.
    """
    p = base_point
    voxel = ( p,
             [p[0]+mx, p[1], p[2]],
             [p[0], p[1]+my, p[2]],
             [p[0]+mx, p[1]+my, p[2]],
             [p[0], p[1], p[2]+mz],
             [p[0]+mx, p[1], p[2]+mz],
             [p[0], p[1]+my, p[2]+mz],
             [p[0]+mx, p[1]+my, p[2]+mz],
             )
    return np.asarray( voxel )

def subdivide( voxel ):
    """
    subdivide a voxel
    """
    new_pts = [voxel[0]]
    # this creates a new sub-voxel with points arranged in the same
    # order as the original
    for v in voxel[1:]:
        new_pts.append( voxel[0] + 0.5 * (v-voxel[0]) )
    
    # mx = np.mean( [p[0][0], p[1][0]] )
    # my = np.mean( [p[0][1], p[2][1]] )
    # mz = np.mean( [p[0][2], p[4][2]] )
    # new_pts = construct_voxel( p[0], mx, my, mz )
    return np.asarray( new_pts )
    # voxel_list.append( new_pts )
    
    # print "new pts", new_pts
    # print "volume", volume
    # print min_volume

    # while volume > min_volume:
    #     subdivide( new_pts, min_volume, voxel_list=voxel_list )
    # return voxel_list


def run():
    fname = 'octree_cylinder_in_cylinder.vtk'

    reader = vtk.vtkDataSetReader()
    reader.SetFileName( fname )
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    #reader.ReadCellData()

    reader.Update()
    
    data=reader.GetOutput()
    npoints=data.GetNumberOfPoints()
    ncells = data.GetNumberOfCells()
    nscalars=reader.GetNumberOfScalarsInFile()
    nvectors=reader.GetNumberOfVectorsInFile()

    d = data.GetCellData()
    narrays = d.GetNumberOfArrays()
    #d = data.GetPointData()
    cell_data = d.GetArray( d.GetArrayName( 0 ) )
    comp = cell_data.GetNumberOfComponents()

    # grab POINTS
    pdata = data.GetPoints().GetData()
    psize = pdata.GetSize()
    points = np.asarray( [ pdata.GetValue( i ) for i in xrange( psize ) ] )
    points.resize( (psize/3, 3) )

    # cell values (scalars)
    d = data.GetCellData()
    narrays = d.GetNumberOfArrays()
    csize = cell_data.GetNumberOfTuples()
    cell_data = d.GetArray( d.GetArrayName( 0 ) )
    cellvals = np.asarray( [ cell_data.GetValue( i ) for i in xrange( csize ) ] )

    # grab the cells themselves -- defined by indices of points
    c = data.GetCells()
    cellsize = c.GetSize()
    cell_arr = c.GetData()
    # first element of every row lists the number of vertices of the voxel
    vsize = cell_arr.GetValue( 0 ) + 1
    # ignore first entry in cell deinition
    cells = np.asarray( [ cell_arr.GetValue( i ) for i in xrange( cellsize) ] )
    cells.resize( (cellsize/vsize, vsize) )

    # form a equispaced grid and map cellvals to grid elements
    px,py = points.shape
    nx, ny, nz = [ np.unique( points[:,i] ) for i in range( py ) ]
    ixgrid = np.ix_( nx, ny, nz )

    # map cell corners to scalars values
    coords = []
    for idx, voxel in enumerate( cells ):
        coords.append( points[ voxel[1] ].tolist() + [cellvals[ idx ]] )
        #coords = np.asarray( coords )

    # grid size
    smidge = 0.0000001
    dx, dy, dz = [ np.diff( g ) for g in [nx,ny,nz] ]
    gsize = dx.min()
    cell_threshold = gsize + 0.001
    min_vol = gsize**3 + smidge

    # For each cell defined by the corners [v0,...,v7], see if it needs to
    # be refined.  If not, make the lower left corner (v0, see VTK
    # documentation, "voxel") a point in the cubicle complex, and give it
    # the scalar value of the cell in question. If the cell is too large,
    # first refine it, making sure that each subcell has the same scalar
    # value as its parent. Then include each of the four lower left
    # corners in the cubicle complex.
    refined_cells = []
    for idx, c in enumerate( cells ):
        p = points[ c[1:] ]
        # points defining length, width, height of voxel
        volume = find_volume( p )
        if volume > min_vol:
            new_pts = subdivide( p )
            for pt in new_pts:
                refined_cells.append( pt.tolist() + [cellvals[idx]] )
        else:
            refined_cells.append( p[0].tolist() + [cellvals[idx]] )
    refined_grid = np.asarray( refined_cells )
    return refined_grid, cellvals, dx[0], dy[0], dz[0]

if __name__ == "__main__":

    import chomp_betti as cb

    # test_voxel = np.array( [[0,0,0], [1,0,0], [0,1,0], [1,1,0],
    #                          [0,0,1], [1,0,1], [0,1,1], [1,1,1]], dtype=np.float )
    # test_voxel *= 1./3
    # vol = find_volume( test_voxel )
    # pts = subdivide( test_voxel )
    
    # print "voxel\n", test_voxel
    # print "voxel vol", vol
    # print ""
    # print "sub voxel\n", pts
    # print "vol", find_volume( pts )
    # print ""

    # dx = np.absolute( pts[1][0]-pts[0][0] )
    # print "dx", dx

    # scaled  = pts * 1./dx
    # print "scaled\n", scaled

    pts, scalars, dx, dy, dz = run()
    pts[:,:-1] *= 1./dx
    # pts[:,0] *= 1./dx
    # pts[:,1] *= 1./dy
    # pts[:,2] *= 1./dz

    uniq_scalars = np.unique( scalars )
    colors_idx = [ np.where( pts[:,-1] == s )[0]
                  for s in uniq_scalars ]
    cubset = [ np.array( pts[ colors_idx[i], :-1 ], dtype=np.uint )
                 for i,s in enumerate( uniq_scalars ) ]
    fname = 'cylinder'
    for i,s in enumerate( uniq_scalars ):
        cb.array2chomp( cubset[i], fname+str(s)+'.cub' )
