# encoding: utf-8
# mplFOAM module
# some useful functions to import, plot OpenFoam data with matplotlib

from __future__ import print_function

import os
import numpy as np
from vtk.util import numpy_support as npvtk
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()


def read_openfoam_case(directory=None,
                       filename=None,
                       mesh_regions=None,
                       cell_arrays=None,
                       times=None,
                       verbose=False):
    """
    Read the OpenFOAM case

    Parameters
    ----------
    directory : str, optional
                Directory path of the OpenFOAM case

    filename : str, optional
                 Filename of the .OpenfOAM case for Paraview

    mesh_regions : list of str
                List of the different regions of the mesh to be loaded

    cell_arrays : list of str
                List of the different fields to be loaded

    times : list of str
            List of the different timestep to be loaded

    Return
    ------
    openfoam_case_merged : paraview object
    """

    if directory is None:
        directory = os.getcwd().split('/')[-1]

    if filename is None:
        filename = directory + '.foam'

    case_file = open(filename, 'w')
    case_file.close()

    # Define the OpenFOAM data source
    openfoam_case = OpenFOAMReader(FileName=filename)

    # Decomposepolyhedra into tetrahedra and pyramids
    # Doesn't work to plot the contour with matplotlib 
    # openfoam_case.Decomposepolyhedra = 0

    # Import all the mesh regions if none are specified
    if mesh_regions is None:
        mesh_regions = openfoam_case.MeshRegions.Available
    openfoam_case.MeshRegions = mesh_regions

    # Import all the cell arrays if none are specified
    if cell_arrays is None:
        cell_arrays = openfoam_case.CellArrays.Available
    openfoam_case.CellArrays = cell_arrays

    # Read the specified time steps
    if times is None:
        # Get the latest timestep if none is specified
        openfoam_case.TimestepValues = openfoam_case.TimestepValues[-1]
    else:
        openfoam_case.TimestepValues = times

    # Print some basic informations on the loaded case
    if verbose:
        print("Case filename: ", filename)
        print("Mesh parts: ", mesh_regions)
        print("Volume fields: ", cell_arrays)
        print("Timestep values: ", openfoam_case.TimestepValues)

    openfoam_case_merged = MergeBlocks(openfoam_case)
    # openfoam_data = servermanager.Fetch(openfoam_case_merged)

    return openfoam_case_merged


def extractDataInPatch(case, verbose=False):
    "Extract the data from a patch and"

    # Select the active source
    # SetActiveSource(case)

    # Define the data
    Tetrahedralize1 = Tetrahedralize()
    data = servermanager.Fetch(Tetrahedralize1)

    # Define the number of points, cells and arrays
    nb_points = case.GetNumberOfPoints()
    nb_cells = case.GetNumberOfCells()
    nb_arrays = case.GetPointData().GetNumberOfArrays()

    # Get Paraview's own triangulation
    # cells = data.GetPolys()
    # triangles = cells.GetData()
    # nb_triangles = triangles.GetNumberOfTuples()/4

    nb_triangles = nb_cells
    tri = np.zeros((nb_triangles, 3))

    # for i in xrange(0, nb_triangles):
    #      tri[i, 0] = triangles.GetTuple(4*i + 1)[0]
    #      tri[i, 1] = triangles.GetTuple(4*i + 2)[0]
    #      tri[i, 2] = triangles.GetTuple(4*i + 3)[0]

    # Display some informations on the slice
    if verbose:

        # print 'Patch name :', patchName
        print("Number of cells:"), nb_cells
        print("Number of triangles:"), nb_triangles
        print("Number of points:"), nb_points
        print("Number of arrays:"), nb_arrays

        for i in range(nb_arrays):
            print("Array [', i, '] name:"), case.GetPointData().GetArrayName(i)

    # Put the poi(ts coordinates in x, y and z arrays
    x = []
    y = []
    z = []
    # points=zeros((nb_points,3))

    for i in range(nb_points):
        coord = case.GetPoint(i)
        xx, yy, zz = coord[:3]
        x.append(xx)
        y.append(yy)
        z.append(zz)

        # points[i,0] = xx
        # points[i,1] = yy
        # points[i,2] = zz

    x = np.array(x)
    y = np.array(y)
    z = np.array(z)

    # print("Points: "),points

    # Define the velocity components U=(u,v,w)
    U = npvtk.vtk_to_numpy(case.GetPointData().GetArray('U'))
    u = U[:, 0]
    v = U[:, 1]
    w = U[:, 2]

    # Define the cinematique pressure p => p/\rho
    p = npvtk.vtk_to_numpy(case.GetPointData().GetArray('p'))

    # Return the values extracted from the slice
    return x, y, z, u, v, w, p, tri


def extract_plane(case,
                  slice_origin=[0, 0, 0],
                  slice_normal=[0, 0, 1],
                  time=None,
                  verbose=False):
    """
    Extract Paraview slice object from OpenFOAM case `case`

    Parameters
    ----------
    case : Paraview object
            OpenFOAM case

    slice_origin : list of float
                    Slice point origin [x, y, z] in cartesian coordinates

    slice_normal : list of float
                    Slice point normal [nx, ny, nz] in cartesian coordinates

    Return
    ------
    plane_data: paraview slice object
    view = GetActiveView()
    >>> view.ViewTime = 3000
    AnimateReader
    () starts at the beginning and runs to
    the end with a fixed increment. You can change
    that and do your own start, end, and time increment.
    tsteps= reader.TimestepValues
    start = 2
    incr= 3
    end = 7
    for i in tsteps[start:end:incr]:
    view.ViewTime = tsteps[i]
    view.StillRender()
    ‏
    imgfile = “image.%03d.png” % (start+i*incr)
    view.WriteImage(imgfile“vtkPNGWriter”, 1)
    """
    # Select the active source
    SetActiveSource(case)

    # Define the plane
    Slice1 = Slice(Input=case)
    Slice1.SliceType = 'Plane'
    Slice1.SliceType.Origin = slice_origin
    Slice1.SliceType.Normal = slice_normal
    # Slice1.Triangulatetheslice = 0
    Slice1.UpdatePipeline(time=time)
    # Select the current time
    # animation_scene = GetAnimationScene()
    # view = GetActiveView()
    # if time is None:
    #     time = animation_scene.EndTime()
    print("Time=",time)
    # view.ViewTime = time
    # animation_scene.AnimationTime = time

    plane_data = servermanager.Fetch(Slice1)

    # Define the number of points, cells and arrays
    nb_points = plane_data.GetNumberOfPoints()
    nb_cells = plane_data.GetNumberOfCells()
    nb_arrays = plane_data.GetPointData().GetNumberOfArrays()

    # Display some informations on the slice
    if verbose:
        print("Slice origin: ", slice_origin)
        print("Slice normal: ", slice_normal)
        print("Number of cells: ", nb_cells)
        print("Number of points: ", nb_points)
        print("Number of arrays: ", nb_arrays)

        for i in range(nb_arrays):
            print("Array [", i, "] name:", plane_data.GetPointData().GetArrayName(i))

    return plane_data


def extract_plane_triangulation(plane_data, verbose=False):
    """Extract the triangulation of the plane

    Parameters
    ----------
    plane_data: Paraview slice object

    Return
    ------
    tri:
    triangles:
    """

    # Get Paraview's own triangulation
    cells = plane_data.GetPolys()
    triangles = cells.GetData()
    nb_triangles = triangles.GetNumberOfTuples() / 4

    tri = np.zeros((nb_triangles, 3))

    for i in xrange(0, nb_triangles):
        tri[i, 0] = triangles.GetTuple(4 * i + 1)[0]
        tri[i, 1] = triangles.GetTuple(4 * i + 2)[0]
        tri[i, 2] = triangles.GetTuple(4 * i + 3)[0]

    # Display some informations on the triangulation
    if verbose:
        print("Number of triangles: "), nb_triangles

    # Return the triangulation of the plane
    return tri, triangles


def extract_plane_points(plane_data, verbose=False):
    """Extract the points cartesian coordinates (x,y,z) of the plane

    Parameters
    ----------
    plane_data: Paraview slice object

    Return
    ------
    x,y,z: arrays of floats
            cartesian coordinates of the points in the plane
    """

    # Get the number of points
    nb_points = plane_data.GetNumberOfPoints()

    # Put the points coordinates in x, y and z arrays
    x = []
    y = []
    z = []

    for i in range(nb_points):
        coord = plane_data.GetPoint(i)
        xx, yy, zz = coord[:3]
        x.append(xx)
        y.append(yy)
        z.append(zz)

    x = np.array(x)
    y = np.array(y)
    z = np.array(z)

    # Display some informations on the triangulation points
    if verbose:
        print("Number of points:", nb_points)
        print("xrange: [%5.3f, %5.3f]") % (x.min(), x.max())
        print("yrange: [%5.3f, %5.3f]") % (y.min(), y.max())
        print("zrange: [%5.3f, %5.3f]") % (z.min(), z.max())

    # Return the slice points coordinates
    return x, y, z


def extract_plane_vector_field(plane_data, field_name, verbose=False):
    """Extract the vector field components in the plane `plane_data`

    Parameters
    ----------
    plane_data: Paraview slice object

    field_name: str
                Name of the vector field

    Return
    ------
    u,v,w: arrays of floats
           Vector field components
    """

    # Define the velocity components U=(u,v,w)
    U = npvtk.vtk_to_numpy(plane_data.GetPointData().GetArray(field_name))
    u = U[:, 0]
    v = U[:, 1]
    w = U[:, 2]

    # Return the vector field components
    return u, v, w


def extract_plane_scalar_field(plane_data, field_name, verbose=False):
    """Extract the scalar field `field_name` in the plane `plane_data`

    Parameters
    ----------
    plane_data: Paraview slice object

    field_name: str
                Name of the scalar field

    Return
    ------
    scalar: array of floats
            Scalar field
    """

    # Put the scalar field field_name in array scalar
    scalar = npvtk.vtk_to_numpy(plane_data.GetPointData().GetArray(field_name))

    # Return the scalar field values extracted from the plane
    return scalar

def cart2pol(x, y):
    """Compute the polar coordinates `field_name` in the plane `plane_data`

    Parameters
    ----------
    x and y cartesian coordinates

    Return
    ------
    r and theta polar coordinates
    """
    r = np.sqrt(x**2+y**2)
    #phi = arctan(y/x)+pi*heav(-x)*sgn(y)
    theta = np.arctan2(y,x)
    return r, theta
    
def cart2pol_vector(u,v,theta):
    uradial  = np.cos(theta)*u + np.sin(theta)*v
    utheta = -np.sin(theta)*u + np.cos(theta)*v
    return utheta, uradial

def pol2cart(theta, rho):
    x = rho * np.cos(theta)
    y = rho * np.sin(theta)
    return x, y

def cart2sph(x, y, z):
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    el = np.arctan2(z, hxy)
    az = np.arctan2(y, x)
    return az, el, r

def sph2cart(az, el, r):
    rcos_theta = r * np.cos(el)
    x = rcos_theta * np.cos(az)
    y = rcos_theta * np.sin(az)
    z = r * np.sin(el)
    return x, y, z