# encoding: utf-8
# Extract data from a plane
import matplotlib.pyplot as plt
import math
import numpy as np
from mplfoam import *

# Open the foam case
case = read_openfoam_case()

# Define the plane at position z_pos from which the data are extracted
z_pos = 6.0
origin = [0, 0, z_pos]
normal = [0, 0, 1]
time = 3000

# Extract the defined plane
plane = extract_plane(case, slice_origin=origin, slice_normal=normal, time=time)

# Define the velocity components U=(u,v,w)
u, v, w = extract_plane_vector_field(plane, 'U')

# Extract the cartesian coordinates of the triangulation points
x,y,z = extract_plane_points(plane, verbose=False)

# Extract the triangulation of the plane
tri,triangles = extract_plane_triangulation(plane)

# Extract the pressure
p = extract_plane_scalar_field(plane, 'p')

# Compute the velocity magnitude
u_magnitude = np.sqrt(u**2 + v**2 + w**2)

# Define geometric and kinematic parameters
a = 0.095
b = 0.0475
b_a = b/a
h = (1-b_a)/(1+b_a)
S = math.pi*a*b
E = math.pi*(a + b)*(3 - math.sqrt(4-h**2))
Dh = 4*S/E
print('Diamètre hydraulique Dh = ',Dh)

U0 = 0.1323

# Draw the mesh
plt.figure("Mesh")
plt.gca().set_aspect('equal')
plt.triplot(x, y,tri, lw=0.5, color='black')
plt.scatter(x, y, s=20, c='red', alpha=0.9)
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.title(r"Mesh at $z=$" + str(z_pos))
plt.grid(linestyle=':')

# Draw the velocity contour with mesh
plt.figure("Velocity")
plt.gca().set_aspect('equal')
plt.tripcolor(x, y, tri, w, shading='flat', edgecolors='')
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.title(r"Velocity contours at $z=$" + str(z_pos))
plt.grid(linestyle=':')
cbar = plt.colorbar()
cbar.set_label('Velocity', labelpad=10)

# Draw the kinematic pressure contour with mesh
plt.figure("Pressure")
plt.gca().set_aspect('equal')
plt.tripcolor(x, y, tri, p, shading='flat', edgecolors='')
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.title(r"Kinematic pressure at $z=$" + str(z_pos))
plt.grid(linestyle=':')
cbar = plt.colorbar()
cbar.set_label('Kinematic Pressure', labelpad=10)

# Draw the velocity vector field
plt.figure("Velocity Vectors")
plt.gca().set_aspect('equal')
plt.tripcolor(x, y, tri, np.zeros(np.size(x)), cmap=plt.cm.gray_r)
plt.quiver(x, y, u, v)
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.title(r"Velocity vectors at $z=$" + str(z_pos))
plt.grid(linestyle=':')

#
# PlotOverLine curve plots
#

# Define points coordinates in cylindrical coordinates
r1 = 0.0
r2 = b

theta1 = np.deg2rad(90)
theta2 = np.deg2rad(90)

z1 = 0
z2 = 0

# Transform points coordinates from cylindrical to cartesian
x1 = r1*math.cos(theta1)
x2 = r2*math.cos(theta2)

y1 = r1*math.sin(theta1)
y2 = r2*math.sin(theta2)

plt.figure("Velocity plot")
plt.figure("Temperature plot")

for i in range(0,11,2):

    # Define the source data
    SetActiveSource(case)

    # Define the line
    PlotOverLine1 = PlotOverLine(Source="High Resolution Line Source")
    zPos = i*Dh
    PlotOverLine1.Source.Point1 = [x1, y1, z1 + zPos]
    PlotOverLine1.Source.Point2 = [x2, y2, z2 + zPos]

    # Transform the paraview arrays into numpy ones
    dataLine = servermanager.Fetch(PlotOverLine1)
    ULine = npvtk.vtk_to_numpy(dataLine.GetPointData().GetArray('U'))
    T = npvtk.vtk_to_numpy(dataLine.GetPointData().GetArray('T'))
    arcLength = npvtk.vtk_to_numpy(dataLine.GetPointData().GetArray('arc_length'))

    # Plots
    plt.figure("Velocity plot")
    plt.figure
    plt.plot(ULine[:,2]/max(ULine[:,2]),arcLength/max(arcLength),label='z/d='+str(i))
    # xlim(0,1)
    # ylim(0,1.1)
    plt.xlabel(r'$Uz/Uz_{max}}$')
    plt.ylabel(r'$r/R$')
    plt.title("Radial velocity distribution")
    plt.legend()
    plt.grid(linestyle=':')

    plt.figure("Temperature plot")
    plt.plot(T, arcLength/max(arcLength), label='z/d='+str(i))
    # xlim(0,1)
    # ylim(0,1.1)
    plt.xlabel(r'$T$')
    plt.ylabel(r'$r/R$')
    plt.title("Radial temperature distribution")
    plt.legend()
    plt.grid(linestyle=':')

#
# Plot On Intersection Curves
#

plt.figure()

for i in range(1,11,2):

    # Define the source data
    SetActiveSource(case)

    # Define the line
    plotOnIntersectionCurves1 = PlotOnIntersectionCurves(Input=case)
    plotOnIntersectionCurves1.SliceType = 'Plane'
    z_pos = i*Dh
    plotOnIntersectionCurves1.SliceType.Origin = [0.0, 0.0, z_pos]
    plotOnIntersectionCurves1.SliceType.Normal = [0.0, 0.0, 1.0]

    # Transform the paraview arrays into numpy ones
    dataLine = servermanager.Fetch(plotOnIntersectionCurves1)
    # Need this to access the point data
    data = dataLine.GetBlock(0)
    T = npvtk.vtk_to_numpy(data.GetPointData().GetArray('T'))
    tau_w = npvtk.vtk_to_numpy(data.GetPointData().GetArray('wallShearStress'))
    # Compute the velocity magnitude
    tau_w_magnitude = np.sqrt(tau_w[:,0]**2 + tau_w[:,1]**2 + tau_w[:,2]**2)
    arcLength = npvtk.vtk_to_numpy(data.GetPointData().GetArray('arc_length'))
    # Get the points coordinates
    points = data.GetPoints()
    pos = npvtk.vtk_to_numpy(points.GetData())
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]
    # Compute the polar coordinates
    r, theta = cart2pol(x,y)
    # Plot
    # plt.plot(theta, T, label='z/d='+str(i))
    plt.polar(theta, tau_w_magnitude, label='z/d='+str(i))
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'${\tau}_{w}$')
    plt.title("Azimuthal wall shear stress distribution")
    plt.legend()
    plt.grid(linestyle=':')

# Show all the figures
plt.show()
