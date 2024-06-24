import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from parameter_list import *
from function_library import *

# Use 'TkAgg' as the backend for Matplotlib
matplotlib.use('TkAgg')


# Main function
def planet_intensity(Orbital_angle):

    # angle = np.linspace(0, np.pi / 2, 100)
    # I = np.zeros(len(angle))
    # for i in range(len(angle)):
    #     I[i] = I_dist_afterreflect(R1,a,angle[i],camera,camera)
    # plt.plot(angle, I)

    # Create meshgrid for the planet
    phiP_list = np.linspace(-np.pi / 2, np.pi / 2, 180)
    thetaP_list = np.linspace(0, 2 * np.pi, 360)
    Intensity = np.zeros((len(phiP_list), len(thetaP_list)))

    # Loop through all points on the planet's surface
    for i, phiP in enumerate(phiP_list):
        for j, thetaP in enumerate(thetaP_list):
            # Calculate the normal vector and position
            nv, Pos, r = normal_vec(phiP, thetaP, Orbital_angle, a, e, R2)
            # Calculate the reflected vector
            RV = reflect(Pos, nv)
            
            # Check if the reflection direction is towards the camera
            if check_direction(RV, nv, camera, Pos):
                # Calculate the angle between the camera and the reflected vector
                angle = angle_between(camera, RV)
                # Calculate the intensity of the reflected light
                Intensity[i, j] = I_dist_afterreflect(R1, r, angle, nv, RV, Pos)
            else:
                Intensity[i, j] = 0

    # Create a plot
    planet_drawer()


    def planet_drawer():
        # Create a sphere plot
        phiP, thetaP = np.meshgrid(phiP_list, thetaP_list)
        x = R2 * np.cos(phiP) * np.cos(thetaP)
        y = R2 * np.cos(phiP) * np.sin(thetaP)
        z = R2 * np.sin(phiP)

        # Plotting the sphere
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Plot the surface with intensity as color
        mappable = ax.plot_surface(x, y, z, facecolors=plt.cm.gray(Intensity.T), rstride=1, cstride=1, antialiased=False)

        # Plot the incident and reflected vectors
        ax.quiver(-(2000 + R2) * np.cos(Orbital_angle), -(2000 + R2) * np.sin(Orbital_angle), 0, np.cos(Orbital_angle), np.sin(Orbital_angle), 0, color='r', length=2000.0, normalize=True)
        ax.quiver(R2 * camera[0], R2 * camera[1], R2 * camera[2], camera[0], camera[1], camera[2], color='g', length=2000.0, normalize=True, linestyle='dashed')

        # Set axis labels
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.set_zlabel('Z (km)')

        # Set the view angle 
        elev = np.arcsin(camera[2])
        if camera[2] == 1:
            azim = 0
        else:
            azim = np.arccos(camera[0]/np.cos(elev))
        ax.view_init(elev=np.rad2deg(elev) , azim=np.rad2deg(azim))

        # Show the plot
        plt.show()

    #A line passes point Pos, with the direction vector Camera. Calculate the distance between this line and the origin.
    #The line is defined by the equation: r = Pos + t*Camera
    # #The distance between the line and the origin is given by: |Pos x Camera|/|Camera|
    # print(np.linalg.norm(np.cross(Pos, camera))/np.linalg.norm(camera))
    # #The distance between the line and the origin is given by: |Pos| sin(theta)
    # print(np.linalg.norm(Pos)*np.sin(angle_between(Pos, camera)))

    


