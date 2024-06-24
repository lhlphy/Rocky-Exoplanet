import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from parameter_list import *
from function_library import *
import multiprocessing
import time

# global Intensity
# Intensity = np.zeros((len(phiP_list), len(thetaP_list)))
# Use 'TkAgg' as the backend for Matplotlib
#matplotlib.use('TkAgg')

# Main function


# angle = np.linspace(0, np.pi / 2, 100)
# I = np.zeros(len(angle))
# for i in range(len(angle)):
#     I[i] = I_dist_afterreflect(R1,a,angle[i],camera,camera)
# plt.plot(angle, I)


def main(i, j, Intensity, diffuse_ratio):
    
    phiP = phiP_list[i]
    thetaP = thetaP_list[j]
    # Calculate the normal vector and position
    nv, Pos, r = normal_vec(phiP, thetaP, Theta, a, e, R2)
    # Calculate the reflected vector
    RV = reflect(Pos, nv)
    
    # Check if the reflection direction is towards the camera
    if check_direction(RV, nv, camera, Pos):
        # Calculate the angle between the camera and the reflected vector
        #angle = angle_between(camera, RV)
        # Calculate the intensity of the reflected light
        Diffuse = Oren_Nayar_BRDF(R1, r, nv, Pos,camera, 6000,1e-6)
        SR  = specular_reflection(0.5, RV, camera, nv, r, 6000, 1e-6)
        with Intensity.get_lock():
            Intensity[SIZE[1]*i+j] = (Diffuse + SR) #* blackbody_radiation(6000, 1e-6)

        with diffuse_ratio.get_lock():
            if Diffuse == 0:
                diffuse_ratio[SIZE[1]*i+j] = 0
            else:
                diffuse_ratio[SIZE[1]*i+j] = Diffuse / (Diffuse + SR)
    # else:
    #     Intensity[i, j] = 0
    #print(Intensity[SIZE[1]*i+j])


if __name__ == '__main__':
    t1 = time.time()
    processes = []
    Intensity = multiprocessing.Array('d', SIZE[0]*SIZE[1])   
    diffuse_ratio = multiprocessing.Array('d', SIZE[0]*SIZE[1])

    # Loop through all points on the planet's surface
    for i, phiP in enumerate(phiP_list):
        for j, thetaP in enumerate(thetaP_list):
            process = multiprocessing.Process(target=main, args=(i, j, Intensity, diffuse_ratio))
            processes.append(process)
            process.start()

        for process in processes:
            process.join()
        
        processes = []
        print(i)


    Intensity = Intensity* blackbody_radiation(6000, 1e-6)
    Intensity = np.array(Intensity[:]).reshape(SIZE[0], SIZE[1])

    print(Intensity)
    #Intensity = Intensity / np.max(Intensity)
    # Create a sphere plot
    phiP, thetaP = np.meshgrid(phiP_list, thetaP_list)
    x = R2 * np.cos(phiP) * np.cos(thetaP)
    y = R2 * np.cos(phiP) * np.sin(thetaP)
    z = R2 * np.sin(phiP)

    # Plotting the sphere
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot the surface with intensity as color
    mappable = ax.plot_surface(x, y, z, facecolors=plt.cm.gray(Intensity.T /np.max(Intensity)), rstride=1, cstride=1, antialiased=False)

    # Plot the incident and reflected vectors
    ax.quiver(-(2000 + R2) * np.cos(Theta), -(2000 + R2) * np.sin(Theta), 0, np.cos(Theta), np.sin(Theta), 0, color='r', length=2000.0, normalize=True)
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
    
    t2 = time.time()
    print("Program run time:",t2-t1,'s')
    

#A line passes point Pos, with the direction vector Camera. Calculate the distance between this line and the origin.
#The line is defined by the equation: r = Pos + t*Camera
# #The distance between the line and the origin is given by: |Pos x Camera|/|Camera|
# print(np.linalg.norm(np.cross(Pos, camera))/np.linalg.norm(camera))
# #The distance between the line and the origin is given by: |Pos| sin(theta)
# print(np.linalg.norm(Pos)*np.sin(angle_between(Pos, camera)))

