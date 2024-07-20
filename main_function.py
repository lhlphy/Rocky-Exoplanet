import numpy as np
import matplotlib.pyplot as plt
from parameter_list import *
from function_library import *
import multiprocessing
import os
from matplotlib import rcParams
from matplotlib.pylab import mpl
from scipy.interpolate import interp1d
import time
#import time

# Use 'TkAgg' as the backend for Matplotlib
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


def BRDF(i, j, Intensity, I_diffuse, Theta, Coarse, Model='Lambert', semaphore = None):
    """
    calculate the reflection and diffusion intensity  (divided by B(T,lam)) of the planet surface at the point (i,j)
    i: the index of phiP
    j: the index of thetaP
    Intensity: the shared memory of the Intensity
    I_diffuse: the shared memory of the I_diffuse
    Theta: orbital phase angle
    Coarse: the coarse of the surface, the standard derivation of incline angle of the surface (0-pi/2)
    Temperature: the temperature of the star
    Model: the model of the reflection and diffusion
        switch:
        1. Lambert: the Lambertian model Coarse = 0
        2. Oren_Nayar: the Oren-Nayar model
        3. Gaussian_wave: the Gaussian wave model (depends on the surfave wind speed)

    """
    semaphore.acquire()
    phiP = phiP_list[i]
    thetaP = thetaP_list[j]
    # Calculate the normal vector and position
    nv, Pos, r = normal_vec(phiP, thetaP, Theta, a, e, R2)

    # if check_intersection_with_star(Pos, camera):  # Check if the line intersects with the star--Check block
    #     with Intensity.get_lock():
    #         Intensity[SIZE[1]*i+j] = 0

    #     return
    
    # Calculate the reflected vector
    RV = reflect(Pos, nv)
    
    # Check if the reflection direction is towards the camera
    if check_direction(RV, nv, camera, Pos):
        # Calculate the angle between the camera and the reflected vector
        # angle = angle_between(camera, RV)
        # Calculate the intensity of the reflected light
        # Model Choice 
        if Model == 'Lambert':   #Coarse = 0
            Diffuse = Oren_Nayar_BRDF(R1, r, nv, Pos, camera, 0 )
            SR  = specular_reflection(RV, camera, nv, r)
            # SR is the reflected light intensity divided by B(T,lam)
        elif Model == 'Oren_Nayar':
            Diffuse = Oren_Nayar_BRDF(R1, r, nv, Pos, camera, Coarse)
            SR  = specular_reflection(RV, camera, nv, r)
        elif Model == 'Gaussian_wave':  # In this model Diffuse and RF are considered together
            Diffuse = Wave_reflect(R1, r, nv, Pos, camera )
            SR = 0

    
        with Intensity.get_lock():
            Intensity[SIZE[1]*i+j] = (Diffuse + SR) #* blackbody_radiation(6000, 1e-6)
            # Intensity : total intensity divided by B(T,lam)   

        with I_diffuse.get_lock():
            I_diffuse[SIZE[1]*i+j] = Diffuse
    # else:
    #     Intensity[i, j] = 0
    #print(Intensity[SIZE[1]*i+j])
    semaphore.release()

@decorator_timer('global_intensity')
def global_intensity(Theta, Coarse = Coarse_g, id=0, Model = 'Lambert', mode = 'geo'):
    """
    Calculate the intensity map of the reflection and diffusion of the planet surface
    mode: the mode of the calculation
        switch: "geo" or "phy"
        "geo": It's a gemotry problem, don't need to consider the thermal radiation and albedo
            return:
            Intensity: the intensity of the reflection and diffusion/B(T,lam)
            I_diffuse: the intensity of the diffusion/B(T,lam)
            Intensity - I_diffuse: the intensity of the reflection/B(T,lam)

        "phy": It's a physical problem, need to consider the thermal radiation and albedo
            return:
            Intensity: the intensity of the reflection and diffusion/B(T,lam)
            I_diffuse: the intensity of the diffusion/B(T,lam)
            Intensity - I_diffuse: the intensity of the reflection/B(T,lam)

            Wavelength and Temperature are from the parameter_list.py
        
    """
    # if mode == 'geo':  # only consider the geometry problem
    #     SPE_REF = 1
    #     DIF_REF = 1
    max_processes = 900000
    semaphore = multiprocessing.Semaphore(max_processes)
    processes = []
    Intensity = multiprocessing.Array('d', SIZE[0]*SIZE[1])   # diffusion + reflection
    I_diffuse = multiprocessing.Array('d', SIZE[0]*SIZE[1]) # diffusion  intensity

    # Loop through all points on the planet's surface
    #calculate the intensity of the reflect and diffusion using the BRDF function
    for i, phiP in enumerate(phiP_list):
        for j, thetaP in enumerate(thetaP_list):
            process = multiprocessing.Process(target = BRDF, args=(i, j, Intensity, I_diffuse, Theta, Coarse, Model, semaphore))
            processes.append(process)
            process.start()

    for process in processes:
        process.join()

    if mode == 'phy':
        Intensity = Intensity* blackbody_radiation(Temperature, Wavelengh)
        I_diffuse = I_diffuse* blackbody_radiation(Temperature, Wavelengh)

    Intensity = np.array(Intensity[:]).reshape(SIZE[0], SIZE[1])
    I_diffuse = np.array(I_diffuse[:]).reshape(SIZE[0], SIZE[1])

    # print(Intensity)
    # Intensity = Intensity / np.max(Intensity)
    # Create a sphere plot
    phiP, thetaP = np.meshgrid(phiP_list, thetaP_list)
    x = R2 * np.cos(phiP) * np.cos(thetaP)
    y = R2 * np.cos(phiP) * np.sin(thetaP)
    z = R2 * np.sin(phiP)

    # Plotting the sphere
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot the surface with intensity as color
    temp = np.max(Intensity)
    if temp == 0:
        mappable = ax.plot_surface(x, y, z, facecolors=plt.cm.gray(Intensity.T ), rstride=1, cstride=1, antialiased=False)
    else:
        mappable = ax.plot_surface(x, y, z, facecolors=plt.cm.gray(Intensity.T /temp), rstride=1, cstride=1, antialiased=False)

    # Plot the incident and reflected vectors
    # ax.quiver(-(2000 + R2) * np.cos(Theta), -(2000 + R2) * np.sin(Theta), 0, np.cos(Theta), np.sin(Theta), 0, color='r', length=2000.0, normalize=True)
    # ax.quiver(R2 * camera[0], R2 * camera[1], R2 * camera[2], camera[0], camera[1], camera[2], color='g', length=2000.0, normalize=True, linestyle='dashed')

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
    #plt.show()
    #save the plot to temp/ folder
    os.makedirs(f'temp/R{id}/plots', exist_ok=True)
    name = f'temp/R{id}/plots/plot_'+str(int(Coarse*180/np.pi))+'_'+str(int(Theta*180/np.pi))+'.png'
    plt.savefig(name)
    plt.close()
    
    return Intensity, I_diffuse, Intensity - I_diffuse
    #print("Program run time:",t2-t1,'s')
    

#A line passes point Pos, with the direction vector Camera. Calculate the distance between this line and the origin.
#The line is defined by the equation: r = Pos + t*Camera
# #The distance between the line and the origin is given by: |Pos x Camera|/|Camera|
# print(np.linalg.norm(np.cross(Pos, camera))/np.linalg.norm(camera))
# #The distance between the line and the origin is given by: |Pos| sin(theta)
# print(np.linalg.norm(Pos)*np.sin(angle_between(Pos, camera)))
    
@decorator_timer('thermal_spectrum')
def thermal_spectrum(wavelength_bound, Temperature= Temperature , id=0, Ntheta = 5, NWavelength = 1):
    """
    Calculate the blackbody radiation spectrum

    Ratio: [size]: NWavelength * (2*Ntheta)
    Theta_list: [size]: (2*Ntheta)
    """ 
    t10 = time.time()
    # Planck's constant
    h = 6.62607015e-34
    # Speed of light
    c = 299792458
    # Boltzmann constant
    k = 1.380649e-23
    # Wavelengths
    if Ntheta == 1:
        Theta_list = np.array([np.pi])
    else:
        Theta_list = np.linspace(0, np.pi, Ntheta)  # 0-pi 与 pi-2pi 重复
    Wave_list = np.linspace(wavelength_bound[0], wavelength_bound[1], NWavelength)
    TMAP0 = Tmap(0, id)
    t11 = time.time()
    # processes = []
    # spectrum_P = multiprocessing.Array('d', len(Wavelength))   
    # spectrum_S = multiprocessing.Array('d', len(Wavelength))
    SP = np.zeros([len(Theta_list), len(Wave_list)])
    SS = np.zeros([len(Theta_list), len(Wave_list)])
    RAT = np.zeros([len(Theta_list), len(Wave_list)])

    for i, Theta in enumerate(Theta_list):
        if e == 0:
            TMAP = TMAP0
        else:
            TMAP = Tmap(Theta, id)

        spectrum_P = np.zeros(len(Wave_list)) 
        spectrum_S = np.zeros(len(Wave_list)) 
        for j, wavelength in enumerate(Wave_list):
            # Calculate the blackbody radiation spectrum
            if Theta < 1e-6:
                spectrum_P[j] = 0
            else:
                spectrum_P[j] = Radiation_cal(TMAP, Theta, camera, Temperature, wavelength)
            spectrum_S[j] = Cal_star_flux(Theta, wavelength, Temperature)
        # Plot the spectrum

        ratio = spectrum_P/spectrum_S
        SP[i,:] = spectrum_P 
        SS[i,:] = spectrum_S 
        RAT[i,:] = ratio

        if NWavelength > 1:
            ratio_plotter(Wave_list, spectrum_S, spectrum_P, ratio, id, Theta)

    t12 = time.time()
    # save RAT to temp/ folder
    Ratio = sym_complete(RAT, 0)  # 补全对称部分
    Ratio = Ratio.T
    Spectrum_S = sym_complete(SS, 0)
    Spectrum_S = Spectrum_S.T
    Th_list = sym_complete(Theta_list, 0)
    Th_list[Ntheta: 2*Ntheta] = 2*np.pi - Th_list[Ntheta: 2*Ntheta]

    np.save(f'temp/R{id}/Results/Ratio.npy', Ratio)
    np.save(f'temp/R{id}/Results/Star_flux.npy',  Spectrum_S)
    np.save(f'temp/R{id}/Results/Theta.npy', Th_list)

    print(t11-t10, ' s')
    print(t12-t11, ' s')


    ## 写一个自动画Ratio- Th_list图的程序

    # print(RAT)
    # if Ntheta > 2 :   # Ntheta = 1, 2 too less to plot; the main intention is to plot the contrast ratio
    #     x = np.linspace(0, 2*np.pi, 400)
    #     f = interp1d(Theta_list2, ratio, kind='cubic')
    #     y = f(x)

    #     plt.figure(figsize=(8, 8))
    #     plt.plot(x , y * 1e6 , 'k-')
    #     plt.xlabel('Orbital Phase Angle (rad)')
    #     plt.ylabel('Contrast Ratio (ppm)')
    #     plt.title('LHS 3844 b')
    #     plt.savefig(f'temp/R{id}/Results/contrast_ratio.png')
    #     plt.close()
    #     np.save(f'temp/R{id}/Results/ratio.npy', ratio)



    ##写一个自动画ratio图的程序

def ratio_plotter(Wavelength, spectrum_S, spectrum_P, ratio, id, Theta):
            
    #默认字体
    # config = {
    #     "font.family":'Times New Roman',
    #     "font.size": 16,
    #     "mathtext.fontset":'stix',
    #     "font.serif": ['Times New Roman'],
    # }
    # rcParams.update(config)
    #mpl.rcParams['font.sans-serif'] = ['SimHei']   #显示中文
    mpl.rcParams['axes.unicode_minus']=False       #显示负号

    a = 0.45
    fig, ax = plt.subplots(figsize=(26*a,13*a))
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.3)
    plt.rcParams['ytick.direction'] = 'in'# 刻度线显示在内部
    plt.rcParams['xtick.direction'] = 'in'# 刻度线显示在内部

    axpos = [0.1, 0.15, 0.7, 0.7]
    bwith = 2
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.set_position(axpos)
    #ax.axhline(y=np.average(Tc[3:]), color='gray', ls='-', )
    ax.plot(Wavelength* 1e6, spectrum_S, 'k-')
    ax.set_ylim(ymin=0, ymax= np.max(spectrum_S)*1.1)
    ax.set_xlabel('$\mathrm{Wavelength \; (\mu m)}$', fontsize=18)
    ax.set_ylabel('$\mathrm{Spectrum\; of\; Star \; (W \cdot sr^{-1}\cdot nm^{-1})}$', fontsize=18)
    ax.tick_params(length=6, width=2)
    ax.spines['right'].set_visible(False)

    lambda_color = 'blue'
    labmda_ax = ax.twinx()
    labmda_ax.set_position(axpos)
    labmda_ax.plot(Wavelength* 1e6, spectrum_P, '--', color=lambda_color)
    labmda_ax.set_ylim(ymin=0, ymax= np.max(spectrum_P)*1.1)
    labmda_ax.set_xlabel('$\mathrm{Wavelength \; (\mu m)}$', fontsize=18)
    labmda_ax.tick_params(length=6, width=2, color=lambda_color, labelcolor=lambda_color)
    labmda_ax.set_ylabel('$\mathrm{Spectrum\; of\; Planet \; (W \cdot sr^{-1}\cdot nm^{-1})}$', fontsize=18, color=lambda_color)
    labmda_ax.spines['right'].set(color=lambda_color, linewidth=2.0, linestyle=':')

    omglog_color = 'red'
    omglog_ax = ax.twinx()
    # 使用科学计数法的刻度
    omglog_ax.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
    # 获取 y 轴 OffsetText 对象
    offset_text = omglog_ax.yaxis.get_offset_text()
    # 调整位置示例，偏移 (dx, dy) 单位是像素 (points)
    offset_text.set_position((1.12, 0))
    # 调整字体大小
    offset_text.set_size(18)  # 或者使用 offset_text.set_fontsize(12)
    omglog_ax.spines['right'].set_position(('data', np.max(Wavelength)*1.15e6))
    omglog_ax.set_ylim(0, np.max(ratio)*1.1)
    omglog_ax.set_position(axpos)
    omglog_ax.plot(Wavelength* 1e6, ratio, '-.', color=omglog_color)
    omglog_ax.set_ylabel('$\mathrm{Relative\; spectrum\; of\; planet}$', fontsize=18, color=omglog_color)
    omglog_ax.tick_params(length=6, width=2, color=omglog_color, labelcolor=omglog_color)
    omglog_ax.spines['right'].set(color=omglog_color, linewidth=2.0, linestyle='-.')


    plt.title(f'$\mathrm{{\Theta = {Theta}}}$', fontsize=18)
    # save the plot to temp/ folder
    os.makedirs(f'temp/R{id}/Results', exist_ok=True)
    name = f'temp/R{id}/Results/spectrum_'+str(int(Theta*180/np.pi))+'.png'
    plt.savefig(name)
    plt.close()

@decorator_timer('vertify_radiation')
def vertify_radiation(wavelength_bound, Temperature= Temperature, id=0, Ntheta = 5, NWavelength = 1):
    # Calculate the blackbody radiation spectrum
    os.makedirs(f'temp/V{id}/plots', exist_ok=True)
    os.makedirs(f'temp/V{id}/Results', exist_ok=True)
    if Ntheta == 1:
        Theta_list = np.array([np.pi])
    else:
        Theta_list = np.linspace(0, np.pi, Ntheta)  # 0-pi 与 pi-2pi 重复
    Wavelength = np.linspace(wavelength_bound[0], wavelength_bound[1], NWavelength)

    RAT = np.zeros([len(Theta_list), len(Wavelength)])

    for i, Theta in enumerate(Theta_list):
        for j, wavelength in enumerate(Wavelength):
            # Calculate the blackbody radiation spectrum
            if Theta < 1e-6:
                RAT[i, j] = 0
            else:
                RAT[i, j] = para_rad(Theta, wavelength, Temperature)
        # Plot the spectrum

    #save RAT to temp/ folder
    np.save(f'temp/V{id}/Results/RAT.npy', RAT)
    ratio = np.zeros(len(Theta_list)*2 - 1)

    print(RAT)
    if Ntheta > 2:   # Ntheta = 1, 2 too less to plot; the main intention is to plot the contrast ratio
        ratio[0:len(Theta_list)] = RAT[:,0].flatten()
        ratio[len(Theta_list)-1:] = RAT[::-1].flatten()
        Theta_list2 = np.linspace(0, 2*np.pi, 2*len(Theta_list)-1)
        x = np.linspace(0, 2*np.pi, 400)
        f = interp1d(Theta_list2, ratio, kind='cubic')
        y = f(x)

        plt.figure(figsize=(8, 8))
        plt.plot(x , y * 1e6 , 'k-')
        plt.xlabel('Orbital Phase Angle (rad)')
        plt.ylabel('Contrast Ratio (ppm)')
        plt.title('LHS 3844 b')
        plt.savefig(f'temp/V{id}/Results/contrast_ratio.png')
        plt.close()
        np.save(f'temp/V{id}/Results/ratio.npy', ratio)

    Tmap = np.zeros((SIZE[0], SIZE[1]))
    phiP_list = np.linspace(0, 2* np.pi, SIZE[1])
    thetaP_list = np.linspace(0, np.pi, SIZE[0])
    r = orbit_calculator(a, e, Theta)

    def Tp(cos_Psi):
        #calculate the temperature of planet surface
        if cos_Psi > 0:
            return  Temperature* np.sqrt(R1 / r) * cos_Psi**(1/4)
        else:
            return 0

    for i, thetaP in enumerate(thetaP_list):
        for j, phiP in enumerate(phiP_list):
            cos_Psi = np.sin(thetaP) * np.cos(phiP)
            Tmap[i, j] = Tp(cos_Psi)

    # plot Tmap, xlabel: "$\theta$", ylabel: "$\phi$", using the gray map to show the temperature distribution
    plt.imshow(Tmap, cmap='gray')
    plt.xlabel('$\phi$')
    plt.ylabel('$\theta$')
    plt.colorbar()
    plt.savefig(f'temp/V{id}/Results/Tmap{int(Theta*180/np.pi)}.png')
    plt.close()
    np.save(f'temp/V{id}/Results/Tmap{int(Theta*180/np.pi)}.npy', Tmap)


