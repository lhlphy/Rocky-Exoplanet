import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import dblquad

AU = 1.496e8 # Astronomical unit (km)
RE = 6.371e3 # Earth radius (km)
RSun = 6.96e5 # Solar radius (km) 695,500 
semi_axis = 0.00709 * AU # semi-major axis of the orbit (km)
Rs = 0.458 * RSun    # 0.458 Solar radius (km)
Rp = 0.699 * RE    # 0.699 Earth raidus (km)

def int_func_smooth(theta, phi, Theta):  # specular reflection with perfect smooth surface
    # theta integration from 0 to pi, phi integration from pi/2 to 3pi/2
    sint = np.sin(theta)
    sinT = np.sin(Theta)
    cosT = np.cos(Theta)
    sinp = np.sin(phi)
    cosp = np.cos(phi)
    cond1 = (sint * cosp * cosT + sint * sinp * sinT > 0) 
    cond2 = (cosT *(1 - 2*sint **2 *cosp **2) - 2* sint **2 * sinp
             * cosp *sinT - np.sqrt(1 - Rs **2/ semi_axis **2) > 0)
    
    # print(cond1, cond2)
    if cond1 and cond2:
        return sint **2 * cosp
    else:
        return 0
    
def P_slope(slope, sigma): # guassian distribution of the slope of the surface
    sigma = sigma /180 *np.pi  # sigma is the slope of the surface
    P = 1/(np.sqrt(2 * np.pi) * sigma) * np.exp(-(slope **2) / (2 * sigma **2))
    return P

def int_func_rough(theta, phi, Theta, sigma):  # specular reflection with rough surface
    # theta integration from 0 to pi, phi integration from pi/2 to 3pi/2
    # theta integration from 0 to pi, phi integration from pi/2 to 3pi/2
    sint = np.sin(theta)
    sinT = np.sin(Theta)
    cosT = np.cos(Theta)
    sinp = np.sin(phi)
    cosp = np.cos(phi)
    cond1 = (sint * cosp * cosT + sint * sinp * sinT > 0) 
    slope = np.arccos(cosT *(1 - 2 * sint **2 * cosp **2) - 2 * sint **2 * sinp * cosp * sinT)
    Ps = P_slope(slope, sigma)
    
    # print(cond1, cond2)
    if cond1 :
        return sint **2 * cosp * Ps
    else:
        return 0

    
def Specular_vertify_cal(NTheta, Nsigma):
    Theta_list = np.linspace(0, np.pi, NTheta)
    SR_smooth = np.zeros(NTheta) # specular reflection of the smooth surface
    
    if Nsigma != 0:
        roughness_list = np.linspace(0, 2, Nsigma+1)
        roughness_list = roughness_list[1:]  # remove roughness = 0, because it is the same as smooth surface, and will cause error in guassian distribution
        SR_rough = np.zeros((Nsigma, NTheta)) # specular reflection of the rough surface
    
    for i, Theta in enumerate(Theta_list):
        print(Theta)
        SR_smooth[i] = dblquad(int_func_smooth, np.pi/2,  3*np.pi/2, 0, np.pi, args = (Theta,), epsabs=1e-6, epsrel=1e-6)[0]
        if Nsigma != 0:
            for j, sigma in enumerate(roughness_list):       
                SR_rough[j][i] = dblquad(int_func_rough, np.pi/2,  3*np.pi/2, 0, np.pi, args = (Theta, sigma), epsabs=1e-6, epsrel=1e-6)[0]
        
    SR_smooth = np.concatenate((SR_smooth, SR_smooth[::-1])) # expand 0-pi to 0-2*pi
    SR_smooth = SR_smooth *(- Rp**2 / Rs **2 / np.pi) # convert to contrast ratio
    if Nsigma != 0:
        SR_rough = np.concatenate((SR_rough, SR_rough[:,::-1]), axis = 1)
        SR_rough = SR_rough *(- Rp**2 / Rs **2 / np.pi)
    Theta_list = np.concatenate((Theta_list, np.pi* 2 - Theta_list[::-1])) # expand 0-pi to 0-2*pi
    
    print("SR_smooth: ",SR_smooth * 1e6)
    plt.figure()
    # 绘制模拟数据和理论值对比图
    # sim_theta = np.load('sim_theta.npy')
    # sim_SR = np.load('sim_SR.npy')
    
    # plt.plot(sim_theta, sim_SR, label = 'Numercial result')
    
    plt.plot(Theta_list/ (2* np.pi), SR_smooth * 1e6, label = 'Analytical result for smooth surface')
    if Nsigma != 0:
        print("SR_rough: ",SR_rough * 1e6)
        for k, roughness in enumerate(roughness_list):
            plt.plot(Theta_list/ (2* np.pi), SR_rough[k] * 1e6, label = f'rough surface, slope={roughness:.2f}') 
    plt.xlabel('Orbital Phase')
    plt.ylabel(r'$F_{Specular}$/$F_{*}$ (ppm)')
    plt.legend()
    plt.show()
            
if __name__ == '__main__':
    Specular_vertify_cal(6, 10)
    
    ## compare the numerical result and analytical result
    # p1 = np.loadtxt('simed.txt', delimiter = ',')
    # p2 = np.loadtxt('simplified.txt', delimiter = ',')
    # p3 = np.loadtxt('Analysis.txt', delimiter = ',')
    # plt.plot(p1[:, 0] * 2*np.pi, p1[:, 1], label = 'Simulated')
    
    # plt.plot(p3[:, 0] * 2* np.pi, p3[:, 1], label = 'Analytical result') 
    # plt.plot(p2[:, 0], p2[:, 1], label = 'Simplified')
    # plt.xlabel('Orbital Phase')
    # plt.ylabel(r'$F_{Specular}$/$F_{*}$ (ppm)')
    # plt.legend()
    # plt.show()







    
    