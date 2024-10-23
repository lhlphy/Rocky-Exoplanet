import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import dblquad

AU = 1.496e8 # Astronomical unit (km)
RE = 6.371e3 # Earth radius (km)
RSun = 6.96e5 # Solar radius (km) 695,500 
semi_axis = 0.00709 * AU # semi-major axis of the orbit (km)
Rs = 0.458 * RSun    # 0.458 Solar radius (km)
Rp = 0.699 * RE    # 0.699 Earth raidus (km)

def int_func(theta, phi, Theta): 
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
    
def Specular_vertify_cal(NTheta):
    Theta_list = np.linspace(0, np.pi, NTheta)
    SR = np.zeros(NTheta)
    for i, Theta in enumerate(Theta_list):
        print(Theta)
        SR[i] = dblquad(int_func, np.pi/2,  3*np.pi/2, 0, np.pi, args = (Theta,), epsabs=1e-6, epsrel=1e-6)[0]
        
    SR = np.concatenate((SR, SR[::-1]))
    SR = SR *(- Rp**2 / Rs **2 / np.pi)
    Theta_list = np.concatenate((Theta_list, np.pi* 2 - Theta_list[::-1]))
    
    print(SR * 1e6)
    
    plt.plot(Theta_list/ (2* np.pi), SR * 1e6)
    plt.xlabel('Orbital Phase')
    plt.ylabel(r'$F_{Specular}$/$F_{*}$ (ppm)')
    plt.show()
            
if __name__ == '__main__':
    Specular_vertify_cal(10)
    
    