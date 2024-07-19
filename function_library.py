import numpy as np
from parameter_list import *
from scipy.integrate import dblquad, quad
from scipy import interpolate
from scipy.optimize import root  
import matplotlib.pyplot as plt
import os


def orbit_calculator(a, e, Theta):
    """
    Calculate the distance r from the Sun to the Earth at a given orbital angle Theta.
    
    Parameters:
    a (float): Semi-major axis of the Earth's orbit.
    e (float): Eccentricity of the Earth's orbit.
    Theta (float): Orbital angle.
    
    Returns:
    float: Distance r.
    """
    return a * (1 - e**2) / (1 + e * np.cos(Theta))


def I_dist(R1, r, angle):
    """
    Calculate the initial intensity of sunlight at a given angle on the Earth's surface.
    
    Parameters:
    R1 (float): Radius of the Sun.
    r (float): Distance from the Sun to the Earth.
    angle (float): Angle of incidence.
    
    Returns:
    float: Intensity I.
    """
    C = 1  # Coefficient of proportionality
    I = C * 2 * np.pi * (r * np.sin(2 * angle) / 2 + np.sin(angle) * np.sqrt(R1**2 - r**2 * np.sin(angle)**2))
    return I


def check_normalization(vec):
    """
    Normalize the given vector if it is not already normalized.
    
    Parameters:
    vec (array): Input vector.
    
    Returns:
    array: Normalized vector.
    """
    norm = np.linalg.norm(vec)
    return vec / norm if norm != 1 else vec


def normal_vec(phiP, thetaP ,Theta ,a, e, R2):
    """
    Calculate the normal vector on the Earth's surface and its position.
    
    Parameters:
    phiP (float): Latitude angle.
    thetaP (float): Longitude angle.
    Theta (float): Orbital angle.
    
    Returns:
    tuple: Normal vector, position vector, and distance r.
    """
    r = orbit_calculator(a, e, Theta)

    xP = R2 * np.cos(phiP) * np.cos(thetaP)
    yP = R2 * np.cos(phiP) * np.sin(thetaP)
    zP = R2 * np.sin(phiP)
    x = r * np.cos(Theta) + xP
    y = r * np.sin(Theta) + yP
    z = zP
    nv = np.array([xP, yP, zP])
    Pos = np.array([x, y, z])
    return check_normalization(nv), Pos, r


def reflect(Pos, normal):
    """
    Calculate the reflected vector.
    
    Parameters:
    Pos (array): Incident vector. (Unnormalized)
    normal (array): Normal vector at the point of reflection.
    
    Returns:
    array: Reflected vector.
    """
    rv = Pos - 2 * np.dot(Pos, normal) * normal
    return check_normalization(rv)


def angle_between(v1, v2):
    """
    Calculate the angle between two vectors.
    
    Parameters:
    v1 (array): First vector.
    v2 (array): Second vector.
    
    Returns:
    float: Angle between the vectors.
    """
    return np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))

def vec2Euler(vec):
    """
    Convert a vector to Euler angles.
    
    Parameters:
    vec (array): Input vector.
    
    Returns:
    array: Euler angles.
    """
    x, y, z = vec
    if z == 0:
        theta = np.pi / 2
    else:
        theta = np.arctan(np.sqrt(x**2 + y**2)/z)

    if x == 0:
        phi = np.pi/2
    else:
        phi = np.arctan(y/x)
    return  theta, phi

def Euler2vec(theta, phi):
    """
    Convert Euler angles to a vector.
    
    Parameters:
    theta (float): Theta angle.
    phi (float): Phi angle.
    
    Returns:
    array: Vector.
    """
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return np.array([x, y, z])


def check_direction(RV, normal, camera, Pos):
    """
    Check if the reflection direction is towards the camera.
    
    Parameters:
    RV (array): Reflected vector.
    normal (array): Normal vector at the point of reflection.
    camera (array): Camera/observation vector.
    Pos (array): Position vector.
    
    Returns:
    bool: True if the reflection direction is towards the camera, False otherwise.
    """
    i = angle_between(RV, normal)
    j = angle_between(camera, normal)
    #k = angle_between(camera, RV)
    return j < np.pi / 2 and np.dot(Pos, normal) < 0  and 0 <= i < np.pi / 2

def check_intersection_with_star(Pos,Camera):
    """
    Check if the line intersects with the star.
    """
    return np.linalg.norm(np.cross(Pos, Camera))/np.linalg.norm(Camera) < R1

def sym_complete(Var,axis=0):
    """
    Complete the symmetric part of the 2D array.
    
    Parameters:
    Var (array): Input 2D array.
    axis (int): Axis of symmetry.
    
    Returns:
    array: Symmetric 2D array.
    """
    print(Var.ndim)
    if Var.ndim == 1:
        return np.hstack((Var, Var[::-1]))
    
    if axis == 0:
        return np.vstack((Var, Var[::-1]))
    elif axis == 1:
        return np.hstack((Var, Var[:, ::-1]))
    else:
        raise ValueError("The axis of symmetry must be 0 or 1")
    


def blackbody_radiation(T, lam, B=1):
    """
    Calculate the blackbody radiation intensity at a given temperature and wavelength.
    !!!Totally same with function: B(lam,T)!!!
    
    Parameters:
    T (float): Temperature.
    lam (float): Wavelength.
    
    Returns:
    float: Intensity.
    """
    h = 6.626e-34  # Planck's constant
    c = 3.0e8  # Speed of light
    k = 1.38e-23  # Boltzmann constant
    return  2* B * h * c**2 / lam**5 / (np.exp(h * c / lam / k / T) - 1)

#print(blackbody_radiation(6000, 1e-6))


# def I_dist_afterreflect(R1, r, angle, normal, RV, Pos): #old version
#     """
#     Calculate the intensity of reflected sunlight considering reflectivity and diffusion.
    
#     Parameters:
#     R1 (float): Radius of the Sun.
#     r (float): Distance from the Sun to the Earth.
#     angle (float): Angle of incidence.
#     normal (array): Normal vector at the point of reflection.
#     RV (array): Reflected vector.
#     Pos (array): Position vector.
    
#     Returns:
#     float: Intensity I after reflection.
#     """
#     C = 1  # Coefficient of proportionality
#     R = 0.1  # Reflectivity of the planet
#     M = 0.2  # Diffusion coefficient of the planet

#     # Ensure the angle is within the valid range  {angle = angle_between(camera, RV)}
#     # if angle > np.pi / 2:
#     #     raise ValueError("The angle is larger than pi/2")

#     # Calculate the total intensity
#     S_hat = 1     #2 * np.pi * R1**2 / r * np.sqrt(r**2 - R1**2) * (1 - R1 / r)
#     I_tot = C * S_hat/r**2

#     # Calculate the diffused intensity
#     ID = I_tot * M * max(0, np.dot(normal, RV))

#     if check_intersection_with_star(Pos, camera):
#         return 0
    
#     # Calculate the final intensity considering the reflectivity
#     if r * np.sin(angle) < R1 and angle < np.pi/2:
#         cosbeta = np.sqrt(1 - (r/R1*np.sin(angle))**2)
#         sinS = (r/R1 * np.sin(2 * angle) / 2 + np.sin(angle) * cosbeta)
#         if angle == 0:
#             dis = r - R1
#         else:
#             dis = R1/np.sin(angle)*sinS

#         I =  cosbeta /dis**2 + ID
#         #print(cosbeta /dis**2,',',ID)
#     else:
#         I = ID

    
#     return I *a**2 *2
def Wave_reflect(R1, r, normal, Pos, camera ):
    """
    Calculate the intensity of reflected sunlight considering reflectivity and diffusion.
    
    Parameters:
    R1 (float): Radius of the Sun.
    r (float): Distance from the Sun to the Earth.
    angle (float): Angle of incidence.
    normal (array): Normal vector at the point of reflection.
    RV (array): Reflected vector.
    Pos (array): Position vector.
    
    Returns:
    float: Intensity I after reflection divided by B(T,lam).
    """
    #change the coordinate system from global to local

    uz = normal
    uy = check_normalization(np.cross(np.array([0,0,1]), uz))
    ux = np.cross(uy, uz)

    cx = np.dot(camera, ux)
    cy = np.dot(camera, uy)
    cz = np.dot(camera, uz)
    camera_local = np.array([cx, cy, cz])

    Px = np.dot(Pos, ux)
    Py = np.dot(Pos, uy)
    Pz = np.dot(Pos, uz)
    Pos_local = np.array([Px, Py, Pz])

    theta_c = angle_between(camera, normal)
    t , phi_camera = vec2Euler(camera_local)
    theta_P, phi_P = vec2Euler(Pos_local)

    
    def integrate_func(theta_i, phi):

        tvec = Euler2vec(theta_i, phi)
        tvec = rotate_vector(tvec, 'y', theta_P)
        tvec = rotate_vector(tvec, 'z', phi_P)
        

        theta_i, phi = vec2Euler(tvec)

        #theta_i = np.pi - angle_between(Pos, normal)
        phi_diff = phi - phi_camera
        mWi = tvec
        angle = angle_between(mWi, -Pos_local)
        angle_max = np.arcsin(R1/r)

        if angle > angle_max:   #超出star的范围
            return 0
       
        theta_i = angle_between(mWi, camera_local)/2

        ref_normal = check_normalization(check_normalization(mWi) + check_normalization(camera_local))   #mid vector
        titl = angle_between(ref_normal, normal)   #the titl of surface reflect needed

        return wave_dist(titl)* np.sin(theta_i) *Fresnel(theta_i, N1, N2)
    # np.sin(theta_i) is from the derivative of solid angle

    theta_max = np.arcsin(R1/r)
    Integ = dblquad(integrate_func, 0, 2*np.pi, 0, theta_max )
    #print(Integ)
    Dtheta = np.pi/SIZE[0]
    Dphi = 2*np.pi/SIZE[1]
    theta = angle_between(normal, np.array([0,0,1]))
    DA = R2**2 *np.sin(theta) *Dtheta *Dphi 
    return Integ[0] *DA *np.cos(theta_c) #* blackbody_radiation(Temperature, Wavelength) 



def Oren_Nayar_BRDF(R1, r, normal, Pos, camera, Coarse = 0):
    """
    Calculate the intensity of  diffusion.
    
    Parameters:
    R1 (float): Radius of the Sun.
    r (float): Distance from the Sun to the Earth.
    angle (float): Angle of incidence.
    normal (array): Normal vector at the point of reflection.
    RV (array): Reflected vector.
    Pos (array): Position vector.
    Coarse (float): standard deviation of the Gaussian distribution of the surface inclination (0-pi/2)
    
    Returns:
    float: the intensity of the diffusion divided by the B_s(T,lam)
    """
    #change the coordinate system from global to local
    uz = normal
    uy = check_normalization(np.cross(np.array([0,0,1]), uz))
    ux = np.cross(uy, uz)

    cx = np.dot(camera, ux)
    cy = np.dot(camera, uy)
    cz = np.dot(camera, uz)
    camera_local = np.array([cx, cy, cz])

    Px = np.dot(Pos, ux)
    Py = np.dot(Pos, uy)
    Pz = np.dot(Pos, uz)
    Pos_local = np.array([Px, Py, Pz])

    theta_c = angle_between(camera, normal)  #the angle between the camera and the normal vector
    t , phi_camera = vec2Euler(camera_local)


    def fr(theta_i, theta_c, sigma, phi_diff):
        #https://zhuanlan.zhihu.com/p/500809166
        A = 1 - 0.5 * sigma**2 / (sigma**2 + 0.33)
        B = 0.45 * sigma**2 / (sigma**2 + 0.09)
        alpha = max(theta_i, theta_c)
        beta = min(theta_i, theta_c)

        max_cosphi = max(0, np.cos(phi_diff))
        fr = 1 / np.pi *(A + B * max_cosphi * np.sin(alpha) * np.tan(beta))
        return fr
    
    def integrate_func(theta_i, phi):
        #theta_i = np.pi - angle_between(Pos, normal)
        phi_diff = phi - phi_camera
        mWi = Euler2vec(theta_i, phi)
        angle = angle_between(mWi, -Pos_local)
        angle_max = np.arcsin(R1/r)

        if angle > angle_max:
            return 0
        else:
            return fr(theta_i, theta_c, Coarse, phi_diff) * np.cos(theta_i) *np.sin(theta_i)
        #Bug Repaired in 7/2: np.sin(theta_i) -> np.cos(theta_i) *np.sin(theta_i)
        # the first cos is from the Equation, the second sin is from the expression of solid angle dOmega = sin(theta)*dTheta*dPhi
            
        #phi_diff = np.pi - angle_between(np.cross(Pos, normal), np.cross(camera, normal))
        # res = fr(theta_i, theta_c, sigma, phi_diff) * blackbody_radiation(6000, 1e-6) * np.sin(theta_i)
        # return res

    #hemi-sphere integral *  #put the blackbody radiation here(out of the integral)
    # def INT(func, phi_min, phi_max, theta_min, theta_max, epsabs=1e-4):
    #     theta_list = np.linspace(theta_min, theta_max, 30)
    #     phi_list = np.linspace(phi_min, phi_max, 30)
    #     res = 0
    #     for theta in theta_list:
    #         for phi in phi_list:
    #             res += func(theta, phi)

    #     return res
    
    theta_max = np.arcsin(R1/r)
    Integ = dblquad(integrate_func, 0, 2*np.pi,0, theta_max ,epsrel=1e-2,epsabs=1e-3)
    #print(Integ)
    Dtheta = np.pi/SIZE[0]
    Dphi = 2*np.pi/SIZE[1]
    theta = angle_between(normal, np.array([0,0,1]))
    DA = R2**2 *np.sin(theta) *Dtheta *Dphi 
    return Integ[0] *DA *np.cos(theta_c) #* blackbody_radiation(Temperature, Wavelength)   


def specular_reflection(RV, camera, normal, r, Temperature= Temperature):
    """
    Calculate the intensity of specular reflection.
    
    Parameters:
    RV (array): Reflected vector.
    camera (array): Camera/observation vector.
    normal (array): Normal vector at the point of reflection.
    
    Returns:
    float: Intensity I after reflection divided by B(T,lam).
    """
    # Calculate the angle between the camera and the reflected vector
    angle = angle_between(camera, RV)
    # Ensure the angle is within the valid range
    angle_max = np.arcsin(R1/r)
    if angle > angle_max:
        return 0
    
    theta_c = angle_between(camera, normal)
    theta = angle_between(normal, np.array([0,0,1]))
    # Calculate the intensity of the reflected light
    Dtheta = np.pi/SIZE[0]
    Dphi = 2*np.pi/SIZE[1]
    DA = R2**2 *np.sin(theta)*Dtheta*Dphi
    if specular_coefficent > 1:
        specular_coefficent = Fresnel(theta_c , N1, N2)        # Change here for different model of the reflection

    if specular_coefficent < 0:
        specular_coefficent = REF_fit(theta_c)
    
    return specular_coefficent * DA * np.cos(theta_c) #* blackbody_radiation(Temperature, Wavelength)
   #Bug repaired in 7/1: * np.cos(theta_c)


def REF_fit(theta_r):
    """
    Using experimental data (R-theta_r) to fit the reflection coefficient.
    theta_r: the angle between the camera and the reflected vector (unit: rad).
    """
# 将数据导入 numpy.array
    data = np.array([
        [0.5001124098736187, 0.032118887679963626],
        [0.5736634148850275, 0.033408595176071754],
        [0.6427189168144851, 0.06128842237195531],
        [0.6756905530154491, 0.061634824962347956],
        [0.7072081925499654, 0.14880020720629672],
        [0.7377997362818145, 0.07614928554558476],
        [0.7657341464915512, 0.05993461959404889],
        [0.8190099736814429, 0.0539751987302266],
        [0.8658608940248835, 0.08729094197427156],
        [0.9063547987087033, 0.128163086839302],
        [0.939014585557566, 0.2742055067313157],
        [0.9655296573641872, 0.32861747980467393],
        [0.9843427928509834, 0.42644452403552857],
        [0.9956399613573477, 0.4808017640592614],
        [0.9991605507473171, 0.6134471897978909]
    ])

    # 拆分数据为 x 和 y
    x = data[:, 0] 
    y = data[:, 1] 

    # 使用 scipy 的线性插值Linear拟合数据
    spl = interpolate.interp1d(x, y, kind='linear')
    
    xq = np.cos(theta_r)
    if xq < 0.5:
        yq = 0.032 
    else:
        yq = spl(xq)
    # 计算对应的 y 值
    return yq


# 按列绘制TOT_Intensity,x轴为Theta_list,绘制在同一副图中
def result_ploter(Var, name, Theta_list, Coarse, id):
    plt.figure()

    plt.plot(Theta_list, Var[:], label = "Coarse = " + str(int(Coarse*180/np.pi)))
    plt.xlabel('Orbit angle')
    plt.ylabel(name)
    plt.legend()

    plt.savefig(f'temp/R{id}/Results/{name}.png')
    plt.close()

def multi_result_plotter(Var, name, Theta_list, Coarse, id):
    plt.figure()
    # 启用 LaTeX 渲染
    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')

    si = np.size(Var, 0)
    for i in range(si):
        plt.plot(Theta_list, Var[i, :], label = name[i])

    plt.legend()
    plt.xlabel('Orbit angle')
    plt.ylabel("Flux_planet/Flux_star")

    plt.savefig(f'temp/R{id}/Results/Result.png')
    plt.close()


def wave_dist(alpha):
    #the distribution of the wave on the ocean obey Gaussian distribution
    #  alpha is the angle between the wave and the normal vector of the ocean
    # CITE: https://arxiv.org/pdf/0801.1852
    #print("HHH")
    
    sigma2 = 0.003 + 0.00512* Wind_speed   #Wind_speed in parameter_list
    return 1/np.sqrt(2*np.pi*sigma2)*np.exp(- np.tan(alpha)**2 /(2*sigma2))


def Cal_intersection_area(d, r1, r2):
    """
    Calculate the area of the intersection of two circles.
    d : float, distance between the centers of the two circles
    r1 : float, radius of the first circle
    r2 : float, radius of the second circle
    """
    if d >= (r1 + r2):  # 两圆相离
        return 0
    if (r1 - r2) >= d:  # 两圆内含，r1 大
        return np.pi * r2 ** 2
    if (r2 - r1) >= d:  # 两圆内含，r2 大
        return np.pi * r1 ** 2
    
    angle1 = np.arccos((r1 ** 2 + d ** 2 - r2 ** 2) / (2 * r1 * d))
    angle2 = np.arccos((r2 ** 2 + d ** 2 - r1 ** 2) / (2 * r2 * d))
    
    s1 = angle1 * r1 ** 2
    s2 = angle2 * r2 ** 2
    s3 = r1 * d * np.sin(angle1)
    s = s1 + s2 - s3
    
    return s


def Cal_star_area(Theta):
    # Calculate the area of the star that radiates light to the Earth, considering the blocking effect of the planet
    r = orbit_calculator(a, e, Theta) # Distance from the Sun to the Earth
    Planet = np.array([r * np.cos(Theta), r * np.sin(Theta), 0])  # Position vector of the planet

    Area = np.pi * R1 ** 2
    if np.dot(camera,Planet) < 0:  # Planet is behind the star, impossible to block the light
        return Area
    else:
        d = np.linalg.norm(np.cross(Planet, camera))
        return Area - Cal_intersection_area(d, R1, R2) 


def Cal_star_flux(Theta, Wavelength = Wavelengh, Temperature = Temperature):
    # Calculate the flux of the star that radiates light to the Earth
    #检查Theta的数据类型是数还是数组？
    Si = np.size(Theta)
    if Si == 1:
        Area = Cal_star_area(Theta)
        return Area * B(Wavelength, Temperature)
    
    Flux = np.zeros(Si)
    
    for i, Th in enumerate(Theta):
        Area = Cal_star_area(Th)
        Flux[i] = Area 

    return Flux * B(Wavelength, Temperature)


def Fresnel(theta_i, n1, n2):
    # Calculate the Fresnel coefficient
    # n1, n2: refractive index of the two media
    # theta_i: angle of incidence
    # no polarized light
    # Reference: https://en.wikipedia.org/wiki/Fresnel_equations

    theta_t = np.arcsin(n1/n2 * np.sin(theta_i))
    Rs = ((n1 * np.cos(theta_i) - n2 * np.cos(theta_t)) / (n1 * np.cos(theta_i) + n2 * np.cos(theta_t)) )**2
    Rp = ((n1 * np.cos(theta_t) - n2 * np.cos(theta_i)) / (n1 * np.cos(theta_t) + n2 * np.cos(theta_i)) )**2
    Reff = (Rs + Rp) / 2

    return Reff

def rotate_vector(vector, axis, angle):
    """
    旋转向量。
    
    Parameters:
    vector (array-like): 待旋转的向量，例如 [x, y, z]
    axis (str): 'x', 'y', 或 'z' 指定绕哪个轴旋转
    angle (float): 旋转角度（以弧度为单位）
    
    Returns:
    np.ndarray: 旋转后的向量
    """
    # 定义旋转矩阵
    if axis == 'x':
        R = np.array([[1, 0, 0],
                      [0, np.cos(angle), -np.sin(angle)],
                      [0, np.sin(angle), np.cos(angle)]])
    elif axis == 'y':
        R = np.array([[np.cos(angle), 0, np.sin(angle)],
                      [0, 1, 0],
                      [-np.sin(angle), 0, np.cos(angle)]])
    elif axis == 'z':
        R = np.array([[np.cos(angle), -np.sin(angle), 0],
                      [np.sin(angle), np.cos(angle), 0],
                      [0, 0, 1]])
    else:
        raise ValueError("轴必须是 'x', 'y' 或 'z'")

    # 旋转向量
    rotated_vector = R @ np.array(vector)
    return rotated_vector

# # 示例使用
# if __name__ == "__main__":
#     vector = [1, 0, 0]  # 待旋转的向量
#     axis = 'z'  # 旋转轴
#     angle = np.pi / 2  # 旋转90度
#     rotated_vector = rotate_vector(vector, axis, angle)
#     print("原向量:", vector)
#     print("旋转后的向量:", rotated_vector)

def B(lam,T):
    # !!!Totally same with function: blackbody_radiation(T,lam)!!!
    h = 6.626e-34  # Planck's constant
    c = 3.0e8  # Speed of light
    k = 1.38e-23  # Boltzmann constant
    A = (np.exp(h * c / lam / k / T) - 1)
    B = 2* h * c**2 / lam**5 / A
    return B

def Temperature_cal(ksi, Theta, Int_B):
    ## calculate the temperature distribution of the planet
    ## ksi is the angle between the normal vector and the vector from the star to the planet
    r = orbit_calculator(a, e, Theta)

    ksi_m1 = np.pi/2 - np.arcsin((R1+R2)/r)
    ksi_m2 = np.pi/2 + np.arcsin((R1-R2)/r)
    if ksi_m1 > ksi_m2:
        #报错信息： ksi_m1 > ksi_m2
        print("ksi_m1 > ksi_m2")
        raise ValueError("ksi_m1 > ksi_m2")
    
    if ksi <= ksi_m2:
        T0 = Temperature * np.sqrt(R1/r)

    else:
        return 0

    if ksi < ksi_m1:
        rP = np.sqrt(r**2 + R2**2 - 2*r*R2*np.cos(ksi))
        zeta = np.arcsin(R2/rP*np.sin(ksi))
        Phi = zeta + ksi

        LHS = Int_B *(R1/r)**2 *np.cos(Phi)
        def equation(T):
            func = lambda lam: B(lam, T) * (1-Albedo(lam))
            return quad(func, 0, np.inf)[0] - LHS
        
        sol = root(equation, T0)
        T = sol.x[0]

    elif ksi_m1 <= ksi <= ksi_m2: # the planet is in the shadow of the star

        def integrate_func(thetas, phis): # the integral function
            normalP = np.array([-R2* np.cos(ksi), 0, R2* np.sin(ksi)])
            normalS = np.array([np.sin(thetas)*np.cos(phis), np.sin(thetas)*np.sin(phis), np.cos(thetas)])* R1
            P = np.array([r, 0 ,0])
            Pos = P + normalP - normalS

            TH = angle_between(Pos, normalS)
            PSI = angle_between(Pos, -normalP)

            if TH < np.pi/2 and PSI < np.pi/2:
                I = (R1/np.linalg.norm(Pos))**2 * np.sin(thetas) *np.cos(TH) * np.cos(PSI)
            else:
                I = 0

            return I
        
        Int = dblquad(integrate_func, 0, 2* np.pi, 0, np.pi)
        LHS = Int[0] * Int_B/ np.pi

        def equation(T): # the integral function
            
            func = lambda lam: B(lam, T) * (1-Albedo(lam))
            return quad(func, 0, np.inf)[0] - LHS

        sol = root(equation, T0)
        T = sol.x[0]

    else:
        T = 0

    print(T)
    return T

# def Temperature_cal(ksi, Theta, Albedo = 0, Temperature = Temperature):
#     ## calculate the temperature distribution of the planet
#     ## ksi is the angle between the normal vector and the vector from the star to the planet
#     r = orbit_calculator(a, e, Theta)

#     ksi_m1 = np.pi/2 - np.arcsin((R1+R2)/r)
#     ksi_m2 = np.pi/2 + np.arcsin((R1-R2)/r)
#     if ksi_m1 > ksi_m2:
#         #报错信息： ksi_m1 > ksi_m2
#         print("ksi_m1 > ksi_m2")
#         raise ValueError("ksi_m1 > ksi_m2")

#     if ksi < ksi_m1:
#         rP = np.sqrt(r**2 + R2**2 - 2*r*R2*np.cos(ksi))
#         zeta = np.arcsin(R2/rP*np.sin(ksi))
#         Phi = zeta + ksi
#         T = Temperature * np.sqrt(R1/r) *((1-Albedo)*np.cos(Phi))**(1/4)

#     elif ksi_m1 <= ksi <= ksi_m2: # the planet is in the shadow of the star

#         def integrate_func(thetas, phis): # the integral function
#             normalP = np.array([-R2* np.cos(ksi), 0, R2* np.sin(ksi)])
#             normalS = np.array([np.sin(thetas)*np.cos(phis), np.sin(thetas)*np.sin(phis), np.cos(thetas)])* R1
#             P = np.array([r, 0 ,0])
#             Pos = P + normalP - normalS

#             TH = angle_between(Pos, normalS)
#             PSI = angle_between(Pos, -normalP)

#             if TH < np.pi/2 and PSI < np.pi/2:
#                 I = (R1/np.linalg.norm(Pos))**2 * np.sin(thetas) *np.cos(TH) * np.cos(PSI)
#             else:
#                 I = 0

#             return I
        
#         Int = dblquad(integrate_func, 0, 2* np.pi, 0, np.pi)
#         T = Temperature * (Int[0] /np.pi* (1- Albedo))**(1/4)

#     else:
#         T = 0

#     return T

def Tmap(Theta, id = 0 ):
    ## calculate the temperature map of the planet
    ## the map is a 2D array of thetaP and phiP
    def int_func(lam):
        return B(lam, Temperature) * (1- Albedo(lam))

    Int_B = quad(int_func, 0 , np.inf)[0]

    Tmap_1D = np.zeros(SIZE[0])  #Use ksi as parameter to calculate the temperature map, according to the rotation symmetry, we can get the temperature map of the planet
    ksi_list = np.linspace(0, np.pi, SIZE[0])

    for i in range(SIZE[0]):
        Tmap_1D[i] = Temperature_cal(ksi_list[i], Theta, Int_B)
    ## interpolate the 1D array to 2D array using the spline interpolation
    spl = interpolate.interp1d(ksi_list , Tmap_1D, kind='linear') #spline interpolation

    Tmap = np.zeros((SIZE[0], SIZE[1]))
    phiP_list = np.linspace(0, 2* np.pi, SIZE[1])
    thetaP_list = np.linspace(0, np.pi, SIZE[0])

    for i, thetaP in enumerate(thetaP_list):
        for j, phiP in enumerate(phiP_list):
            ksi = np.arccos(np.sin(thetaP) * np.cos(phiP)) 
            Tmap[i, j] = spl(ksi)

    # plot Tmap, xlabel: "$\theta$", ylabel: "$\phi$", using the gray map to show the temperature distribution
    os.makedirs(f'temp/R{id}/plots', exist_ok=True)
    plt.imshow(Tmap, cmap='gray')
    plt.xlabel('$\phi$')
    plt.ylabel('$\theta$')
    plt.colorbar()
    plt.savefig(f'temp/R{id}/plots/Tmap{int(Theta*180/np.pi)}.png')
    plt.close()
    np.save(f'temp/R{id}/plots/Tmap{int(Theta*180/np.pi)}.npy', Tmap)

    return Tmap


def Radiation_cal(Tmap, Theta, camera, Temperature, Wavelength = 0):
    ## calculate the radiation distribution of the planet
    ## the map is a 2D array of thetaP and phiP
    Rad = 0
    Dtheta = np.pi/SIZE[0]
    Dphi = 2*np.pi/SIZE[1]
    thetaP_list = np.linspace(0, np.pi, SIZE[0])
    phiP_list = np.linspace(0, 2* np.pi, SIZE[1])

    for i, thetaP in enumerate(thetaP_list):
        for j, phiP in enumerate(phiP_list):
            normalP = np.array([-np.sin(thetaP)*np.cos(phiP), -np.sin(thetaP)*np.sin(phiP), np.cos(thetaP)])
            normalP = rotate_vector(normalP, 'z', Theta)  #calculate the normal vector of the planet in global coordinate system

            angle = angle_between(normalP, camera)
            if angle >= np.pi/2:
                continue
            else:
                T = Tmap[i, j]
                if T < 1e-2:  # the temperature is too low, the radiation is negligible
                    continue

                dA = R2**2 * np.sin(thetaP) * Dtheta * Dphi
                if Wavelength == 0:
                    sigma = 5.670373 * 1e-8
                    Rad += sigma * Temperature**4/np.pi * np.cos(angle) * dA
                else:
                    Rad += (1-Albedo(Wavelength)) * B(Wavelength, T) * np.cos(angle) * dA

    return Rad



def para_rad(Theta, lam = 0, Temperature = Temperature):
    r = orbit_calculator(a, e, Theta)
    h = 6.626e-34  # Planck's constant
    c_const = 3.0e8  # Speed of light
    k = 1.38e-23  # Boltzmann constant
    Co = h * c_const/ k

    def Tp(cos_Psi):
        #calculate the temperature of planet surface
        if cos_Psi > 0:
            return  Temperature* np.sqrt(R1 / r) * cos_Psi**(1/4)
        else:
            return 0
        
    def Bf(lam,T):
        # planck function for radiation spectrum calculation
        if T < 0.01 :  # T = 0 cause the division by zero
            return 0
        
        A = (np.exp(Co / lam / T) - 1)
        B = 1 / lam**5 / A
        return B
    
    def Fp_func(theta, phi):
        # Radiation integral function
        cos_Psi = np.cos(theta) * np.cos(phi)
        if cos_Psi < 0:
            return 0
        
        cos_Psi_obs = np.cos(phi) * np.cos(theta + Theta - np.pi)
        
        if lam == 0: #全谱辐射强度， 相当于对波长进行积分
            return Tp(cos_Psi)**4 * np.cos(phi) *cos_Psi_obs
        else:   #特定波长附近的辐射功率密度
            return Bf(lam, Tp(cos_Psi)) * np.cos(phi) *cos_Psi_obs
        
    
    Int = dblquad(Fp_func,-np.pi/2, np.pi/2, np.pi/2-Theta, 3/2*np.pi- Theta)
    if lam == 0:
        Fp = Sigma * R2**2/ np.pi * Int[0]
        Fs = Sigma * Temperature**4 * np.pi* R1**2
    else:
        Fp = 2* h *c_const**2 * R2**2 * Int[0]
        Fs = Cal_star_flux(Theta, lam, Temperature)   #2* h *c_const**2 * Bf(lam, Temperature) * np.pi * R1**2

    return Fp/Fs


