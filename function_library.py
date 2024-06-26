import numpy as np
from parameter_list import *
from scipy.integrate import dblquad
from scipy import interpolate



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

def blackbody_radiation(T, lam, B=1):
    """
    Calculate the blackbody radiation intensity at a given temperature and wavelength.
    
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


def Oren_Nayar_BRDF(R1, r, normal, Pos, camera, Coarse = 0, DIF_REF = 0.5, Temperature = 6000, Wavelengh = 1e-6):
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
    float: Intensity I after reflection.
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


    def fr(theta_i, theta_c, sigma, rho, phi_diff):
        #https://zhuanlan.zhihu.com/p/500809166
        A = 1 - 0.5 * sigma**2 / (sigma**2 + 0.33)
        B = 0.45 * sigma**2 / (sigma**2 + 0.09)
        alpha = max(theta_i, theta_c)
        beta = min(theta_i, theta_c)

        max_cosphi = max(0, np.cos(phi_diff))
        fr = rho / np.pi *(A + B * max_cosphi * np.sin(alpha) * np.tan(beta))
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
            return fr(theta_i, theta_c, Coarse, DIF_REF, phi_diff) * np.sin(theta_i)
            
        #phi_diff = np.pi - angle_between(np.cross(Pos, normal), np.cross(camera, normal))
        # res = fr(theta_i, theta_c, sigma, rho, phi_diff) * blackbody_radiation(6000, 1e-6) * np.sin(theta_i)
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
    

    Integ = dblquad(integrate_func, 0, 2*np.pi,0, np.pi/2 ,epsrel=1e-2,epsabs=1e-3)
    #print(Integ)
    Dtheta = np.pi/SIZE[0]
    Dphi = 2*np.pi/SIZE[1]
    theta = angle_between(normal, np.array([0,0,1]))
    return Integ[0] * R2**2 *np.sin(theta) *Dtheta *Dphi *np.cos(theta_c) #* blackbody_radiation(Temperature, Wavelengh)   


def specular_reflection(specular_coefficent ,RV, camera, normal, r, Temperature= 6000, Wavelengh = 1e-6):
    """
    Calculate the intensity of specular reflection.
    
    Parameters:
    RV (array): Reflected vector.
    camera (array): Camera/observation vector.
    normal (array): Normal vector at the point of reflection.
    
    Returns:
    float: Intensity I after reflection.
    """
    # Calculate the angle between the camera and the reflected vector
    angle = angle_between(camera, RV)
    # Ensure the angle is within the valid range
    angle_max = np.arcsin(R1/r)
    if angle > angle_max:
        return 0
    
    #theta_c = angle_between(camera, normal)
    theta = angle_between(normal, np.array([0,0,1]))
    # Calculate the intensity of the reflected light
    Dtheta = np.pi/SIZE[0]
    Dphi = 2*np.pi/SIZE[1]
    DA = R2**2 *np.sin(theta)*Dtheta*Dphi
    if specular_coefficent > 1:
        specular_coefficent = REF_fit(angle_between(normal,camera))
    
    return specular_coefficent * DA #* blackbody_radiation(Temperature, Wavelengh)


def REF_fit(theta_r):
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


