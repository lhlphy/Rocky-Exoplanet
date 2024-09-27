import numpy as np
from scipy.integrate import simps, quad
from parameter_list import PPs

def integrate_spectrum_with_interpolation(wave_list, IS, Wavelength_min, Wavelength_max):

    # 进行插值，计算Wavelength_min和Wavelength_max处的光谱数据
    IS_min = np.interp(Wavelength_min, wave_list, IS)
    IS_max = np.interp(Wavelength_max, wave_list, IS)

    # 确保波长范围在给定的范围内
    mask = (wave_list > Wavelength_min) & (wave_list < Wavelength_max)
    
    # 获取在范围内的波长和光谱数据
    wave_list_range = wave_list[mask]
    IS_range = IS[mask]

    # 将插值的Wavelength_min和Wavelength_max以及对应的光谱数据添加到范围内的数据
    wave_list_range = np.concatenate(([Wavelength_min], wave_list_range, [Wavelength_max]))
    IS_range = np.concatenate(([IS_min], IS_range, [IS_max]))

    # 使用辛普森积分公式进行积分
    integral = simps(IS_range, wave_list_range)
    
    return integral


def bond_albedo_calculator(Wavelength_min, Wavelength_max, name , Nmark = -1):
    # load data, I_diffuse,I_specular is contrast ratio 
    I_diffuse = np.load(f'temp/{name}/variables/I_diffuse.npy')
    I_specular = np.load(f'temp/{name}/variables/I_specular.npy')
    Star_flux = np.load(f'temp/{name}/variables/Star_flux.npy')
    wave_list = np.load(f'temp/{name}/variables/wave_list.npy')
    Thermal = np.load(f'temp/{name}/variables/Thermal.npy')

    if Wavelength_max > wave_list[-1] or Wavelength_min < wave_list[0]:
        raise ValueError('Wavelength exceed the bounds! ')

    # convert I_diffuse,I_specular (contrast ratio) to absolute intensity
    N1 = Star_flux.shape[0]
    N2 = Star_flux.shape[1]
    if Nmark == -1:
        Nmark = N2//2
        
    ID = I_diffuse[:, Nmark]
    IS = I_specular[:, Nmark]
    Star_flux = Star_flux[:, Nmark]
    Thermal = Thermal[:, Nmark]

    ID = ID *Star_flux
    IS = IS *Star_flux
    IT = Thermal *Star_flux

    # integrate over the wavelength range
    TS = integrate_spectrum_with_interpolation(wave_list, IS, Wavelength_min, Wavelength_max) # specular
    TD = integrate_spectrum_with_interpolation(wave_list, ID, Wavelength_min, Wavelength_max) # diffuse
    TT = integrate_spectrum_with_interpolation(wave_list, IT, Wavelength_min, Wavelength_max) # Thermal radiation
    Tstar = integrate_spectrum_with_interpolation(wave_list, Star_flux, Wavelength_min, Wavelength_max) # Star radiation flux

    # Tin = quad(B, Wavelength_min, Wavelength_max, args=(PPs.Stellar_T,)) * np.pi * R2**2  # Energy flux that the planet received

    Spectral_Contrast_ratio_S = TS / Tstar
    Spectral_Contrast_ratio_D = TD / Tstar
    Spectral_Contrast_ratio_T = TT / Tstar

    print('Spectral Contrast ratio of Specular reflection in [',Wavelength_min, Wavelength_max,'] $\mu$m is:', Spectral_Contrast_ratio_S* 1e6,' ppm' )
    print('Spectral Contrast ratio of diffuse reflection in [',Wavelength_min, Wavelength_max,'] $\mu$m is:', Spectral_Contrast_ratio_D* 1e6, ' ppm' )
    print('Spectral Contrast ratio of Thermal radiation in [',Wavelength_min, Wavelength_max,'] $\mu$m is:', Spectral_Contrast_ratio_T* 1e6,' ppm' )
    
    return (Spectral_Contrast_ratio_S + Spectral_Contrast_ratio_T) *1e6, (Spectral_Contrast_ratio_D + Spectral_Contrast_ratio_T) *1e6

bond_albedo_calculator(3.9e-6, 4.1e-6, 'R5')



    


    


