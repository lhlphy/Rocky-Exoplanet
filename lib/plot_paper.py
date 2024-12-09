import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from parameter_list import PPs


def data_loader(name, Obs_wave):
    # load data
    I_specular = np.load(f'temp/{name}/variables/I_specular.npy')
    I_intensity = np.load(f'temp/{name}/variables/I_intensity.npy')
    I_diffuse = np.load(f'temp/{name}/variables/I_diffuse.npy')
    Theta_list = np.load(f'temp/{name}/variables/Theta.npy')
    Wave_list = np.load(f'temp/{name}/variables/wave_list.npy')
    Thermal = np.load(f'temp/{name}/variables/Thermal.npy')

    if type(Obs_wave) == float:
        Obs_wave = np.array([Obs_wave])

    # interpolate
    # print(Wave_list.shape, I_specular.T.shape)
    I_specular = np.nan_to_num(I_specular, nan=0) 
    I_intensity = np.nan_to_num(I_intensity, nan=0) 
    I_diffuse = np.nan_to_num(I_diffuse, nan=0) 
    spls = interp1d(Wave_list, I_specular.T, kind='linear', fill_value='extrapolate')
    spli = interp1d(Wave_list, I_intensity.T, kind='linear', fill_value='extrapolate')
    spld = interp1d(Wave_list, I_diffuse.T, kind='linear', fill_value='extrapolate')
    splt = interp1d(Wave_list, Thermal.T, kind='linear', fill_value='extrapolate')
    I_s = spls(Obs_wave).T
    I_i = spli(Obs_wave).T
    I_d = spld(Obs_wave).T
    I_t = splt(Obs_wave).T
    
    
    if np.size(Theta_list) == 2:
        print('Wavelength: ', Obs_wave[0]*1e6, 'um')
        print('Secondary eclipse depth difference: ',(I_d-I_s) *1e6,' ppm')
        return 

    theta = np.linspace(0, 2*np.pi, 200)
    Theta_list[np.size(Theta_list)//2] += 1e-10
    spls = interp1d(Theta_list, I_s, kind='linear', fill_value='extrapolate')
    spli = interp1d(Theta_list, I_i, kind='linear', fill_value='extrapolate')
    spld = interp1d(Theta_list, I_d, kind='linear', fill_value='extrapolate')
    splt = interp1d(Theta_list, I_t, kind='linear', fill_value='extrapolate')
    I_s = spls(theta)
    I_i = spli(theta)
    I_d = spld(theta)
    I_t = splt(theta)
    return I_s, I_i, I_d, I_t, theta
    
    
def specular_diffuse_plot(name_specular, name_diffuse, Obs_wave, transit = 'off'):
    # 绘制diffuse和specular的phase curve(*not* include thermal emission), 绘图中包含2曲线，
    # 分别是specular和diffuse的phase curve, 不区分表面材料性质的low与high albedo, 仅用于对比高分辨率的diffuse and specular_only model
    Is_specular, Ii_specular, Id_specular, It_specular, theta = data_loader(name_specular, Obs_wave)
    Is_diffuse, Ii_diffuse, Id_diffuse, It_diffuse, theta = data_loader(name_diffuse, Obs_wave)
    
    i = 0  # specular 与 diffuse与波长无关，仅取决于几何结构，所以i可以任意选取
    theta = theta/(2 *np.pi)

    fig, ax = plt.subplots(figsize=(9,6))
    up_lim = 0  # ylim的上限
    if transit == 'off':  # 不考虑transit的情况
        ax.plot(theta, Is_specular[i,:] *1e6, label='specular', color='b', linewidth=2)
        ax.plot(theta, Id_diffuse[i,:] *1e6, label='diffuse', color='k', linewidth=2)
        up_lim = np.max([np.max(Is_specular[i,:]), np.max(Id_diffuse[i,:])]) *1e6
    else: # 考虑transit的情况, transit修正存储在一个数组中 It_specular, It_diffuse
        ax.plot(theta, (Is_specular[i,:] + It_specular[i,:]) *1e6, label='specular', color='b', linewidth=2)
        ax.plot(theta, (Id_diffuse[i,:] + It_diffuse[i,:]) *1e6, label='diffuse', color='k', linewidth=2)
        up_lim = np.max([np.max(Is_specular[i,:] + It_specular[i,:]), np.max(Id_diffuse[i,:] + It_diffuse[i,:])]) *1e6
        
    ax.set_xlabel('Orbital phase', fontsize=18)
    ax.set_ylabel(r'$F_p/F_*$ (ppm)', fontsize=18)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, up_lim *1.1)
    ax.spines['bottom'].set_linewidth(2)    ###设置底部坐标轴的粗细
    ax.spines['left'].set_linewidth(2)  ####设置左边坐标轴的粗细
    ax.spines['right'].set_linewidth(2) ###设置右边坐标轴的粗细
    ax.spines['top'].set_linewidth(2)   ####设置上部坐标轴的粗细
    #刻度值字体大小设置（x轴和y轴同时设置）
    plt.tick_params(labelsize=16)
    plt.legend(fontsize=17, frameon=False)
    plt.savefig(f"temp/{name_specular}/specular_diffuse_{Obs_wave[0]*1e6}_{transit}.png")
    plt.savefig(f"temp/{name_diffuse}/specular_diffuse_{Obs_wave[0]*1e6}_{transit}.png")
    plt.close()
 
def analytical_theory_cal():
    Theta = np.linspace(0, 2*np.pi, 200)
    RESULT = np.zeros((200, 2))
    RESULT[:,0] = Theta
    
    alpha = np.arcsin(PPs.Rs / (PPs.semi_axis - PPs.Rp))
    deltaphi = alpha
    phi0 = (Theta + np.pi)/2
    theta0 = 0
    deltatheta = np.arccos(1-2*(np.cos(alpha)-1)/(np.cos(Theta)-1))
    F = -PPs.Rp**2 / PPs.Rs**2 /8 *(deltatheta + np.sin(deltatheta)) * (np.sin(phi0+ deltaphi/2) - np.sin(phi0- deltaphi/2))
    
    RESULT[:,1] = F
    return RESULT
   
def specular_diffuse_plot_theory(name_specular, name_diffuse, Obs_wave, transit = 'off'):
     # 绘制模拟的diffuse和specular的phase curve, 以及解析近似理论和光学理论结果，绘图中包含4曲线
     # 分别是： 模拟的specular和diffuse的phase curve, 解析近似理论和光学理论结果
    Is_specular, Ii_specular, Id_specular, It_specular, theta = data_loader(name_specular, Obs_wave)
    Is_diffuse, Ii_diffuse, Id_diffuse, It_diffuse, theta = data_loader(name_diffuse, Obs_wave)
    
    i = 0
    theta = theta/(2 *np.pi)

    fig, ax = plt.subplots(figsize=(9,6))
    if transit == 'off':
        ax.plot(theta, Is_specular[i,:] *1e6, label='specular', color='b', linewidth=2)
        ax.plot(theta, Id_diffuse[i,:] *1e6, label='diffuse', color='k', linewidth=2)
    else:
        # thermal同时包含了thermal emission 和transit的修正项， 前者为正值或0，后者为负值
        # 对于specualr_only and lambert_only模型, 这一步并不必要，因为本来就没有计算thermal emission，整个thermal 都是transit的修正项
        # 但如果不小心使用了非only的模型，那么这一步就是必要的， 需要将thermal emission 去除，只保留transit的修正项
        It_diffuse[It_diffuse > 0] =0   # 当然，直接置为零会带来一定的误差，但考虑到(thermal emission << transit)，这个误差是可以接受的
        It_specular[It_specular > 0] =0
        ax.plot(theta, (Is_specular[i,:] + It_specular[i,:]) *1e6, label='specular', color='b', linewidth=2)
        ax.plot(theta, (Id_diffuse[i,:] + It_diffuse[i,:]) *1e6, label='diffuse', color='k', linewidth=2)
        
    theory1 = analytical_theory_cal()
    # set all NAN to 0
    theory1 = np.nan_to_num(theory1, nan=0)
    # print(theory1)
    # theory1 = np.loadtxt('theory1.txt', delimiter = ',')
    ax.plot(theory1[:,0]/(2*np.pi), theory1[:,1] * 1e6, label='Our model', color='r', linewidth=2, linestyle='--')
    print('theory1:', theory1[(theory1.shape[0])//2, 1] *1e6)
    
    # 绘制一条平行于x轴的直线，颜色为'gray'，线宽为1
    theory2 = (PPs.Rp/2/PPs.semi_axis)**2 *1e6 # 21.1722  #21.3234
    ax.axhline(y=theory2, color='gray', linestyle='--', linewidth=2, label = 'optical')
    # ax.plot((0, theory2), (1, theory2), color='gray', linestyle='--', linewidth=2, label = 'virtual image')
    
    ax.set_xlabel('Orbital phase', fontsize=18)
    ax.set_ylabel(r'$F_p/F_*$ (ppm)', fontsize=18)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, np.max((Id_diffuse[i,:] + It_diffuse[i,:])) *1e6 *1.1)
    ax.spines['bottom'].set_linewidth(2)    ###设置底部坐标轴的粗细
    ax.spines['left'].set_linewidth(2)  ####设置左边坐标轴的粗细
    ax.spines['right'].set_linewidth(2) ###设置右边坐标轴的粗细
    ax.spines['top'].set_linewidth(2)   ####设置上部坐标轴的粗细
    #刻度值字体大小设置（x轴和y轴同时设置）
    plt.tick_params(labelsize=16)
    plt.legend(fontsize=17, frameon=False)
    plt.savefig(f"temp/{name_specular}/specular_diffuse_{Obs_wave[0]*1e6}_{transit}_theory.png")
    plt.savefig(f"temp/{name_diffuse}/specular_diffuse_{Obs_wave[0]*1e6}_{transit}_theory.png")
    plt.savefig(f"temp/{name_specular}/specular_diffuse_{Obs_wave[0]*1e6}_{transit}_theory.pdf")
    plt.savefig(f"temp/{name_diffuse}/specular_diffuse_{Obs_wave[0]*1e6}_{transit}_theory.pdf")
    plt.close()
      
    
def surface_model_compare(name_specular, name_diffuse, name_Fresnel, Obs_wave, transit = 'off'):
     # 绘制模拟的diffuse和specular的phase curve, 以及解析近似理论和光学理论结果，绘图中包含4曲线
     # 分别是： 模拟的specular和diffuse的phase curve, 解析近似理论和光学理论结果
    Is_specular, Ii_specular, Id_specular, It_specular, thetas = data_loader(name_specular, Obs_wave)
    Is_diffuse, Ii_diffuse, Id_diffuse, It_diffuse, thetad = data_loader(name_diffuse, Obs_wave)
    Is_Fresnel, Ii_Fresnel, Id_Fresnel, It_Fresnel, thetaf = data_loader(name_Fresnel, Obs_wave)
    
    i = 0
    thetas = thetas/(2 *np.pi)
    thetad = thetad/(2 *np.pi)
    thetaf = thetaf/(2 *np.pi)

    fig, ax = plt.subplots(figsize=(9,6))
    if transit == 'off':
        ax.plot(thetas, Is_specular[i,:] *1e6, label='specular', color='b', linewidth=2)
        ax.plot(thetad, Id_diffuse[i,:] *1e6, label='diffuse', color='k', linewidth=2)
        ax.plot(thetaf, Is_Fresnel[i,:] *1e6, label='Fresnel', color='g', linewidth=2)
    else:
        # thermal同时包含了thermal emission 和transit的修正项， 前者为正值或0，后者为负值
        # 对于specualr_only and lambert_only模型, 这一步并不必要，因为本来就没有计算thermal emission，整个thermal 都是transit的修正项
        # 但如果不小心使用了非only的模型，那么这一步就是必要的， 需要将thermal emission 去除，只保留transit的修正项
        It_diffuse[It_diffuse > 0] =0   # 当然，直接置为零会带来一定的误差，但考虑到(thermal emission << transit)，这个误差是可以接受的
        It_specular[It_specular > 0] =0
        It_Fresnel[It_Fresnel > 0] =0
        ax.plot(thetas, (Is_specular[i,:] + It_specular[i,:]) *1e6, label='specular', color='b', linewidth=2)
        ax.plot(thetad, (Id_diffuse[i,:] + It_diffuse[i,:]) *1e6, label='diffuse', color='k', linewidth=2)
        ax.plot(thetaf, (Is_Fresnel[i,:] + It_Fresnel[i,:]) *1e6, label='Fresnel', color='g', linewidth=2)
        
    theory1 = analytical_theory_cal()
    # set all NAN to 0
    theory1 = np.nan_to_num(theory1, nan=0)
    # print(theory1)
    # theory1 = np.loadtxt('theory1.txt', delimiter = ',')
    ax.plot(theory1[:,0]/(2*np.pi), theory1[:,1] * 1e6, label='Our model', color='r', linewidth=2, linestyle='--')
    print('theory1:', theory1[(theory1.shape[0])//2, 1] *1e6)
    
    # 绘制一条平行于x轴的直线，颜色为'gray'，线宽为1
    theory2 = (PPs.Rp/2/PPs.semi_axis)**2 *1e6 # 21.1722  #21.3234
    ax.axhline(y=theory2, color='gray', linestyle='--', linewidth=2, label = 'optical')
    # ax.plot((0, theory2), (1, theory2), color='gray', linestyle='--', linewidth=2, label = 'virtual image')
    
    ax.set_xlabel('Orbital phase', fontsize=18)
    ax.set_ylabel(r'$F_p/F_*$ (ppm)', fontsize=18)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, np.max((Id_diffuse[i,:] + It_diffuse[i,:])) *1e6 *1.1)
    ax.spines['bottom'].set_linewidth(2)    ###设置底部坐标轴的粗细
    ax.spines['left'].set_linewidth(2)  ####设置左边坐标轴的粗细
    ax.spines['right'].set_linewidth(2) ###设置右边坐标轴的粗细
    ax.spines['top'].set_linewidth(2)   ####设置上部坐标轴的粗细
    #刻度值字体大小设置（x轴和y轴同时设置）
    plt.tick_params(labelsize=16)
    plt.legend(fontsize=17, frameon=False)
    plt.savefig(f"temp/{name_specular}/specular_diffuse_{Obs_wave[0]*1e6}_{transit}_theory.png")
    plt.savefig(f"temp/{name_diffuse}/specular_diffuse_{Obs_wave[0]*1e6}_{transit}_theory.png")
    plt.savefig(f"temp/{name_specular}/specular_diffuse_{Obs_wave[0]*1e6}_{transit}_theory.pdf")
    plt.savefig(f"temp/{name_diffuse}/specular_diffuse_{Obs_wave[0]*1e6}_{transit}_theory.pdf")
    plt.close()
    
    
if __name__ == "__main__":
    # specular_diffuse_plot("R8copy", "R6copy", np.array([3]) * 1e-6, transit='off')
    # 在使用transit='on'时，注意'R1'和'R2'位置上的PC必须经过 transit_cal.py 的计算；应该为'R1copy'和'R2copy'的形式
    # specular_diffuse_plot_theory("specular_copy", "lambert_copy", np.array([3]) * 1e-6, transit='on')
    surface_model_compare("specular_copy", "lambert_copy", "Fresnel_copy", np.array([3]) * 1e-6, transit='on')
    
