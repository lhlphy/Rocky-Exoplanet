# 该脚本绘制模型对比图, 用于比较不同模型计算的结果 model_comparison.pdf
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from parameter_list import PPs

def data_loader(name, Obs_wave):
    '''
    load data, transform .file to .npy array
    name: midle name of temp folder
    Obs_wave: wavelength in um, used for interpolation (transform a matrix into vector)
    '''
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

    theta = np.linspace(0, 2*np.pi, 2000)
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
    '''
    绘制diffuse和specular的phase curve(*not* include thermal emission), 绘图中包含2曲线,
    分别是specular和diffuse的phase curve, 不区分表面材料性质的low与high albedo, 仅用于对比高分辨率的diffuse and specular_only model
    '''
    Is_specular, Ii_specular, Id_specular, It_specular, theta = data_loader(name_specular, Obs_wave)
    Is_diffuse, Ii_diffuse, Id_diffuse, It_diffuse, theta = data_loader(name_diffuse, Obs_wave)
    
    i = 0  # specular 与 diffuse与波长无关, 仅取决于几何结构, 所以i可以任意选取
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
 
def analytical_theory_cal(transit = 'on', An = 1, Fresnel = 'on', alpha = 0, Rp2Rs = 0):
    '''
    计算解析近似理论的phase 
    analytical_theory_cal中的A默认为1, 需要乘以由An决定的Fresnel系数
    transit: 是否考虑transit修正, 默认为'on'
    alpha: 可以自定义, 如果不指定, 则调用PPs计算
    Rp2Rs: Rp/Rs 可以自定义, 如果不指定, 则调用PPs计算
    
    attention: analytical model 仅alpha, Rp2Rs涉及到系统参数, 若不自行指定, 则调用PPs计算
    '''
    from parameter_list import PPs
    if (alpha == 0) and (Rp2Rs == 0):
        alpha = np.arcsin(PPs.Rs / (PPs.semi_axis - PPs.Rp))
        Rp2Rs = PPs.Rp / PPs.Rs
    
    Theta = np.linspace(0, 2*np.pi, 2000)
    RESULT = np.zeros((2000, 2))
    RESULT[:,0] = Theta
    
    deltaphi = alpha
    phi0 = (Theta + np.pi)/2
    theta0 = 0
    deltatheta = np.arccos(1-2*(np.cos(alpha)-1)/(np.cos(Theta)-1))
    print('alpha:', alpha)
    
    # F = -Rp2Rs**2 /8 *(deltatheta + np.sin(deltatheta)) * (np.sin(phi0+ deltaphi/2) - np.sin(phi0- deltaphi/2)) # original解析近似理论(未完全化简)
    F = Rp2Rs**2 *alpha**2 /4 *(1 - alpha**2 /24 *(2-np.cos(Theta)) /np.sin(Theta/2)**2) # 化简后的解析近似理论
    
    if Fresnel == 'on':  # 对于Fresnel模型, 需要乘以Fresnel系数
        F = F * PPs.A_Fresnel(I_angle= (np.pi-Theta)/2, A_normal= An)
    else:      # 若不考虑Fresnel效应, 直接乘以恒定的反射率常数 An
        F = F * An
    
    if transit == 'on': # transit修正, 将 pi-alpha<Theta<pi+alpha 的值置为0
        F[(np.pi - alpha <Theta) & (Theta < np.pi + alpha)] =0
        RESULT[:,1] = F
    elif transit == 'off':
        RESULT[:,1] = F 
        
    return RESULT 
   
def specular_diffuse_plot_theory(name_specular, name_diffuse, Obs_wave, transit = 'off', An = 1):
    '''
    绘制模拟的diffuse和specular的phase curve, 以及解析近似理论和光学理论结果, 绘图中包含4曲线
    分别是： 模拟的specular和diffuse的phase curve, 解析近似理论和光学理论结果
    An: 正入射的反射率, 本模型使用的Lambert模型默认反射率为1, 因此I_lambert需要乘以An; 
        specular模型若使用“Fresnel 0.2copy"这类, 其中已经包含了An=0.2的修正, 因此不需要乘以An;
        analytical_theory_cal中已经计算了Fresnel系数
    '''
    from parameter_list import PPs
    Is_specular, Ii_specular, Id_specular, It_specular, theta = data_loader(name_specular, Obs_wave)
    Is_diffuse, Ii_diffuse, Id_diffuse, It_diffuse, theta = data_loader(name_diffuse, Obs_wave)
    
    # smooth the phase curve 平滑修正
    alpha = np.arcsin(PPs.Rs / (PPs.semi_axis - PPs.Rp))
    dalpha = np.arcsin(PPs.Rp / (PPs.semi_axis - PPs.Rp))
    from scipy.signal import savgol_filter
    Is_specular[0, (theta < np.pi - alpha) & (theta > 2 *alpha)] = savgol_filter(Is_specular[0, (theta < np.pi - alpha) & (theta > 2 *alpha)], window_length=300, polyorder=3) # 左侧平滑化
    Is_specular[0, (theta > np.pi + alpha) & (theta < 2*np.pi -2 *alpha)] = savgol_filter(Is_specular[0, (theta > np.pi + alpha) & (theta < 2*np.pi -2 *alpha)], window_length=300, polyorder=3) # 右侧平滑化
    Is_specular[0, (theta > np.pi - alpha +dalpha) & (theta < np.pi + alpha-dalpha)] =0 # 中间transit 修正项置为0
    
    
    i = 0
    theta = theta/(2 *np.pi)

    fig, ax = plt.subplots(figsize=(9,6))
    if transit == 'off':
        ax.plot(theta, Id_diffuse[i,:]* 3/2 *1e6 *An, label='Lambert', color='k', linewidth=2)
        ax.plot(theta, Is_specular[i,:] *1e6, label='Numerical', color='b', linewidth=2)
    else:
        # thermal同时包含了thermal emission 和transit的修正项,  前者为正值或0, 后者为负值
        # 对于specualr_only and lambert_only模型, 这一步并不必要, 因为本来就没有计算thermal emission, 整个thermal 都是transit的修正项
        # 但如果不小心使用了非only的模型, 那么这一步就是必要的,  需要将thermal emission 去除, 只保留transit的修正项
        It_diffuse[It_diffuse > 0] = 0   # 当然, 直接置为零会带来一定的误差, 但考虑到(thermal emission << transit), 这个误差是可以接受的
        It_specular[It_specular > 0] = 0
        ax.plot(theta, (Id_diffuse[i,:] + It_diffuse[i,:])*3/2 *1e6 *An, label='Lambert', color='k', linewidth=2)
        ax.plot(theta, (Is_specular[i,:] + It_specular[i,:]) *1e6, label='Numerical', color='b', linewidth=2)
        # print((Is_specular[i,:] + It_specular[i,:]) *1e6)
        
    theory1 = analytical_theory_cal(An=0.2)
    # set all NAN to 0
    theory1 = np.nan_to_num(theory1, nan=0)
    # print(theory1)
    # theory1 = np.loadtxt('theory1.txt', delimiter = ',')
    ax.plot(theory1[:,0]/(2*np.pi), theory1[:,1] * 1e6, label='Analytical', color='r', linewidth=2, linestyle='--')
    print('theory1:', theory1[(theory1.shape[0])//2, 1] *1e6)
    
    # 绘制一条平行于x轴的直线, 颜色为'gray', 线宽为1
    theory2 = (PPs.Rp/2/PPs.semi_axis)**2 *1e6 # 21.1722  #21.3234
    ax.axhline(y=theory2 * An, color='gray', linestyle='--', linewidth=2, label = 'Optical')
    # ax.plot((0, theory2), (1, theory2), color='gray', linestyle='--', linewidth=2, label = 'virtual image')
    
    ax.set_xlabel('Orbital phase', fontsize=18)
    ax.set_ylabel(r'$F_p/F_*$ (ppm)', fontsize=18)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, np.max([np.max((Id_diffuse[i,:] + It_diffuse[i,:])*3/2 *An), np.max((Is_specular[i,:] + It_specular[i,:])), np.max(theory1[:,1])]) *1e6 *1.05)
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

    
def polarization_ploter(FR0 = 0.1, Obs_wave = np.array([3]) * 1e-6, transit = 'off', vertify = False, P_color = 'r'):
    '''
    plot the comparison of phase curve of P-pol and S-pol, and (non-polarization = P-pol + S-pol)
    
    vertify: 用另一组非极化计算结果验证, 正常False即可
    P_color: 不同FR的着色, 需要与Mathmetica Fresnel_R.pdf, surface_model_compare Fresnel{FR}.pdf 着色保持一致
             使用一组sequential colormap
    '''
    P_pol = f'P_pol_{FR0}_copy'
    S_pol = f'S_pol_{FR0}_copy'
    Is_P, Ii_P, Id_P, It_P, theta = data_loader(P_pol, Obs_wave)
    Is_S, Ii_S, Id_S, It_S, theta = data_loader(S_pol, Obs_wave)
    
    i = 0
    theta =theta/(2 *np.pi)
    # smooth the data and curve, extract the array data from matrix
    from scipy.signal import savgol_filter
    Is_P = Is_P[i,:]
    Is_S = Is_S[i,:]
    It_S = It_S[i,:]
    It_P = It_P[i,:]
    
    fig, ax = plt.subplots(figsize=(9,6))
    if transit == 'off':
        ax.plot(theta, Is_P *1e6, label='P-polarization', color=P_color, linewidth=2.5, linestyle='solid')
        ax.plot(theta, Is_S *1e6, label='S-polarization', color=P_color, linewidth=2.5, linestyle='dashed')
        ax.plot(theta, (Is_S + Is_P) *1e6, label='Net-polarization', color='k', linewidth=2.5)
    else:
        It_P[It_P > 0] =0   # 当然, 直接置为零会带来一定的误差, 但考虑到(thermal emission << transit), 这个误差是可以接受的
        It_S[It_S > 0] =0
        ax.plot(theta, (Is_P + It_P) *1e6, label='P-polarization', color=P_color, linewidth=2.5, linestyle='solid')
        ax.plot(theta, (Is_S + It_S) *1e6, label='S-polarization', color=P_color, linewidth=2.5, linestyle='dashed')
        ax.plot(theta, (Is_S + Is_P + It_P) *1e6, label='Net-polarization', color='k', linewidth=2.5)
        if (It_P != It_S).all():
            raise ValueError('It_P and It_S are not equal')
        
    # 绘制非极化计算结果 --polarization None, 用于验证P,S极化计算结果
    if vertify:  
        no_pol = f'No_pol_{FR0}_copy'
        Is_no, Ii_no, Id_no, It_no, _ = data_loader(no_pol, Obs_wave)
        
        if transit == 'off':
            ax.plot(theta, Is_no *1e6, label='Net-polarization', color='b', linewidth=2.5, linestyle='dashed')
        else:
            It_no[It_no > 0] =0   # 当然, 直接置为零会带来一定的误差, 但考虑到(thermal emission << transit), 这个误差是可以接受的
            ax.plot(theta, (Is_no + It_no) *1e6, label='Net-polarization', color='b', linewidth=2.5, linestyle='dashed')
            
    ax.set_xlabel('Orbital phase', fontsize=20)
    ax.set_ylabel(r'$F_p/F_*$ (ppm)', fontsize=20)
    ax.set_xlim(0, 1)
    # ax.set_ylim(-440, -430)
    if transit == 'off':
        ax.set_ylim(0, np.max((Is_S+Is_P)*1e6 *1.05))
    else:
        ax.set_ylim(0, np.max((Is_S+Is_P+It_P)*1e6 *1.05))
    ax.spines['bottom'].set_linewidth(2)    ###设置底部坐标轴的粗细
    ax.spines['left'].set_linewidth(2)  ####设置左边坐标轴的粗细
    ax.spines['right'].set_linewidth(2) ###设置右边坐标轴的粗细
    ax.spines['top'].set_linewidth(2)   ####设置上部坐标轴的粗细
    #刻度值字体大小设置（x轴和y轴同时设置）
    plt.tick_params(labelsize=18)
    if FR == 0.1:  # 仅在FR=0.1时绘制legend, 4个子图只需要一个legend
        plt.legend(fontsize=20, frameon=False)
    
    # set different info for .pdf name
    if transit == 'on':
        info = 'Transit'
    elif transit == 'off':
        info = 'noTransit'
    # plt.savefig(f"temp/{P_pol}/Polarization_{Obs_wave[0]*1e6}_{FR0}_{info}.png")
    # plt.savefig(f"temp/{S_pol}/Polarization_{Obs_wave[0]*1e6}_{FR0}_{info}.png")
    plt.savefig(f"temp/{P_pol}/Polarization_{Obs_wave[0]*1e6}_{FR0}_{info}.pdf")
    plt.savefig(f"temp/{S_pol}/Polarization_{Obs_wave[0]*1e6}_{FR0}_{info}.pdf")
    plt.close()
    
    ### 绘制偏振度PC
    fig, ax = plt.subplots(figsize=(9,6))
    Degree_polarization = np.abs((Is_P - Is_S) /(Is_P + Is_S))
    ax.plot(theta, Degree_polarization, linewidth=2.5)
    ax.set_xlabel('Orbital phase', fontsize=18)
    ax.set_ylabel('Degree of polarization', fontsize=18)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.spines['bottom'].set_linewidth(2)    ###设置底部坐标轴的粗细
    ax.spines['left'].set_linewidth(2)  ####设置左边坐标轴的粗细
    ax.spines['right'].set_linewidth(2) ###设置右边坐标轴的粗细
    ax.spines['top'].set_linewidth(2)   ####设置上部坐标轴的粗细
    #刻度值字体大小设置（x轴和y轴同时设置）
    plt.tick_params(labelsize=16)
    
    plt.savefig(f"temp/{P_pol}/Degree_polarization_{Obs_wave[0]*1e6}_{FR0}.pdf")
    plt.savefig(f"temp/{S_pol}/Degree_polarization_{Obs_wave[0]*1e6}_{FR0}.pdf")
    plt.close()
    

def surface_model_compare(name_specular, name_diffuse, name_Fresnel, Obs_wave, transit = 'off', FR = PPs.std_FR, F_color = 'r'):
    '''
    plot the comparison of phase curve of specular, diffuse and Fresnel, 绘图中包含3曲线
    emphasize the difference of Fresnel model from specular
    Parameters:
    FR: Lambert and Specular input model have R=1, so when compare with Fresnel, multiply FR
    F_color: 不同FR的着色, 需要与 Mathmetica Fresnel_R.pdf, polarization_ploter 着色保持一致
    '''
    Is_specular, Ii_specular, Id_specular, It_specular, thetas = data_loader(name_specular, Obs_wave)
    Is_diffuse, Ii_diffuse, Id_diffuse, It_diffuse, thetad = data_loader(name_diffuse, Obs_wave)
    Is_Fresnel, Ii_Fresnel, Id_Fresnel, It_Fresnel, thetaf = data_loader(name_Fresnel, Obs_wave)
    
    i = 0
    thetas = thetas/(2 *np.pi)
    thetad = thetad/(2 *np.pi)
    thetaf = thetaf/(2 *np.pi)

    fig, ax = plt.subplots(figsize=(9,6))
    if transit == 'off':
        ax.plot(thetad, Id_diffuse[i,:] *1e6, label='Lambert', color='k', linewidth=2.5)
        ax.plot(thetas, Is_specular[i,:] *1e6, label='Specular', color='b', linewidth=2.5)
        ax.plot(thetaf, Is_Fresnel[i,:] *1e6, label='Fresnel', color=F_color, linewidth=2.5)
    else:
        # thermal同时包含了thermal emission 和transit的修正项,  前者为正值或0, 后者为负值
        # 对于specualr_only and lambert_only模型, 这一步并不必要, 因为本来就没有计算thermal emission, 整个thermal 都是transit的修正项
        # 但如果不小心使用了非only的模型, 那么这一步就是必要的,  需要将thermal emission 去除, 只保留transit的修正项
        It_diffuse[It_diffuse > 0] =0   # 当然, 直接置为零会带来一定的误差, 但考虑到(thermal emission << transit), 这个误差是可以接受的
        It_specular[It_specular > 0] =0
        It_Fresnel[It_Fresnel > 0] =0
        ax.plot(thetad, (Id_diffuse[i,:] + It_diffuse[i,:])*FR *3/2 *1e6, label='Lambert', color='k', linewidth=2.5)
        ax.plot(thetas, (Is_specular[i,:] + It_specular[i,:])*FR *1e6, label='Specular', color='b', linewidth=2.5)
        ax.plot(thetaf, (Is_Fresnel[i,:] + It_Fresnel[i,:]) *1e6, label='Fresnel', color=F_color, linewidth=2.5)
        
    # ### 绘制解析近似理论结果 Analytical result
    # theory1 = analytical_theory_cal() # 解析近似理论结果
    # # set all NAN to 0
    # theory1 = np.nan_to_num(theory1, nan=0)
    # # print(theory1)
    # # theory1 = np.loadtxt('theory1.txt', delimiter = ',')
    # ax.plot(theory1[:,0]/(2*np.pi), theory1[:,1]* FR * 1e6, label='Analytical', color='r', linewidth=2, linestyle='--')
    # print('theory1:', theory1[(theory1.shape[0])//2, 1] *1e6)
    
    # ### 绘制光学理论结果 Simple optical result
    # theory2 = (PPs.Rp/2/PPs.semi_axis)**2 *1e6 * FR # 21.1722  #21.3234  # calculate the simple optical result
    # ax.axhline(y=theory2, color='gray', linestyle='--', linewidth=2, label = 'Optical')
    # # ax.plot((0, theory2), (1, theory2), color='gray', linestyle='--', linewidth=2, label = 'virtual image')
    
    ax.set_xlabel('Orbital phase', fontsize=20)
    ax.set_ylabel(r'$F_p/F_*$ (ppm)', fontsize=20)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, max (np.max(Id_diffuse[i,:] + It_diffuse[i,:]) *3/2 * FR, np.max(Is_Fresnel[i,:] + It_Fresnel[i,:])) *1e6 *1.05)  # *1.4 for Fresnel0.1
    ax.spines['bottom'].set_linewidth(2)    ###设置底部坐标轴的粗细
    ax.spines['left'].set_linewidth(2)  ####设置左边坐标轴的粗细
    ax.spines['right'].set_linewidth(2) ###设置右边坐标轴的粗细
    ax.spines['top'].set_linewidth(2)   ####设置上部坐标轴的粗细
    #刻度值字体大小设置（x轴和y轴同时设置）
    plt.tick_params(labelsize=18)
    if FR == 0.1:  # 仅在FR=0.1时绘制legend, 4个子图只需要一个legend
        plt.legend(fontsize=20, frameon=False)
    plt.savefig(f"temp/{name_Fresnel}/Fresnel{FR}.png")
    plt.savefig(f"temp/{name_Fresnel}/Fresnel{FR}.pdf")
    plt.close()
    
if __name__ == "__main__":
    # specular_diffuse_plot("R8copy", "R6copy", np.array([3]) * 1e-6, transit='off')
    # 在使用transit='on'时, 注意'R1'和'R2'位置上的PC必须经过 transit_cal.py 的计算；应该为'R1copy'和'R2copy'的形式
    # specular_diffuse_plot_theory("specular_copy", "lambert_copy", np.array([3]) * 1e-6, transit='on')
    specular_diffuse_plot_theory("Trappist1b_Fresnel 0.2copy", "Trappist1b_Fresnel 0.2copy", np.array([3]) * 1e-6, transit='on', An = 0.2)
    # specular_diffuse_plot_theory("R1copy", "R1copy", np.array([3]) * 1e-6, transit='on', FR0= 0.1)
    
    # ### 从plasma色图中均匀取出4个颜色
    # cmap = plt.get_cmap('plasma')
    # color_list = [cmap(i) for i in [0.5, 0.4, 0.2, 0]]
    # print(color_list)
    # ### 绘制不同FR的Fresnel模型与理论结果的对比 Appendix:Fresnel
    # FR_list = [0.1, 0.2, 0.4, 0.8]
    # # # color_list = [(247/255, 193/255, 198/255), (240/255, 141/255, 149/255), (232/255, 71/255, 85/255), (199/255, 25/255, 40/255)]
    # for FR, color in zip(FR_list, color_list):
    #     surface_model_compare("specular_copy", "lambert_copy", f"Fresnel {FR}copy", np.array([3]) * 1e-6, transit='on', FR=FR, F_color=color)
    
    # ## 绘制不同FRnormal下的 P,S偏振光以及非偏振光的phase curve, 并绘制偏振度PC Appendix:Pol
    # FR_list = [0.1, 0.2, 0.4, 0.8]
    # # color_list = [(247/255, 193/255, 198/255), (240/255, 141/255, 149/255), (232/255, 71/255, 85/255), (199/255, 25/255, 40/255)]
    # for FR, color in zip(FR_list, color_list):
    #     polarization_ploter(FR0 =FR, vertify = False, transit= 'on', P_color=color)
    #     polarization_ploter(FR0 =FR, vertify = False, transit= 'off', P_color=color)


