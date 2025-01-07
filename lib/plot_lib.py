import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os
os.environ['roughness'] = '0'
import argparse
from bond_albedo_calculator import bond_albedo_calculator
from function_library import chi2_cal
from save_txt import save_txt
from parameter_list import PPs

def IDS_plot(name, Obs_wave):
    """
    plot the phase curve of total intensity, diffusion, specular reflection in one figure
    name: the middle name of the file
    Obs_wave: the wavelength that should be plotted
    """
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
    
    i = 0
    N  = I_s.size
    print('Wavelength: ', Obs_wave[0]*1e6, 'um')
    print('Secondary eclipse depth difference: ',(I_d[i,N//2]-I_s[i,N//2])*1e6,' ppm')
    
    # plot
    plt.rcParams['font.size'] = 12
    
    # plt.figure(figsize=(8, 5))

    # plt.plot(theta, I_s[0,:] * 1e6, label='Specular')
    # plt.plot(theta, I_d[0,:] * 1e6, label='Diffuse')
    # plt.plot(theta, I_t[0,:] * 1e6, label='Thermal')

    # plt.xlabel('Orbital phase')
    # plt.ylabel('Contrast ratio (ppm)')
    # plt.legend()
    # plt.show()
    # os.makedirs(f'temp/P0', exist_ok= True)
    # plt.savefig(f'temp/P0/compare.png')
    # plt.close()
    theta = theta/(2 *np.pi)

    fig, ax = plt.subplots(figsize=(9,6))
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
    ax.plot(theta, I_s[i,:] * 1e6, label='Specular', color='k', linewidth=2)
    ax.set_ylim(ymin=np.min(I_s[i,:]), ymax= np.max(I_s[i,:])*1.1e6)
    ax.set_xlabel('Orbital phase', fontsize=18)
    ax.set_ylabel(r'$F_p/F_*$ (ppm)', fontsize=18)
    ax.tick_params(length=6, width=2)
    ax.spines['right'].set_visible(False)
    np.save(f'temp/{name}/variables/sim_theta.npy', theta)
    np.save(f'temp/{name}/variables/sim_SR.npy', I_s[i,:])

    lambda_color = 'blue'
    labmda_ax = ax.twinx()
    labmda_ax.set_position(axpos)
    labmda_ax.plot(theta, I_d[i,:] * 1e6, label='Diffuse', color=lambda_color, linewidth=2)
    labmda_ax.set_ylim(ymin=np.min(I_s[i,:]), ymax= np.max(I_d[i,:])*1.1e6)
    labmda_ax.set_xlabel('Orbital phase', fontsize=18)
    labmda_ax.tick_params(length=6, width=2, color=lambda_color, labelcolor=lambda_color)
    # labmda_ax.set_ylabel('Contrast ratio (ppm)', fontsize=18, color=lambda_color)
    labmda_ax.spines['right'].set(color=lambda_color, linewidth=2.0, linestyle=':')

    omglog_color = 'red'
    omglog_ax = ax.twinx()
    # 使用科学计数法的刻度
    omglog_ax.ticklabel_format(style='plain', axis='y', scilimits=(0,0))
    # 获取 y 轴 OffsetText 对象
    offset_text = omglog_ax.yaxis.get_offset_text()
    # 调整位置示例，偏移 (dx, dy) 单位是像素 (points)
    offset_text.set_position((1.12, 0))
    # 调整字体大小
    offset_text.set_size(18)  # 或者使用 offset_text.set_fontsize(12)
    omglog_ax.spines['right'].set_position(('data', np.max(theta)*1.1))
    omglog_ax.set_ylim(np.min(I_s[i,:]), np.max(I_t[i,:])*1.1e6)
    omglog_ax.set_position(axpos)
    omglog_ax.plot(theta, I_t[i,:] * 1e6, label='Thermal', color=omglog_color, linewidth=2)
    # omglog_ax.set_ylabel('Contrast ratio (ppm)', fontsize=18, color=omglog_color)
    omglog_ax.tick_params(length=6, width=2, color=omglog_color, labelcolor=omglog_color)
    omglog_ax.spines['right'].set(color=omglog_color, linewidth=2.0, linestyle='-.')
    fig.subplots_adjust(bottom=0.2) 
    fig.legend(['Specular', 'Diffuse', 'Thermal'], ncol = 3, loc = 'lower center',prop={'weight': 'bold', 'size': 18})
    # fig.legend(prop={'weight': 'bold'}) 
    # os.makedirs(f'temp/P3', exist_ok= True)
    plt.savefig(f'temp/{name}/compare_{Obs_wave[0]*1e6}.png')
    plt.close()

    fig, ax = plt.subplots(figsize=(9,6))
    ax.plot(theta, (I_s[i,:]+I_t[i,:])*1e6 , label = 'Specular')
    ax.plot(theta, (I_d[i,:]+I_t[i,:])*1e6 , label = 'Diffuse')
    ax.set_xlabel('Orbital phase', fontsize=18)
    ax.set_ylabel('F_p/F_* (ppm)', fontsize=18)

    plt.errorbar(theta[N//2], (I_s[i,N//2]+I_t[i,N//2])*1e6, yerr = 1.25, capsize= 10)
    plt.plot(theta[N//2], (I_s[i,N//2]+I_t[i,N//2])*1e6,'.')

    plt.legend(['Specular', 'Diffuse'])
    plt.savefig(f'temp/{name}/compare_PC_{Obs_wave[0]*1e6}.png')
    plt.close()


def spectrum_plot(name, wave_range):
    # load data, I_diffuse,I_specular is contrast ratio 
    I_diffuse = np.load(f'temp/{name}/variables/I_diffuse.npy')
    I_specular = np.load(f'temp/{name}/variables/I_specular.npy')
    Star_flux = np.load(f'temp/{name}/variables/Star_flux.npy')
    wave_list = np.load(f'temp/{name}/variables/wave_list.npy')
    Thermal = np.load(f'temp/{name}/variables/Thermal.npy')

    # convert I_diffuse,I_specular (contrast ratio) to absolute intensity
    N1 = Star_flux.shape[0]
    N2 = Star_flux.shape[1]
    ID = I_diffuse[:, N2//2]
    IS = I_specular[:, N2//2]
    Star_flux = Star_flux[:, N2//2]
    IT = Thermal[:, N2//2]
    
    # set to the wavelength range
    indx = (wave_list > wave_range[0]) & (wave_list < wave_range[1])
    wave_list = wave_list[indx]
    IT = IT[indx]
    ID = ID[indx]
    IS = IS[indx]
    
    plt.plot(wave_list *1e6, (IT + ID) *1e6)
    plt.plot(wave_list *1e6, (IT + IS) *1e6)
    plt.legend(['Lambert surface', 'Specular surface'])
    plt.xlabel('Wavelength (μm)')
    plt.ylabel('Contrast ratio (ppm)')
    
    plt.show()
    plt.savefig(f"temp/{name}/spectrum")
    
def compare_phase_curve_plot(name_list, wave_range, instrument = '  ', legend = 'below', xlabel = 'on', ylabel = 'on', errorbar = 0):
    '''
    绘制主要full phase curve, 默认情况下name_list只有两个name, 分别代表low albedo和high albedo, 绘制4条曲线
    分别是: low albedo & Lambert, low albedo & Specular, high albedo & Lambert, high albedo & Specular
    可绘制 OpticalFrame == 'Full_cal' or 'Non_Fresnel' 两种情况下的phase_curve_comp plot
    
    legend: 'below', 'insert', 'off'
        'below' means legend below center the plot
        'insert' means legend in the plot
        'off' means no legend
        
    ylabel, xlabel: 'on', 'off': on or off the y and x label, because when four pics merge together, only the boundary axis label is needed
    
    errorbar: if errorbar == 0, no errorbar; Otherwise, draw a errorbar
    '''
    # load data, I_diffuse,I_specular is contrast ratio 
    sim_obs = np.loadtxt(f"telescope_measure/sim_obs (14).txt", delimiter=' ')
    x = sim_obs[:,0]
    y = sim_obs[:,3] *1e6
    Bin = np.sqrt(1/np.sum(1/y**2))
    
    pallet = ['b','r','k']
    plotarr  = [0] * 8
    fig, ax = plt.subplots()
    
    ebar = np.array([Bin, Bin * np.sqrt(1/0.629), Bin])
    tbar = np.array([0.0647, 0.0407, 0.0647])
    
    xloc = np.array([0.3948, 0.5, 0.6052])
    # Period: 7.72614 h
    Theta_list = np.load(f'temp/{name_list[0]}/variables/Theta.npy')
    Theta_list = Theta_list / (2 *np.pi)   # 将相位角归一化
    data = np.zeros([np.size(Theta_list), 6])
    
    up_bound = 0  # control xlim up_lim
    for i, name in enumerate(name_list):
        Theta_list = np.load(f'temp/{name}/variables/Theta.npy')
        Theta_list = Theta_list / (2 *np.pi)   # 将相位角归一化
        Nt = np.size(Theta_list)
        CR_S = np.zeros([Nt])
        CR_D = np.zeros([Nt])
        for j, theta in enumerate(Theta_list):
            CR_S[j], CR_D[j] = bond_albedo_calculator(wave_range[0], wave_range[1], name, j)
        
        plotarr[i*2], = plt.plot(Theta_list, CR_D, '-', color = pallet[i], linewidth = 2)
        plotarr[i*2+1], = plt.plot(Theta_list, CR_S, '--', color = pallet[i], linewidth = 2)
        data[:,0] = Theta_list
        up_bound = np.max([np.max(CR_D), np.max(CR_S), up_bound]) # find the max value for ylim
        if i != 2:
            data[:,2 + i*2] = CR_S
            data[:,3 + i*2] = CR_D
        elif i == 2:
            data[:,1] = CR_S
        
        spl = interp1d(Theta_list, CR_D, kind='linear')
        yloc = spl(xloc)
        for k, xl in enumerate(xloc):
            if i != 2:
                print(' ')
                # plt.errorbar(xloc[k], yloc[k], yerr = ebar[k], xerr=tbar[k], fmt='o', color = pallet[i], ecolor=pallet[i], linestyle='None')
            # plt.errorbar(xloc[k], yloc[k] + ebar[k] * np.random.random(), yerr = ebar[k], xerr=tbar[k], fmt='o', color = pallet[i], ecolor=pallet[i], linestyle='None')
        
        # 调整布局以便为图例腾出空间  
    fig.subplots_adjust(bottom=0.25) 
    save_txt(data, "Orbital Phase(normalized)\tFp/F*(Blackbody ppm)\tFp/F*(Low albedo & Specular ppm)\tFp/F*(Low albedo & Lambert ppm)\tFp/F*(High albedo & Specular ppm)\tFp/F*(High albedo & Lambert ppm)", f"phase_curve.txt") 

    print(ebar)
    # 设置背景颜色  
    # ax.axvspan(0.4593, 0.5407, color='gray', alpha=0.5)
    # ax.axvspan(0.3303, 0.4593, color='gray', alpha=0.2)
    # ax.axvspan(0.5407, 0.6697, color='gray', alpha=0.2)
    # 设置图例并放置在图窗的正下方  
    plt.ylim([0,up_bound * 1.1])
    plt.xlim([0,1])
    
    # 当 errorbar != 0 时，在坐标系左上角添加一个带误差棒的点
    if errorbar != 0:
        if legend == 'insert':
            ax.errorbar(0.5, up_bound *0.9, yerr=errorbar, fmt='o', color='black', markersize=4)
        else:
            ax.errorbar(0.1, up_bound *0.9, yerr=errorbar, fmt='o', color='black', markersize=4)
    
    # 设置图例位置
    if legend == 'below':  
        plt.legend([plotarr[0],plotarr[2],plotarr[4],plotarr[1],plotarr[3]], ['Low albedo & Lambert', 'High albedo & Lambert', 'Blackbody', 'Low albedo & Specular', 'High albedo & Specular'], loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2)
    elif legend == 'insert':
        plt.legend([plotarr[0],plotarr[2],plotarr[4],plotarr[1],plotarr[3]], [r'Low $A$ & Lambert', r'High $A$ & Lambert', 'Blackbody', r'Low $A$ & Specular', r'High $A$ & Specular'], loc='upper left', bbox_to_anchor=(0, 1.01), fontsize=10, frameon=False)
    elif legend == 'off':
        pass
    else:
        raise ValueError("Invalid legend position specified.")

    # plt.legend(plotarr, ['Low albedo & Lambert', 'Low albedo & Specular',  'Mid albedo & Lambert', 'Mid albedo & Specular',
    #             'High albedo & Lambert', 'High albedo & Specular'], loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2)
    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2)  
    if xlabel == 'on':
        plt.xlabel('Orbital Phase', fontsize = 13)
    if ylabel == 'on':
        plt.ylabel(r'$F_p/F_*$ (ppm)', fontsize = 13)
    # 添加文字和箭头，设置字体透明度  
    # plt.text(0.3948, 0, '1 hour', fontsize=10,fontweight='bold', color='black', ha='center') 
    # # plt.text(0.3948, -5, '', fontsize=12, color='black', ha='center')
    # plt.text(0.5, -4, '0.63 h', fontsize=10,fontweight='bold', color='black', ha='center')
    # # plt.text(0.3948, -3, '1 hour', fontsize=12, color='black', ha='center')
    # plt.text(0.6052, 0, '1 hour', fontsize=10,fontweight='bold', color='black', ha='center')
    # # plt.text(0.3948, -3, '1 hour', fontsize=12, color='black', ha='center')
    
    # plt.text(0.5,-100, 'NIRCam/F444W', fontsize = 10, ha='center', fontweight='bold')
    # plt.text(0.5,-120, '3.9-5 μm', fontsize = 10, ha='center', fontweight='bold')
    plt.text(0.81, up_bound *1 , instrument, fontsize = 10, ha='center', fontweight='bold') # label instrument name and wavelength range
    plt.text(0.81, up_bound *0.92, f'{wave_range[0]*1e6 :.2f}-{wave_range[1]*1e6 :.2f} μm', fontsize = 10, ha='center', fontweight='bold')
    # plt.axis([0, 1, -5, 40])
    Transit_depth = PPs.Rp**2 / PPs.Rs**2 * 1e6
    print(f'Transit depth: {Transit_depth:.2f} ppm')

    # plt.savefig(f"temp/{name}/phase_curve_comp1.pdf", format = 'pdf')
    plt.savefig(f"temp/{name_list[0]}/phase_curve_comp_{instrument.replace('/','_')}.pdf", format = 'pdf')
    plt.savefig(f"temp/{name_list[1]}/phase_curve_comp_{instrument.replace('/','_')}.pdf", format = 'pdf')
    plt.show()
    plt.close()
    
    
   
def real_comp(name_list):
    # load data, I_diffuse,I_specular is contrast ratio 
    ############ load real measured data  ###############
    data_value = np.loadtxt(f"telescope_measure/JWST_data.txt", delimiter=',')
    data_ebar = np.loadtxt(f"telescope_measure/JWST_errorbar.txt", delimiter=',')
    
    DN = data_value.shape[0]
    plotarr  = [0] * (2 * np.size(name_list) +2)
    
    data = np.zeros([DN, 3])
    data[:,0:2] = data_value
    for i in range(DN):
        data[i, 2] = np.abs(data_ebar[i*2,1] - data_ebar[i*2+1,1])/2
        
    # 创建散点图
    save_txt(data)  
    fig, ax = plt.subplots()
    plt.errorbar(data[:,0], data[:,1], yerr=data[:,2], fmt='ok', ecolor='k', linestyle='None', label = 'JWST')  
    
    # 处理离群点，计算chi2是否要去除离群点，store the new data after removing outliers in data2
    DN = DN # the size of data2
    Data2 = np.zeros([DN, 3])
    Data2[0:7] = data[0:7]  # the data point before outlier
    Data2[7:] = data[7:]  # the data point after outlier
    
    
    ## 处理模型数据并绘图
    pallet = ['#033eff', '#ff033e',  'black', '#3eff03', '#03bcff', '#c403ff'] # color library
    chi2 = np.zeros(np.size(name_list)*2 + 2)  # 
    
    # processing error bar
    sim_obs = np.loadtxt(f"telescope_measure/sim_obs (13).txt", delimiter=' ')
    x = sim_obs[:,0]
    y = sim_obs[:,3] *1e6
    Nbin = 6 # the number of bin: 6
    Bin = np.zeros(Nbin)  
    xloc = np.zeros(Nbin)
    NGroup = sim_obs.shape[0] // Nbin # the number of pixels in each bin
    for i in range(Nbin):  # calculate the errorbar of each bins
        xq = x[NGroup*i:NGroup*(i+1)]  # group pixel points to each bin
        yq = y[NGroup*i:NGroup*(i+1)]
        Bin[i] = np.sqrt(1/np.sum(1/yq**2)) # calculate the errorbar of each bins
        xloc[i] = np.mean(xq) # calculate the mean wavelength of each bin
        
    print([f'Pixels: {sim_obs.shape[0]}, bins number: {Nbin}, pixels number in each bin: {NGroup}'])
    
    
    for i, name in enumerate(name_list):  # Plot the spectra of each dataset
        I_diffuse = np.load(f'temp/{name}/variables/I_diffuse.npy')
        I_specular = np.load(f'temp/{name}/variables/I_specular.npy')
        Star_flux = np.load(f'temp/{name}/variables/Star_flux.npy')
        wave_list = np.load(f'temp/{name}/variables/wave_list.npy')
        Thermal = np.load(f'temp/{name}/variables/Thermal.npy')

        # convert I_diffuse,I_specular (contrast ratio) to absolute intensity
        N1 = Star_flux.shape[0]
        N2 = Star_flux.shape[1]
        ID = I_diffuse[:, N2//2]
        IS = I_specular[:, N2//2]
        Star_flux = Star_flux[:, N2//2]
        IT = Thermal[:, N2//2]
        
        if i == 0:
            data1 = np.zeros([np.size(wave_list), 3])
            data1[:,0] = wave_list *1e6
            data1[:,1] = (IT + IS) *1e6
            data1[:,2] = (IT + ID) *1e6
            # save_txt(data, "Wavelength(micron)\tFp/F*(Low albedo_specular ppm)\tFp/F*(Low albedo_Lambert ppm)", "Low albedo.txt")
        elif i == 1:
            data2 = np.zeros([np.size(wave_list), 3])
            data2[:,0] = wave_list *1e6
            data2[:,1] = (IT + IS) *1e6
            data2[:,2] = (IT + ID) *1e6
            # save_txt(data, "Wavelength(micron)\tFp/F*(High albedo_specular ppm)\tFp/F*(High albedo_Lambert ppm)", "High albedo.txt")
        elif i == 2:
            data3 = np.zeros([np.size(wave_list), 3])
            data3[:,0] = wave_list *1e6
            data3[:,1] = (IT + IS) *1e6
            data3[:,2] = (IT + ID) *1e6
            # save_txt(data, "Wavelength(micron)\tFp/F*(Blackbody ppm)\tFp/F*(Blackbody ppm)", "Blackbody.txt")
        # elif i == 3:
        #     data4 = np.zeros([np.size(wave_list), 3])
        #     data4[:,0] = wave_list *1e6
        #     data4[:,1] = (IT + IS) *1e6
        #     data4[:,2] = (IT + ID) *1e6
            # save_txt(data, "Wavelength(micron)\tFp/F*(Blackbody & heat fully redistribution ppm)\tFp/F*(Blackbody & heat fully redistribution ppm)", "Heat.txt")
            
    
        if i != 2:
            plotarr[i], = plt.plot(wave_list *1e6, (IT + ID) *1e6, '--', color = pallet[i]) # lambert diffuse case
            plotarr[i+3], = plt.plot(wave_list *1e6, (IT + IS) *1e6, '-', color = pallet[i]) # specular case
        else:  # blackbody
            plotarr[i], = plt.plot(wave_list *1e6, (IT + ID) *1e6, '-', color = pallet[i]) # lambert diffuse case
        
        if i == 0:
            spl = interp1d(wave_list *1e6,  (IT + ID) *1e6, kind='linear')
            yloc = spl(xloc)
            # plt.errorbar(xloc, yloc, yerr=Bin, fmt='o', ecolor=pallet[i], linestyle='None')
            plt.errorbar(xloc, yloc + Bin* (np.random.random(np.size(Bin)) - 0.5) *2, yerr=Bin, fmt='o', ecolor=pallet[i], linestyle='None')
        
        chi2[i] = chi2_cal(Data2[:,0], Data2[:,1], Data2[:,2], wave_list *1e6, (IT + ID) *1e6)
        chi2[i+3] = chi2_cal(Data2[:,0], Data2[:,1], Data2[:,2], wave_list *1e6, (IT + IS) *1e6)
        
        
    data = np.zeros([np.size(wave_list), 7])
    data[:,0] = data3[:,0]
    data[:,1] = data3[:,1]
    data[:,2] = data1[:,1]
    data[:,3] = data1[:,2]
    data[:,4] = data2[:,1]
    data[:,5] = data2[:,2]
    # data[:,6] = data4[:,1]
    name = "Wavelength(micron)\tFp/F*(Blackbody ppm)\tFp/F*(Low albedo_specular ppm)\tFp/F*(Low albedo_Lambert ppm)\tFp/F*(High albedo_specular ppm)\tFp/F*(High albedo_Lambert ppm)\tFp/F*(Blackbody & heat fully redistribution ppm)"
    save_txt(data, name, "Model_result.txt")
    
    # draw Al2O3 surface
    GJAl2O3 = np.loadtxt('GJAl2O3.txt', skiprows = 1)
    # print(GJAl2O3)
    wave_GJ = GJAl2O3[:,0]
    Al2O3_GJ = GJAl2O3[:,1]
    Black_GJ = GJAl2O3[:,2]
    
    plotarr[5], = plt.plot(wave_GJ, Al2O3_GJ, label = 'Al_2O_3 surface', color = 'orange')
    # plotarr[-1], = plt.plot(wave_GJ, Black_GJ, label = 'Blackbody (XT)')
    chi2[5] = chi2_cal(Data2[:,0], Data2[:,1], Data2[:,2], wave_GJ, Al2O3_GJ)
    # chi2[-1] = chi2_cal(Data2[:,0], Data2[:,1], Data2[:,2], wave_GJ, Black_GJ)
    
        
    fig.subplots_adjust(bottom=0.4)  
    # plt.legend(plotarr, [f'Low albedo & Lambert $\chi^2$={chi2[0]:.2f}', f'Low albedo & Specular $\chi^2$={chi2[1]:.2f}', 
    #             f'High albedo & Lambert $\chi^2$={chi2[2]:.2f}', f'High albedo & Specular $\chi^2$={chi2[3]:.2f}', 
    #             'Blackbody & Lambert', 'Blackbody & Specular', 'Full redistribution & Lambert', 'Full redistribution & Specular',
    #             'O-H & Lambert', 'O-H & Specular'], loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2)
    # plt.legend(plotarr, ['Low albedo & Lambert', 'Low albedo & Specular',  'Mid albedo & Lambert', 'Mid albedo & Specular',
    #             'High albedo & Lambert', 'High albedo & Specular'], loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2)
    plt.legend(plotarr, [ f'Low albedo & Lambert $\chi^2$={chi2[0]:.2f}', f'High albedo & Lambert $\chi^2$={chi2[1]:.2f}', 
                 f'Blackbody $\chi^2$={chi2[2]:.2f}',  f'Low albedo & Specular $\chi^2$={chi2[3]:.2f}', f'High albedo & Specular $\chi^2$={chi2[4]:.2f}', 
                 f'$Al_2O_3$ surface $\chi^2$={chi2[5]:.2f}'], 
               loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=2)
    # plt.legend(plotarr, [ f'Low albedo $\chi^2$={chi2[0]:.2f}', f'High albedo $\chi^2$={chi2[1]:.2f}', 
    #              f'Blackbody $\chi^2$={chi2[2]:.2f}',
    #              f'$Al_2O_3$ surface $\chi^2$={chi2[3]:.2f}'], 
    #            loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=2)
    print(chi2)
    print(Bin)
    # Draw the wavelength range of the instrument
    Lower_bound = 2.4
    Upper_bound = 4.22
    ax.axvspan(Lower_bound, Upper_bound, color='gray', alpha=0.2)
    plt.text((Lower_bound + Upper_bound) /2, 147, 'NIRCam/F322W2', fontsize = 9, ha='center')
    plt.text((Lower_bound + Upper_bound) /2, 135, f'{Lower_bound:.2f}-{Upper_bound:.2f} μm', fontsize = 9, ha='center')
    
    # 添加标题和标签  
    # plt.title('')  
    plt.xlabel('Wavelength (μm)',  fontsize = 12)  
    plt.ylabel('Secondary eclipse depth (ppm)', fontsize = 12) 
    plt.axis([0.5, 12, 0, 160])
    
    plt.savefig('telescope_measure/data_plot12.pdf', format = 'pdf') 

    # 显示图表  
    plt.show()
    plt.close()
   
def real_comp2(name_list):
    # load data, I_diffuse,I_specular is contrast ratio 
    ############ load real measured data  ###############
    data_value = np.loadtxt(f"telescope_measure/JWST_data.txt", delimiter=',')
    data_ebar = np.loadtxt(f"telescope_measure/JWST_errorbar.txt", delimiter=',')
    
    DN = data_value.shape[0]
    plotarr  = [0] * 100
    
    data = np.zeros([DN, 3])
    data[:,0:2] = data_value
    for i in range(DN):
        data[i, 2] = np.abs(data_ebar[i*2,1] - data_ebar[i*2+1,1])/2
        
    # 创建散点图
    save_txt(data)  
    fig, ax = plt.subplots()
    plotarr[4] = plt.errorbar(data[:,0], data[:,1], yerr=data[:,2], fmt='ok', ecolor='k', linestyle='None', label = 'JWST')  
    
    # 处理离群点，计算chi2是否要去除离群点，store the new data after removing outliers in data2
    DN = DN # the size of data2
    Data2 = np.zeros([DN-1, 3])
    Data2[0:6] = data[0:6]  # the data point before outlier
    Data2[6:] = data[7:]  # the data point after outlier
    
    
    ## 处理模型数据并绘图
    pallet = ['red', 'orange',  'black', '#ffff00'] # color library
    chi2 = np.zeros(10)  # the chi2 of JWST MIRI data
    chi2_sim = np.zeros(10) # the chi2 of simulated data
    
    # processing error bar
    # sim_obs = np.loadtxt(f"telescope_measure/sim_obs (24).txt", delimiter=' ')
    sim_obs = np.loadtxt(f"telescope_measure/sim_obs_fin.txt", delimiter=' ')
    x = sim_obs[:,0]
    y = sim_obs[:,3] *1e6
    Nbin = 6 # the number of bin: 6
    Bin = np.zeros(Nbin)  
    xloc = np.zeros(Nbin)
    NGroup = sim_obs.shape[0] // Nbin # the number of pixels in each bin
    for i in range(Nbin):  # calculate the errorbar of each bins
        xq = x[NGroup*i:NGroup*(i+1)]  # group pixel points to each bin
        yq = y[NGroup*i:NGroup*(i+1)]
        Bin[i] = np.sqrt(1/np.sum(1/yq**2)) # calculate the errorbar of each bins
        xloc[i] = np.mean(xq) # calculate the mean wavelength of each bin
        
    Bin = Bin * np.sqrt(2)    # multiply the error by sqrt(2) (comparing in-eclipse vs out-eclipse)
    print([f'Pixels: {sim_obs.shape[0]}, bins number: {Nbin}, pixels number in each bin: {NGroup}'])
    print('Error bar:', Bin)
    
    
    for i, name in enumerate(name_list):  # Plot the spectra of each dataset
        I_diffuse = np.load(f'temp/{name}/variables/I_diffuse.npy')
        I_specular = np.load(f'temp/{name}/variables/I_specular.npy')
        Star_flux = np.load(f'temp/{name}/variables/Star_flux.npy')
        wave_list = np.load(f'temp/{name}/variables/wave_list.npy')
        Thermal = np.load(f'temp/{name}/variables/Thermal.npy')

        # convert I_diffuse,I_specular (contrast ratio) to absolute intensity
        N1 = Star_flux.shape[0]
        N2 = Star_flux.shape[1]
        ID = I_diffuse[:, N2//2]
        IS = I_specular[:, N2//2]
        Star_flux = Star_flux[:, N2//2]
        IT = Thermal[:, N2//2]
        
        if i == 0:
            data1 = np.zeros([np.size(wave_list), 3])
            data1[:,0] = wave_list *1e6
            data1[:,1] = (IT + IS) *1e6
            data1[:,2] = (IT + ID) *1e6
            # save_txt(data, "Wavelength(micron)\tFp/F*(Low albedo_specular ppm)\tFp/F*(Low albedo_Lambert ppm)", "Low albedo.txt")
        elif i == 1:
            data2 = np.zeros([np.size(wave_list), 3])
            data2[:,0] = wave_list *1e6
            data2[:,1] = (IT + IS) *1e6
            data2[:,2] = (IT + ID) *1e6
            # save_txt(data, "Wavelength(micron)\tFp/F*(High albedo_specular ppm)\tFp/F*(High albedo_Lambert ppm)", "High albedo.txt")
        elif i == 2:
            data3 = np.zeros([np.size(wave_list), 3])
            data3[:,0] = wave_list *1e6
            data3[:,1] = (IT + IS) *1e6
            data3[:,2] = (IT + ID) *1e6
            # save_txt(data, "Wavelength(micron)\tFp/F*(Blackbody ppm)\tFp/F*(Blackbody ppm)", "Blackbody.txt")
        # elif i == 3:
        #     data4 = np.zeros([np.size(wave_list), 3])
        #     data4[:,0] = wave_list *1e6
        #     data4[:,1] = (IT + IS) *1e6
        #     data4[:,2] = (IT + ID) *1e6
            # save_txt(data, "Wavelength(micron)\tFp/F*(Blackbody & heat fully redistribution ppm)\tFp/F*(Blackbody & heat fully redistribution ppm)", "Heat.txt")
        if i == 2:
            continue
    
        if i != 2:
            # plotarr[i], = plt.plot(wave_list *1e6, (IT + ID) *1e6, '--', color = pallet[i]) # lambert diffuse case
            plotarr[i], = plt.plot(wave_list *1e6, (IT + IS) *1e6, '-', color = pallet[i]) # specular case
        else:  # blackbody
            plotarr[i], = plt.plot(wave_list *1e6, (IT + ID) *1e6, '-', color = pallet[i]) # lambert diffuse case
        
        if i == 0:
            spl = interp1d(wave_list *1e6,  (IT + IS) *1e6, kind='linear')
            yloc = spl(xloc)
            # plt.errorbar(xloc, yloc, yerr=Bin, fmt='o', ecolor=pallet[i], linestyle='None')
            # Yloc = yloc + Bin* (np.random.random(np.size(Bin)) - 0.5) *2
            # Yloc = np.array([26.58986268 ,31.67471026, 37.53136225, 39.48648747, 39.86244592, 41.44114732])
            # Yloc = np.array([ 33.0184, 33.9478, 39.8372, 36.0629, 45.2019, 46.1313])
            # Bin = np.ones(np.size(Bin)) * 10
            Yloc = np.zeros(np.size(Bin))
            for k, yl in enumerate(yloc):
                Yloc[k] = yl + np.random.normal(loc=0, scale=1, size=1) * Bin[k]
            plotarr[5] = plt.errorbar(xloc, Yloc, yerr=Bin, fmt='o',color=pallet[i], ecolor=pallet[i], linestyle='None', markersize=5, elinewidth=1.5)
        
        # chi2[i] = chi2_cal(Data2[:,0], Data2[:,1], Data2[:,2], wave_list *1e6, (IT + ID) *1e6)
        chi2[i] = chi2_cal(Data2[:,0], Data2[:,1], Data2[:,2], wave_list *1e6, (IT + IS) *1e6)
        chi2_sim[i] = chi2_cal(xloc, Yloc, Bin, wave_list *1e6, (IT + IS) *1e6)
        
        
    data = np.zeros([np.size(wave_list), 7])
    data[:,0] = data3[:,0]
    data[:,1] = data3[:,1]
    data[:,2] = data1[:,1]
    data[:,3] = data1[:,2]
    data[:,4] = data2[:,1]
    data[:,5] = data2[:,2]
    # data[:,6] = data4[:,1]
    name = "Wavelength(micron)\tFp/F*(Blackbody ppm)\tFp/F*(Low albedo_specular ppm)\tFp/F*(Low albedo_Lambert ppm)\tFp/F*(High albedo_specular ppm)\tFp/F*(High albedo_Lambert ppm)\tFp/F*(Blackbody & heat fully redistribution ppm)"
    save_txt(data, name, "Model_result.txt")
    
    # draw Al2O3 surface
    GJAl2O3 = np.loadtxt('GJAl2O3.txt', skiprows = 1)
    # print(GJAl2O3)
    wave_GJ = GJAl2O3[:,0]
    Al2O3_GJ = GJAl2O3[:,1]
    Black_GJ = GJAl2O3[:,2]
    
    plotarr[3], = plt.plot(wave_GJ, Al2O3_GJ, label = 'Al_2O_3 surface', color = '#89cff0')
    # plotarr[-1], = plt.plot(wave_GJ, Black_GJ, label = 'Blackbody (XT)')
    chi2[3] = chi2_cal(Data2[:,0], Data2[:,1], Data2[:,2], wave_GJ, Al2O3_GJ)
    # chi2[-1] = chi2_cal(Data2[:,0], Data2[:,1], Data2[:,2], wave_GJ, Black_GJ)
    chi2_sim[3] = chi2_cal(xloc, Yloc, Bin, wave_GJ, Al2O3_GJ)
    
    # Carbon_W = np.loadtxt('GJAl2O3_Carbon_weathering.txt', skiprows = 1, delimiter = '\t')
    # Iron_W = np.loadtxt('GJAl2O3_Iron_weathering.txt', skiprows = 1, delimiter = '\t')
    # lam = Carbon_W[:,0]
    
    # carb_w = Carbon_W[:,7]
    # chi2[4] = chi2_cal(Data2[:,0], Data2[:,1], Data2[:,2], lam, carb_w)
    # chi2_sim[4] = chi2_cal(xloc, Yloc, Bin, lam, carb_w)
    # plotarr[4], = plt.plot(lam, carb_w, label = f'Carbon weathering 5', linestyle = '-' ,color ='#1a8fc6')
    # for i in range(6):
    #     carb_w = Carbon_W[:,i+2]
    #     chi2[i+4] = chi2_cal(Data2[:,0], Data2[:,1], Data2[:,2], lam, carb_w)
    #     plotarr[i+4], = plt.plot(lam, carb_w, label = f'Carbon weathering {i}', linestyle = '--')
        
    # for i in range(6):
    #     iron_w = Iron_W[:,i+2]
    #     chi2[i+10] = chi2_cal(Data2[:,0], Data2[:,1], Data2[:,2], lam, iron_w)
    #     plotarr[i+10], = plt.plot(lam, iron_w, label = f'Iron weathering {i}', linestyle = '-.')
        
        
    fig.subplots_adjust(bottom=0.3)  
    # plt.legend(plotarr, [f'Low albedo & Lambert $\chi^2$={chi2[0]:.2f}', f'Low albedo & Specular $\chi^2$={chi2[1]:.2f}', 
    #             f'High albedo & Lambert $\chi^2$={chi2[2]:.2f}', f'High albedo & Specular $\chi^2$={chi2[3]:.2f}', 
    #             'Blackbody & Lambert', 'Blackbody & Specular', 'Full redistribution & Lambert', 'Full redistribution & Specular',
    #             'O-H & Lambert', 'O-H & Specular'], loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2)
    # plt.legend(plotarr, ['Low albedo & Lambert', 'Low albedo & Specular',  'Mid albedo & Lambert', 'Mid albedo & Specular',
    #             'High albedo & Lambert', 'High albedo & Specular'], loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2)
    # plt.legend(plotarr, [ f'Low albedo & Lambert $\chi_{{nu}}^2$={chi2[0]:.2f}', f'High albedo & Lambert $\chi_{{nu}}^2$={chi2[1]:.2f}', 
    #              f'Blackbody $\chi_{{nu}}^2$={chi2[2]:.2f}',  f'Low albedo & Specular $\chi_{{nu}}^2$={chi2[3]:.2f}', f'High albedo & Specular $\chi_{{nu}}^2$={chi2[4]:.2f}', 
    #              f'$Al_2O_3$ surface $\chi_{{nu}}^2$={chi2[5]:.2f}'], 
    #            loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=2)
    # plt.legend(plotarr, [ f'Low albedo $\chi_{{nu}}^2$={chi2[0]:.2f}', f'High albedo $\chi_{{nu}}^2$={chi2[1]:.2f}', 
    #              f'Blackbody $\chi_{{nu}}^2$={chi2[2]:.2f}',
    #              f'$Al_2O_3$ surface $\chi_{{nu}}^2$={chi2[3]:.2f}', 'Carbon0.1%','Carbon0.2%','Carbon0.5%','Carbon1%','Carbon2%','Carbon5%',
    #              'Iron0.1%','Iron0.2%','Iron0.5%','Iron1%','Iron2%','Iron5%'], 
    #            loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=4)
    # plt.legend([plotarr[0],  plotarr[2], plotarr[4], plotarr[1], plotarr[3]], [ f'Low albedo $\chi_{{nu}}^2$={chi2[0]:.2f}',  f'Blackbody $\chi_{{nu}}^2$={chi2[2]:.2f}', 
    #             f'Space weathered $\mathrm{{Al_2O_3}}$ $\chi_{{nu}}^2$={chi2[4]:.2f}', f'High albedo $\chi_{{nu}}^2$={chi2[1]:.2f}', f'Pure $\mathrm{{Al_2O_3}}$ $\chi_{{nu}}^2$={chi2[3]:.2f}' ], 
    #            loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=2)
    # plt.legend([plotarr[0],  plotarr[2], plotarr[4], plotarr[1], plotarr[3]], [ f'Low albedo $\chi^2$={chi2[0]:.2f}',  f'Blackbody $\chi^2$={chi2[2]:.2f}', 
    #             f'Space weathered $\mathrm{{Al_2O_3}}$ $\chi^2$={chi2[4]:.2f}', f'High albedo $\chi^2$={chi2[1]:.2f}', f'Pure $\mathrm{{Al_2O_3}}$ $\chi^2$={chi2[3]:.2f}' ], 
    #            loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=2)
    # plt.legend([plotarr[0],  plotarr[2], plotarr[4], plotarr[1], plotarr[3]], [ f'Low albedo',  f'Blackbody', 
    #             'Space weathered $\mathrm{Al_2O_3}$', f'High albedo', 'Pure $\mathrm{Al_2O_3}$' ], 
    #            loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=2)
    # lgd = plt.legend([plotarr[0], plotarr[1], plotarr[3], plotarr[2]], [ f'Low albedo lava', f'High albedo lava', '$\mathrm{Al_2O_3}$',  f'Blackbody' ], 
    #            loc='lower right', ncol=2, fontsize = 10.5)
    lgd = plt.legend([plotarr[5], plotarr[0], plotarr[1], plotarr[4], plotarr[3]], 
                     ['NIRCam data', f'Low albedo lava', f'High albedo lava',  'MIRI data',  '$\mathrm{Al_2O_3}$'], 
               loc='lower right', bbox_to_anchor=(1.01, -0.01), ncol=2, fontsize = 10.5)
    # lgd.get_frame().set_visible(False)
    
    print('chi2 of JWST MIRI data:', chi2)
    print('chi2 of simulated data:', chi2_sim)
    print('reduced chi2 of JWST MIRI data:', chi2/11)
    print('reduced chi2 of simulated data:', chi2_sim/6)
    print('Bin', Bin)
    print('xloc:',xloc)
    print('Yloc:',Yloc)
    # Draw the wavelength range of the instrument
    Lower_bound = 2.4
    Upper_bound = 4.22
    ax.axvspan(Lower_bound, Upper_bound, color='gray', alpha=0.2)
    plt.text((Lower_bound + Upper_bound) /2, 100, 'NIRCam/F322W2', fontsize = 13, ha='center')
    plt.text((Lower_bound + Upper_bound) /2, 89, f'{Lower_bound:.2f}-{Upper_bound:.2f} μm', fontsize = 13, ha='center')
    
    # 添加标题和标签  
    # plt.title('')  
    plt.xlabel('Wavelength (μm)',  fontsize = 13)  
    plt.ylabel('Secondary eclipse depth (ppm)', fontsize = 13) 
    plt.axis([0.5, 12, 0, 110])
    
    plt.savefig('telescope_measure/data_plot17.pdf', format = 'pdf') 

    # 显示图表  
    plt.show()
    plt.close()
 
def Intenstiy_comp(name_list):
    for i, name in enumerate(name_list):  # Plot the spectra of each dataset
        I_diffuse = np.load(f'temp/{name}/variables/I_diffuse.npy')
        I_specular = np.load(f'temp/{name}/variables/I_specular.npy')
        Star_flux = np.load(f'temp/{name}/variables/Star_flux.npy')
        wave_list = np.load(f'temp/{name}/variables/wave_list.npy')
        Thermal = np.load(f'temp/{name}/variables/Thermal.npy')

        # convert I_diffuse,I_specular (contrast ratio) to absolute intensity
        N1 = Star_flux.shape[0]
        N2 = Star_flux.shape[1]
        if i == 0:
            ID1 = I_diffuse[:, N2//2]
            IS1 = I_specular[:, N2//2]
            IT1 = Thermal[:, N2//2]
        elif i == 1:
            ID2 = I_diffuse[:, N2//2]
            IS2 = I_specular[:, N2//2]
            IT2 = Thermal[:, N2//2]          
        
    res = np.zeros([np.size(ID1),2])
    res[:, 0] = np.transpose(wave_list)
    res[:, 1] = IS1 + IT1 - IS2 - IT2
    print(IS1)
    
    np.save(f"temp/{name}/res.npy", res)
    
def phase_curve_plot_withdata(name_list, wave_range, instrument = '  ', model = ' '):
    '''
    绘制phase curve, 并绘制真实datapoint以及errorbar, 并计算chi2 
    # 功能一定程度上与compare_phase_curve_plot类似, 但包含了对data的处理和绘制
    name_list: 可以为一个或两个name; 
        if len(name_list) == 1, <model> is needed to determine albedo model
        if len(name_list) == 2, first is low albedo, second is high albedo
    instrument [wave_range]: For K2-141b, only 'Kepler' and 'Spitzer' are available
        Kepler: [0.43, 0.89]
        Spitzer: [4, 5]
        
    model: 'low' or 'high'. Only len(name_list) == 1 available (仅当只有一个name_list对象时可用,用来确定name的albedo模型)
    '''
    ### load data points and errorbar
    data = np.loadtxt(f'telescope_measure/K2-141b_{instrument}.txt', delimiter=',')
    data_x = data[:,0]
    data_y = data[:,1]
    # data_err = data[:,2]  # load method 1
    # data_err = np.loadtxt(f'telescope_measure/K2-141b_{instrument}_err.txt', delimiter=',')  # load method 2
    if instrument == 'Kepler':  # load method 3
        data_err = np.ones(np.size(data_x)) * 13.5/2  
    elif instrument == 'Spitzer':
        data_err = np.ones(np.size(data_x)) * 75.4/2
    
    pallet = ['b','r','k']
    plotarr  = [0] * 8
    chi2_array = np.zeros([8])
    fig, ax = plt.subplots()
    
    # plt.plot(data_x, data_y,'.k')
    ax.errorbar(data_x, data_y, yerr=data_err, fmt='o', color='black', markersize=4)
    # Period: 7.72614 h
    Theta_list = np.load(f'temp/{name_list[0]}/variables/Theta.npy')
    Theta_list = Theta_list / (2 *np.pi)   # 将相位角归一化
    data = np.zeros([np.size(Theta_list), 6])
    
    up_bound = 0  # control xlim up_lim
    for i, name in enumerate(name_list):
        Theta_list = np.load(f'temp/{name}/variables/Theta.npy')
        Theta_list = Theta_list / (2 *np.pi)   # 将相位角归一化
        Nt = np.size(Theta_list)
        CR_S = np.zeros([Nt])
        CR_D = np.zeros([Nt])
        for j, theta in enumerate(Theta_list):
            CR_S[j], CR_D[j] = bond_albedo_calculator(wave_range[0], wave_range[1], name, j)
        
        plotarr[i*2], = plt.plot(Theta_list, CR_D, '-', color = pallet[i], linewidth = 2)
        plotarr[i*2+1], = plt.plot(Theta_list, CR_S, '--', color = pallet[i], linewidth = 2)
        data[:,0] = Theta_list
        up_bound = np.max([np.max(CR_D), np.max(CR_S), up_bound]) # find the max value for ylim
        if i != 2:
            data[:,2 + i*2] = CR_S
            data[:,3 + i*2] = CR_D
        elif i == 2:
            data[:,1] = CR_S
        
        # calculate chi2
        chi2_array[i*2] = chi2_cal(data_x, data_y, data_err ,Theta_list, CR_D)
        chi2_array[i*2+1] = chi2_cal(data_x, data_y, data_err ,Theta_list, CR_S)
        
    fig.subplots_adjust(bottom=0.25)  # adjust the bottom margin to make room for the legend
    plt.ylim([0,np.max([up_bound, np.max(data_y)]) * 1.1])
    plt.xlim([0,1])
    # set legend based on the size of name_list
    sc = r'$\chi_{\nu}^2$'  # set the chi2 symbol: chi2 (r'$\chi^2$') or reduced chi2 (r'$\chi_{\nu}^2$')
    if len(name_list) == 1:
        plt.legend([plotarr[0],plotarr[1]], ['Lambert, '+sc+f'={chi2_array[0]:.2f}', 'Specular, '+sc+f'={chi2_array[1]:.2f}'])
    elif len(name_list) == 2:
        plt.legend([plotarr[0],plotarr[2],plotarr[1],plotarr[3]], ['Low albedo & Lambert, '+sc+f'={chi2_array[0]:.2f}', 'High albedo & Lambert, '+sc+f'={chi2_array[2]:.2f}', 'Low albedo & Specular, '+sc+f'={chi2_array[1]:.2f}', 'High albedo & Specular, '+sc+f'={chi2_array[3]:.2f}'], loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2)
            
    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2)  
    plt.xlabel('Orbital Phase', fontsize = 13)
    plt.ylabel(r'$F_p/F_*$ (ppm)', fontsize = 13)
    if len(name_list) == 1:
        plt.title(f'{instrument}: {wave_range[0]*1e6:.2f}-{wave_range[1]*1e6:.2f} μm, {model} albedo model')
    elif len(name_list) == 2:
        plt.title(f'{instrument}: {wave_range[0]*1e6:.2f}-{wave_range[1]*1e6:.2f} μm')
    # 添加文字和箭头，设置字体透明度  

    # plt.savefig(f"temp/{name}/phase_curve_comp1.pdf", format = 'pdf')
    for name in name_list:
        plt.savefig(f"temp/{name}/phase_curve_withdata_{instrument}.pdf", format = 'pdf')
    plt.show()
    plt.close()
    

if __name__ =='__main__':
    '''
    instrument list:
    CHEOPS: [0.33, 1.1]
    HST/WFC3/G102: [0.80, 1.15]
    HST/WFC3/G141: [1.075, 1.70]
    JWST/NIRCam/F322W2: [2.7, 4.0]
    Kepler: [0.43, 0.89]
    Spitzer: [4, 5]
    '''   
    # first: low albedo ; second: high albedo
    # compare_phase_curve_plot(['Low_copy', 'High_copy'], np.array([0.33, 1.1])* 1e-6, instrument = 'CHEOPS', legend = 'insert', xlabel = 'off', ylabel='on', errorbar=38.73)
    # compare_phase_curve_plot(['Low_copy', 'High_copy'], np.array([0.80, 1.15])* 1e-6, instrument = 'HST/WFC3/G102', legend = 'off', xlabel='off', ylabel='off', errorbar=7.40)
    # compare_phase_curve_plot(['Low_copy', 'High_copy'], np.array([1.075, 1.70])* 1e-6, instrument = 'HST/WFC3/G141', legend = 'off', xlabel = 'on', ylabel='on', errorbar=6.85)
    # compare_phase_curve_plot(['Low_copy', 'High_copy'], np.array([2.7, 4.0])* 1e-6, instrument = 'JWST/NIRCam/F322W2', legend = 'off', xlabel = 'on', ylabel='off', errorbar=5.14)
    
    phase_curve_plot_withdata(['Low_copy', 'High_copy'], np.array([0.43, 0.89])* 1e-6, instrument = 'Kepler')
    phase_curve_plot_withdata(['Low_copy', 'High_copy'], np.array([4, 5])* 1e-6, instrument='Spitzer')
    
    # phase_curve_plot_withdata(['R1copy'], np.array([0.43, 0.89])* 1e-6, instrument = 'Kepler', model = 'Low')
    # phase_curve_plot_withdata(['R2copy'], np.array([0.43, 0.89])* 1e-6, instrument = 'Kepler', model='High')
    # phase_curve_plot_withdata(['R3copy'], np.array([4, 5])* 1e-6, instrument='Spitzer', model ='Low')
    # phase_curve_plot_withdata(['R4copy'], np.array([4, 5])* 1e-6, instrument='Spitzer', model ='High') 
    
    