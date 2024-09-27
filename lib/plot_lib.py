import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os
import argparse
from bond_albedo_calculator import bond_albedo_calculator
from function_library import chi2_cal

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
    spls = interp1d(Wave_list, I_specular.T, kind='cubic', fill_value='extrapolate')
    spli = interp1d(Wave_list, I_intensity.T, kind='cubic', fill_value='extrapolate')
    spld = interp1d(Wave_list, I_diffuse.T, kind='cubic', fill_value='extrapolate')
    splt = interp1d(Wave_list, Thermal.T, kind='cubic', fill_value='extrapolate')
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
    spls = interp1d(Theta_list, I_s, kind='cubic', fill_value='extrapolate')
    spli = interp1d(Theta_list, I_i, kind='cubic', fill_value='extrapolate')
    spld = interp1d(Theta_list, I_d, kind='cubic', fill_value='extrapolate')
    splt = interp1d(Theta_list, I_t, kind='cubic', fill_value='extrapolate')
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

    # plt.xlabel('Orbital angle (rad)')
    # plt.ylabel('Contrast ratio (ppm)')
    # plt.legend()
    # plt.show()
    # os.makedirs(f'temp/P0', exist_ok= True)
    # plt.savefig(f'temp/P0/compare.png')
    # plt.close()

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
    ax.plot(theta, I_s[i,:] * 1e6, label='Specular', color='k')
    ax.set_ylim(ymin=0, ymax= np.max(I_s[i,:])*1.1e6)
    ax.set_xlabel('Orbital angle (rad)', fontsize=18)
    ax.set_ylabel('Contrast ratio (ppm)', fontsize=18)
    ax.tick_params(length=6, width=2)
    ax.spines['right'].set_visible(False)

    lambda_color = 'blue'
    labmda_ax = ax.twinx()
    labmda_ax.set_position(axpos)
    labmda_ax.plot(theta, I_d[i,:] * 1e6, label='Diffuse', color=lambda_color)
    labmda_ax.set_ylim(ymin=0, ymax= np.max(I_d[i,:])*1.1e6)
    labmda_ax.set_xlabel('Orbital angle (rad)', fontsize=18)
    labmda_ax.tick_params(length=6, width=2, color=lambda_color, labelcolor=lambda_color)
    labmda_ax.set_ylabel('Contrast ratio (ppm)', fontsize=18, color=lambda_color)
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
    omglog_ax.spines['right'].set_position(('data', np.max(theta)*1.15))
    omglog_ax.set_ylim(0, np.max(I_t[i,:])*1.1e6)
    omglog_ax.set_position(axpos)
    omglog_ax.plot(theta, I_t[i,:] * 1e6, label='Thermal', color=omglog_color)
    omglog_ax.set_ylabel('Contrast ratio (ppm)', fontsize=18, color=omglog_color)
    omglog_ax.tick_params(length=6, width=2, color=omglog_color, labelcolor=omglog_color)
    omglog_ax.spines['right'].set(color=omglog_color, linewidth=2.0, linestyle='-.')
    fig.legend(['Specular', 'Diffuse', 'Thermal'], loc='upper left')
    plt.show()
    # os.makedirs(f'temp/P3', exist_ok= True)
    plt.savefig(f'temp/{name}/compare_{Obs_wave[0]*1e6}.png')
    plt.close()

    fig, ax = plt.subplots(figsize=(9,6))
    ax.plot(theta, (I_s[i,:]+I_t[i,:])*1e6 , label = 'Specular')
    ax.plot(theta, (I_d[i,:]+I_t[i,:])*1e6 , label = 'Diffuse')
    ax.set_xlabel('Orbital angle (rad)', fontsize=18)
    ax.set_ylabel('Contrast ratio (ppm)', fontsize=18)

    plt.errorbar(theta[N//2], (I_s[i,N//2]+I_t[i,N//2])*1e6, yerr = 1.25, capsize= 10)
    plt.plot(theta[N//2], (I_s[i,N//2]+I_t[i,N//2])*1e6,'.')

    plt.legend(['Specular', 'Diffuse'])
    plt.show()
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
    plt.xlabel('Wavelength ($\mu$m)')
    plt.ylabel('Contrast ratio (ppm)')
    
    plt.show()
    plt.savefig(f"temp/{name}/spectrum")
    
def compare_phase_curve_plot(name_list, wave_range):
    # load data, I_diffuse,I_specular is contrast ratio 
    sim_obs = np.loadtxt(f"telescope_measure/sim_obs (10).txt", delimiter=' ')
    x = sim_obs[:,0]
    y = sim_obs[:,3] *1e6
    Bin = np.sqrt(1/np.sum(1/y**2))
    
    pallet = ['b','r']
    plotarr  = [0] * 6
    fig, ax = plt.subplots()
    
    ebar = np.array([Bin, Bin * np.sqrt(1/0.629), Bin])
    
    xloc = np.array([0.3948, 0.5, 0.6052])
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
        
        spl = interp1d(Theta_list, CR_D, kind='linear')
        yloc = spl(xloc)
        plt.errorbar(xloc, yloc, yerr = 5, fmt='o', ecolor=pallet[i], linestyle='None')
        
        # 调整布局以便为图例腾出空间  
    fig.subplots_adjust(bottom=0.25)  

    print(ebar)
    # 设置背景颜色  
    ax.axvspan(0.4593, 0.5407, color='gray', alpha=0.5)
    ax.axvspan(0.3303, 0.4593, color='gray', alpha=0.2)
    ax.axvspan(0.5407, 0.6697, color='gray', alpha=0.2)
    # 设置图例并放置在图窗的正下方  
    plt.legend(plotarr, ['Low albedo & Lambert', 'Low albedo & Specular', 
                'High albedo & Lambert', 'High albedo & Specular'], loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2)
    # plt.legend(plotarr, ['Low albedo & Lambert', 'Low albedo & Specular',  'Mid albedo & Lambert', 'Mid albedo & Specular',
    #             'High albedo & Lambert', 'High albedo & Specular'], loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2)
    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2)  
    plt.xlabel('Orbital Phase')
    plt.ylabel('Contrast ratio (ppm)')
    # 添加文字和箭头，设置字体透明度  
    plt.text(0.3948, -5, '1 hour', fontsize=12, color='black', ha='center') 
    # plt.text(0.3948, -5, '', fontsize=12, color='black', ha='center')
    plt.text(0.5, -5, '0.63 hour', fontsize=12, color='black', ha='center')
    # plt.text(0.3948, -3, '1 hour', fontsize=12, color='black', ha='center')
    plt.text(0.6052, -5, '1 hour', fontsize=12, color='black', ha='center')
    # plt.text(0.3948, -3, '1 hour', fontsize=12, color='black', ha='center')
    
    plt.show()
    plt.savefig(f"temp/{name}/phase_curve_comp")
    plt.close()
    
    
   
def real_comp(name_list):
    # load data, I_diffuse,I_specular is contrast ratio 
    ############ load real measured data  ###############
    data_value = np.loadtxt(f"telescope_measure/JWST_data.txt", delimiter=',')
    data_ebar = np.loadtxt(f"telescope_measure/JWST_errorbar.txt", delimiter=',')
    
    DN = data_value.shape[0]
    
    data = np.zeros([DN, 3])
    data[:,0:2] = data_value
    for i in range(DN):
        data[i, 2] = np.abs(data_ebar[i*2,1] - data_ebar[i*2+1,1])/2
        
    # 创建散点图  
    fig, ax = plt.subplots()
    plt.errorbar(data[:,0], data[:,1], yerr=data[:,2], fmt='ok', ecolor='k', linestyle='None')  
    
    DN = DN -1
    data2 = np.zeros([DN, 3])
    data2[0:7] = data[0:7]
    data2[7:] = data[8:]
    
    
    ## 处理模型数据并绘图
    pallet = ['b', 'r']
    plotarr  = [0] * 6
    chi2 = np.zeros(np.size(name_list)*2)
    
    sim_obs = np.loadtxt(f"telescope_measure/sim_obs.txt", delimiter=' ')
    x = sim_obs[:,0]
    y = sim_obs[:,3] *1e6
    Bin = np.zeros(6)
    xloc = np.zeros(6)
    for i in range(6):
        xq = x[214*i:214*(i+1)]
        yq = y[214*i:214*(i+1)]
        Bin[i] = np.sqrt(1/np.sum(1/yq**2))
        xloc[i] = np.mean(xq)
    
  
    for i, name in enumerate(name_list):
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
    
        plotarr[i*2], = plt.plot(wave_list *1e6, (IT + ID) *1e6, '-', color = pallet[i])
        plotarr[i*2+1], = plt.plot(wave_list *1e6, (IT + IS) *1e6, '--', color = pallet[i])
        
        spl = interp1d(wave_list *1e6,  (IT + ID) *1e6, kind='linear')
        yloc = spl(xloc)
        plt.errorbar(xloc, yloc, yerr=Bin, fmt='o', ecolor=pallet[i], linestyle='None')
        
        chi2[i*2] = chi2_cal(data2[:,0], data2[:,1], data2[:,2], wave_list *1e6, (IT + ID) *1e6)
        chi2[i*2 + 1] = chi2_cal(data2[:,0], data2[:,1], data2[:,2], wave_list *1e6, (IT + IS) *1e6)
        
    fig.subplots_adjust(bottom=0.25)  
    plt.legend(plotarr, [f'Low albedo & Lambert $\chi^2$={chi2[0]:.2f}', f'Low albedo & Specular $\chi^2$={chi2[1]:.2f}', 
                f'High albedo & Lambert $\chi^2$={chi2[2]:.2f}', f'High albedo & Specular $\chi^2$={chi2[3]:.2f}'], loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2)
    # plt.legend(plotarr, ['Low albedo & Lambert', 'Low albedo & Specular',  'Mid albedo & Lambert', 'Mid albedo & Specular',
    #             'High albedo & Lambert', 'High albedo & Specular'], loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2)
    
    print(chi2)
    print(Bin)

    # 添加标题和标签  
    # plt.title('')  
    plt.xlabel('Wavelength ($\mu$m)')  
    plt.ylabel('Eclipse depth (ppm)') 
    plt.axis([0.5, 12, 0, 160])
    
    plt.savefig('telescope_measure/data_plot3.png') 

    # 显示图表  
    plt.show()
    plt.close()
    
    
if __name__ =='__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--file', default = 'R1', type = str)
    # args = parser.parse_args()
    # name = 'R0'
    # # IDS_plot(name, np.array([1]) * 1e-6)
    # wave_range = np.array([0.5, 5]) * 1e-6
    # spectrum_plot(name, wave_range)
    
    # compare_spectrum_plot(['R3', 'R4', 'R5'])
    # real_comp(['GJ-367 b low', 'GJ-367 b high'])
    
    compare_phase_curve_plot(['GJ-367 b PC low', 'GJ-367 b PC high'], np.array([2.87, 5.10])* 1e-6)
    
    
    
    