import pandas as pd  
from scipy.constants import h, c, k
import numpy as np

def B_func(T, lam, B=1): # palnk function
    return 2 * B * h * c**2 / lam**5 / (np.exp(h * c / (lam * k * T)) - 1)

def Rank_exo(lam_spec = 1e-6, T_liq = 1600, Nhead = 30, rank_standard = 'Specular_corr' ): 
    '''
    pack exo_rank.py into a function, rank exoplanets by Specular/Thermal at lam_spec
    lam_spec: wavelength of observation
    T_liq: liquid temperature, 作为完全融化区域的温度阈值
    # Nname: number of names of planets
    Nhead: number of planets in the head
    rank_standard: 按照什么参数排序，'Specular_corr', 'Specular', 'Spl_CR' etc.
    '''

    R_E = 6357
    AU = 149_597_870.7 
    Co1 = (R_E/AU)**2
    R_Sun = 695500


    # 读取CSV文件  
    df = pd.read_csv('PS.csv', header=96) 
    df_c = df.copy() 
    #print(df.columns)

    # 确保R、a和name列存在  
    if 'pl_orbsmax' in df_c.columns and 'pl_rade' in df_c.columns and 'pl_name' in df_c.columns:  
        # 计算R/a  
        fl_df1 = df_c[(df_c['pl_rade'] < 1.6)].copy()
        fl_df1['Tsub'] = fl_df1['st_teff'] * np.sqrt(fl_df1['st_rad'] *R_Sun /fl_df1['pl_orbsmax'] /AU)  # sub-stellar temperature
        fl_df1['liq_area'] = 1- (T_liq/fl_df1['Tsub'])**8 # 计算熔融覆盖率
        fl_df = fl_df1[(fl_df1['Tsub'] > T_liq) & (fl_df1['liq_area'] > 0.9)].copy() # 选择Tsub > T_liq, 并且熔融覆盖率>80%的行

        # fl_df['RSM'] = (fl_df['pl_rade'] / fl_df['pl_orbsmax'])**2 * Co1 * 10**(-fl_df['sy_kmag']/5) 

        # fl_df['RSM2'] = (fl_df['pl_rade'] / fl_df['pl_orbsmax'] /2)**2 * Co1 *1e6
        
        fl_df['Specular']  = B_func(fl_df['st_teff'], lam_spec)/B_func(fl_df['pl_eqt'] * (8/3)**0.25, lam_spec) * (fl_df['st_rad'] /fl_df['pl_orbsmax']/2) **2 * (R_Sun/AU)**2
        # 根据R/a降序排序  
        fl_df['Specular_corr'] = fl_df['Specular']  * fl_df['liq_area'] # 融化区域修正后的Specular
        
        fl_df['Spl_CR'] = (fl_df['pl_rade']/fl_df['pl_orbsmax']/2 )**2 *Co1 *1e6 # calculate the ratio of specular F_{specular}/F_*
        sorted_df = fl_df.sort_values(by = rank_standard, ascending=False) 
        
        # 提取前N个恒星的name参数  
        if Nhead ==0:
            top_names = sorted_df[['pl_name','Specular', 'Specular_corr', 'Spl_CR']]
        else:
            top_names = sorted_df[['pl_name','Specular', 'Specular_corr','st_rad', 'pl_orbsmax', 'st_teff', 'pl_eqt', 'Tsub', 'Spl_CR', 'liq_area']].head(Nhead)  
        
        # for i in range(len(top_10_names)):  
        #     print(f" {top_10_names.iloc[i]['pl_name']}")
        return top_names
    else:  
        print("CSV文件中缺少R、a或name列")
        