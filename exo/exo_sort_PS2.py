import pandas as pd  
import numpy as np
from scipy.constants import h, c, k
from scipy.integrate import quad

# 计算光子数N，并计算信噪比1/sqrt(N)
# 常量
R_earth = 6371
R_Jupiter = 69911
R_Sun = 695500
AU = 149597870.7
R_E = 6357
AU = 149_597_870.7 
Co1 = (R_E/AU)**2

Ag = 0.3
kb = 1.38e-23
NA = 6.022e23

# Planck函数
def B_func(T, lam, B=1):
    return 2 * B * h * c**2 / lam**5 / (np.exp(h * c / (lam * k * T)) - 1)

# 积分函数
def Int_func(lam, Ts):
    return B_func(Ts, lam, 1) * lam

# 读取CSV文件  
df = pd.read_csv(r'PS.csv', header=96) 
df_c = df.copy() 
#print(df.columns)

# 确保R、a和name列存在  
if 'pl_orbsmax' in df_c.columns and 'pl_rade' in df_c.columns and 'pl_name' in df_c.columns:  
    # 计算R/a  
    fl_df = df_c[(df_c['pl_rade']<1.6) & (df_c['pl_eqt']>850)]
    # 反射对比度
    fl_df['Reflected_contrast'] = Ag * (fl_df['pl_rade'] / fl_df['pl_orbsmax']) ** 2 * Co1

    # 噪声估计相关参数
    pc = 3.26 * 9.461e12
    D = np.array([6.5, 2.4, 6.5, 0.85, 6.5, 6.5, 0.85, 6.5])
    lam1 = np.array([0.6, 1.1, 1.0, 3.225, 3.12, 9, 21.65, 23.5]) * 1e-6
    lam2 = np.array([1, 1.7, 4.0, 3.975, 4.01, 11, 26.35, 27.5]) * 1e-6
    tal = np.array([0.3, 0.4, 0.4, 0.4, 0.4, 0.36, 0.45, 0.18])
    Dt = 1 * 3600
    d = 20 * pc
    INS = 0

    fl_df['Noise'] = fl_df['Reflected_contrast'] *0
    # 计算光子数和噪声
    for i,Ts in enumerate(fl_df.iloc[:, 'Noise']):
        fl_df.iloc[i, 'Noise'], _ = quad(Int_func, lam1[i], lam2[i], args=(Ts,))
        
    fl_df['Noise'] = np.pi / h / c * tal[i] * Dt * (Rs * D[i] / 2 / d) ** 2 * fl_df['Noise']
    fl_df['Noise'] = 1 / np.sqrt(fl_df['Noise']) / np.sqrt(np.pi)

    fl_df['RMS'] = fl_df['Reflected_contrast'] / fl_df['Noise']

    # fl_df['RSM'] = (fl_df['pl_rade'] / fl_df['pl_orbsmax'])**2 * Co1 * 10**(-fl_df['sy_kmag']/5)

    # fl_df['RSM2'] = (fl_df['pl_rade'] / fl_df['pl_orbsmax'])**2 * Co1 
    
    # 根据R/a降序排序  
    sorted_df = fl_df.sort_values(by='RSM', ascending=False)  
    
    # 提取前10个恒星的name参数  
    top_10_names = sorted_df[['pl_name','Reflected_contrast', 'Noise', 'RMS']].head(10)  
    
    # 输出前10个恒星的name  
    print(top_10_names)  
else:  
    print("CSV文件中缺少R、a或name列")  