import pandas as pd  
from scipy.constants import h, c, k
import numpy as np

def B_func(T, lam, B=1): # palnk function
    return 2 * B * h * c**2 / lam**5 / (np.exp(h * c / (lam * k * T)) - 1)

# rank exoplanets by Specular/Thermal at lam_spec
lam_spec = 4e-6
R_E = 6357
AU = 149_597_870.7 
Co1 = (R_E/AU)**2
R_Sun = 695500
R_Earth = 6371
M_Earth = 5.972e24  
M_Sun = 1.989e30
G = 6.67430e-11

T_liq = 850 *(8/3)**0.25 # liquid temperature, 作为完全融化区域的温度阈值
# 读取CSV文件  
df = pd.read_csv('PS.csv', header=96) 
df_c = df.copy() 
#print(df.columns)

# 确保R、a和name列存在  
if 'pl_orbsmax' in df_c.columns and 'pl_rade' in df_c.columns and 'pl_name' in df_c.columns:  
    # 计算R/a  
    fl_df = df_c[(df_c['pl_rade']<1.6)].copy()
    fl_df = fl_df[(fl_df['pl_eqt']>T_liq)].copy()
    a2 = fl_df['pl_orbsmax'] * (fl_df['pl_eqt']/T_liq) ** 2
    P2 = fl_df['pl_orbper'] * (a2/fl_df['pl_orbsmax'])**1.5
    fl_df['P2'] = P2
    
    t_lock = (AU *fl_df['pl_orbsmax'])**6 /(M_Earth *fl_df['pl_bmasse'])**2 *(R_Earth *fl_df['pl_rade']) /(M_Sun *fl_df['st_mass']) *6e10 *3e10 *1e21
    fl_df['t_lock'] = t_lock
    
    # fl_df = fl_df[(fl_df['P2']<10)].copy()
    sorted_df = fl_df.sort_values(by='t_lock', ascending=False) 
    
    # 提取前10个恒星的name参数  
    top_10_names = sorted_df[['pl_name','P2', 't_lock', 'pl_eqt']].head(100)  
    
    # 输出前10个恒星的name  
    pd.set_option('display.max_rows', None)
    print(top_10_names)
else:  
    print("CSV文件中缺少R、a或name列")