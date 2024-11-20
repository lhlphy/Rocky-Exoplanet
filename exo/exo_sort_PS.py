import pandas as pd  

R_E = 6357
AU = 149_597_870.7 
Co1 = (R_E/AU)**2

# 读取CSV文件  
df = pd.read_csv(r'C:\Users\dell\Desktop\Exo\PS.csv', header=96) 
df_c = df.copy() 
#print(df.columns)

# 确保R、a和name列存在  
if 'pl_orbsmax' in df_c.columns and 'pl_rade' in df_c.columns and 'pl_name' in df_c.columns:  
    # 计算R/a  
    fl_df = df_c[(df_c['pl_rade']<1.6) & (df_c['pl_eqt']>850)]
    # fl_df = df_c

    fl_df['RSM'] = (fl_df['pl_rade'] / fl_df['pl_orbsmax'])**2 * Co1 * 10**(-fl_df['sy_kmag']/5)

    fl_df['RSM2'] = (fl_df['pl_rade'] / fl_df['pl_orbsmax'])**2 * Co1 
    
    # 根据R/a降序排序  
    sorted_df = fl_df.sort_values(by='RSM2', ascending=False)  
    
    # 提取前10个恒星的name参数  
    top_10_names = sorted_df[['pl_name','RSM', 'RSM2']].head(10)  
    
    # 输出前10个恒星的name  
    print(top_10_names)  
else:  
    print("CSV文件中缺少R、a或name列")