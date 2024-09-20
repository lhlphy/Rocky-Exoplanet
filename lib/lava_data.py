# B´arðarbunga basalts, T=1787 K.txt
# MORB, T=861 K.txt
# Fasnia, T=1640 K.txt
# B-Nat, T=1648 K.txt
# B-glass, T=1654 K.txt

# CaFeSiO F3, T=1627 K.txt
# CaFeSiO F4, T=1673 K.txt
# CaFeSiO F5, T=1673 K.txt
# CaFeSiO F6, T=1673 K.txt
# CaFeSiO F7, T=1673 K.txt
# CaFeSiO F10, T=1673 K.txt
# CaFeSiO F11, T=1823 K.txt
# CaFeSiO F12, T=1823 K.txt

#命名中温度有问题


import numpy as np
from scipy.interpolate import interp1d  
import matplotlib.pyplot as plt  

def remove_duplicate(wave_combined, IS_combined):
    # 移除重复的波长值并取平均  
    unique_waves, indices = np.unique(wave_combined, return_index=True)  
    IS_combined_unique = [np.mean(IS_combined[wave_combined == wave]) for wave in unique_waves]  

    return unique_waves, IS_combined_unique

# if __name__ == '__main__':
#     # CaFeSiOF3 = np.loadtxt(f"lava_lib/CaFeSiO F3, T=1627 K.txt", delimiter=',')
#     CaFeSiOF4 = np.loadtxt(f"lava_lib/CaFeSiO F4, T=1673 K.txt", delimiter=',')
#     # CaFeSiOF5 = np.loadtxt(f"lava_lib/CaFeSiO F5, T=1673 K.txt", delimiter=',')
#     CaFeSiOF6 = np.loadtxt(f"lava_lib/CaFeSiO F6, T=1673 K.txt", delimiter=',')
#     # CaFeSiOF7 = np.loadtxt(f"lava_lib/CaFeSiO F7, T=1723 K.txt", delimiter=',')
#     # CaFeSiOF10 = np.loadtxt(f"lava_lib/CaFeSiO F10, T=1673 K.txt", delimiter=',')
#     # CaFeSiOF11 = np.loadtxt(f"lava_lib/CaFeSiO F11, T=1823 K.txt", delimiter=',')
#     # CaFeSiOF12 = np.loadtxt(f"lava_lib/CaFeSiO F12, T=1823 K.txt", delimiter=',')
#     Teide = np.loadtxt(f"lava_lib/Teide, T=1187 K.txt", delimiter = ',')
#     Forsterite = np.loadtxt('lava_lib/Forsterite, T=1473 K.txt', delimiter= ',')
#     Hawaiian = np.loadtxt('lava_lib/Hawaiian basalt, T=1573 K.txt', delimiter= ',')
#     Bglass = np.loadtxt(f"lava_lib/B-glass, T=1654 K.txt", delimiter =',')


#     # 数据提取，第一列是波长，第二列是波长对应的反射率 
#     wave_list1 =  CaFeSiOF4[:,0]
#     IS1 = CaFeSiOF4[:,1]

#     wave_list2 = Bglass[:,0] 
#     IS2 = Bglass[:,1] 

#     # 合并数据  
#     wave_combined = np.concatenate((wave_list1, wave_list2))  
#     IS_combined = np.concatenate((IS1, IS2))

#     wave_combined, IS_combined = remove_duplicate(wave_combined, IS_combined)

#     # 创建插值函数  
#     interp_func = interp1d(wave_combined, IS_combined, kind='cubic')  
    
#     wave_new = np.linspace(np.max(wave_list1), np.min(wave_list2), 100)  

#     # 计算插值光谱数据  
#     IS_new = interp_func(wave_new)
#     # print(wave_new)
#     # print(IS_new)  

#     # 可视化结果  
#     plt.plot(wave_combined, IS_combined, '-', label='Original Data')  
#     plt.plot(wave_new, IS_new, '-', label='Interpolated Data')  
#     plt.xlabel('Wavelength (nm)')  
#     plt.ylabel('Intensity')
#     plt.legend()  
#     plt.show()

#     plt.savefig('lava_lib/spectrum.png') 


class lava_Albedo:
    def __init__(self, type = 'low'):
        self.type = type
        if type =='low':
            CaFeSiOF4 = np.loadtxt(f"lava_lib/CaFeSiO F4, T=1673 K.txt", delimiter=',')
            Bglass = np.loadtxt(f"lava_lib/B-glass, T=1654 K.txt", delimiter=',')

            # 数据提取，第一列是波长，第二列是波长对应的反射率 
            wave_list1 =  CaFeSiOF4[:,0]
            IS1 = CaFeSiOF4[:,1]

            wave_list2 = Bglass[:,0] 
            IS2 = Bglass[:,1] 

            # 合并数据  
            wave_combined = np.concatenate((wave_list1, wave_list2))  
            IS_combined = np.concatenate((IS1, IS2))

            wave_combined, IS_combined = remove_duplicate(wave_combined, IS_combined)

            # 创建插值函数  
            self.interp_func = interp1d(wave_combined, IS_combined, kind='slinear')  
            
        elif type == 'high':
            CaFeSiOF6 = np.loadtxt(f"lava_lib/CaFeSiO F6, T=1673 K.txt", delimiter=',')
            Teide = np.loadtxt(f"lava_lib/Teide, T=1187 K.txt", delimiter = ',')
            Forsterite = np.loadtxt('lava_lib/Forsterite, T=1473 K.txt', delimiter= ',')
            Hawaiian = np.loadtxt('lava_lib/Hawaiian basalt, T=1573 K.txt', delimiter= ',')
            
            wave_list1 =  CaFeSiOF6[:,0]
            IS1 = CaFeSiOF6[:,1]
            
            wave_list2 = Teide[:,0]
            IS2 = Teide[:,1]
            IS2 = IS2[wave_list2 < 2.709]
            wave_list2 = wave_list2[wave_list2 < 2.709]
            
            wave_list3 = Forsterite[:,0]
            IS3 = Forsterite[:,1]
            IS3 = IS3[wave_list3 < 8.55888]
            wave_list3 = wave_list3[wave_list3 < 8.55888]
            
            wave_list4 = Hawaiian[:,0]
            IS4 = Hawaiian[:,1]
            IS4 = IS4[wave_list4 > 8.7017]
            wave_list4 = wave_list4[wave_list4 > 8.7017]
            
            wave_combined = np.concatenate((wave_list1, wave_list2, wave_list3, wave_list4))
            IS_combined = np.concatenate((IS1, IS2, IS3, IS4))
            data_combined = np.array([wave_combined,IS_combined])
            data_combined = data_combined[:, data_combined[0,:].argsort()]
            
            wave_combined = data_combined[0,:]
            IS_combined = data_combined[1,:]
            # print(data_combined)
            # print(wave_combined)
            
            wave_combined, IS_combined = remove_duplicate(wave_combined, IS_combined)
            
            # 创建插值函数
            self.interp_func = interp1d(wave_combined, IS_combined, kind='slinear')
            
            
        self.Wmax = np.max(wave_combined)
        self.Wmin = np.min(wave_combined)
        # print(self.Wmin, self.Wmax)
            
            
    # 计算插值光谱数据  
    def A_interp(self, lam):
        return self.interp_func(lam)


LA = lava_Albedo('high')

