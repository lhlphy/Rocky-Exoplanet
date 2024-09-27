import numpy as np
import os
import sys
# 将B文件夹添加到Python路径  
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../lib'))) 

from function_library import chi2_cal

data_value = np.loadtxt(f"telescope_measure/JWST_data.txt", delimiter=',')
data_ebar = np.loadtxt(f"telescope_measure/JWST_errorbar.txt", delimiter=',')
test_data = np.loadtxt(f"telescope_measure/vertify_temp.txt", delimiter=',')
    
DN = data_value.shape[0]
data = np.zeros([DN, 3])
data[:,0:2] = data_value
for i in range(DN):
    data[i, 2] = np.abs(data_ebar[i*2,1] - data_ebar[i*2+1,1])/2
    
chi2 = chi2_cal(data[:,0], data[:,1], data[:,2], test_data[:,0], test_data[:,1])
print(chi2)
    
