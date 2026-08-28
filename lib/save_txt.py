import numpy as np  

def save_txt(data, header = 'Wavelength(micron)\tFp/F*(ppm)\terror bar(ppm)', filename = 'JWST_MIRI_data.txt'):
    # # 假设你的数据数组如下  
    # data = np.array([  
    #     [1.0, 100.0, 5.0],  
    #     [2.0, 200.0, 10.0],  
    #     [3.0, 300.0, 15.0]  
    # ])  
 
    # 定义列名  
    # header = "Wavelength(micron)\tFp/F*(ppm)\terror bar(ppm)"  

    # 保存为TXT文件  
    np.savetxt(filename, data, header=header, comments='', delimiter='\t', fmt='%.2f')  