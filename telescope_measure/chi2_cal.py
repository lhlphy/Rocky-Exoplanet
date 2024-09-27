import numpy as np  
from scipy.interpolate import interp1d  
import matplotlib.pyplot as plt  

jwst_data = np.loadtxt('JWST_data.txt')
jwst_ebar = np.loadtxt('JWST_errorbar.txt')
jwst_wavelength = jwst_data[:,0]
jwst_spectrum = np.array([...])  # 实测光谱值  
jwst_error = np.array([...])  # 实测误差 
 
model_wavelength = np.array([...])  # 模型波长  
model_spectrum = np.array([...])  # 模型光谱值  

 

# 对模型光谱进行插值，使其与实测数据对齐  
interp_model_spectrum = interp1d(model_wavelength, model_spectrum, kind='linear', fill_value='extrapolate')  
model_spectrum_aligned = interp_model_spectrum(jwst_wavelength)  

# 计算chi2值  
chi2 = np.sum(((jwst_spectrum - model_spectrum_aligned) ** 2) / (jwst_error ** 2))  

# 打印结果  
print(f"Chi-squared (χ²) value: {chi2}")  

# 可选：绘制数据进行可视化分析  
plt.errorbar(jwst_wavelength, jwst_spectrum, yerr=jwst_error, fmt='o', label='JWST Data')  
plt.plot(jwst_wavelength, model_spectrum_aligned, label='Model Spectrum', linestyle='--')  
plt.xlabel('Wavelength')  
plt.ylabel('Spectral Value')  
plt.legend()  
plt.title('Model vs JWST Data')  
plt.show()