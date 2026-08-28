# chi2 to P-value
import numpy as np
from scipy import stats

def chi2_P(chi2, df):
    # degrees of freedom, denoted df in the implementation
    # describes the significance level
    return 1 - stats.chi2.cdf(chi2, df)

def chi2_bias(chi2, df):
    # use chi2 to calculate bias (the number of standard deviations: result = n*sigma)
    # degrees of freedom, denoted df in the implementation
    return np.sqrt(chi2)

def P_chi2(p, df):
    # df: degrees of freedom
    # p: P-value, use p to calculate chi2
    return stats.chi2.ppf(1-p, df)

def P_bias(p, df):  
    """  
    计算偏差与标准差的倍数 n*sigma  
    :param p: 置信度 P  
    :param df: 自由度  
    :return: 偏差与标准差的倍数 n*sigma  
    """  
    # # 计算对应的chi²值  
    # chi2_value = stats.chi2.ppf(1-p, df)  

    # # 计算偏差  
    # n = np.sqrt(chi2_value)  
    n = stats.norm.ppf((1 + p) / 2)  
    return n

def bias_chi2(n, df):
    """  
    get chi2 value from bias  
    :param n: 偏差与标准差的倍数 n*sigma  
    :param df: 自由度  
    :return: 对应的chi2值    
    """  
    return n**2

def bias_P(n, df):  
    """  
    计算模型的置信度 P  
    :param n: 偏差与标准差的倍数 n*sigma  
    :param df: 自由度  
    :return: 置信度 P  
    """  
    # 计算对应的chi²值  
    chi2_value = n**2  

    # 计算置信度  
    p_value = 1 - stats.chi2.cdf(chi2_value, df)  
    return p_value 


if __name__ == '__main__':
    # chi2 = np.array([34.82 , 32.19 , 43.08 , 43.24, 29.27])
    # Ndf = 11
    # P = chi2_P(chi2, Ndf)
    # print(P)

    # chi2 = np.array([1.719, 42.02, 12.20, 2.850, 1.559])
    # Ndf = 6
    # P = chi2_P(chi2, Ndf)
    # print(P)
    
    # chi2 = np.array([36.5389,74.2065,   55.2785,   46.0887,   30.8288])
    # chi2 = np.array([3.50473261, 84.07174115,  3.19704373, 24.45629201,  5.77153684 ])
    # chi2_m = np.abs(np.array([38.52280215, 33.42593126, 20.87121347, 64.04011714]) - 38.52280215)
    # chi2 = np.abs(np.array([0.43397776, 82.2543683,   3.97777042, 12.48808291]) - 0.43397776)
    
    # chi2 = np.array([2.63 ,	108.92  ,	 25.50   ,	13.73  ])
    # Ndf = 6 
    # Pm = chi2_P(chi2_m, 10)
    # P = chi2_P(chi2, 6)
    
    # Bm = P_bias(1- Pm, 10)
    # B = P_bias(1 - P, 6)
    # # print('chi2: ',chi2)
    # # # print('Reduced chi2: ',chi2/Ndf)
    # print('Pm: ',Pm)
    # print('P: ',P)
    # print('Bias of MIRI: ',Bm)
    # print("Bias of NIRCam: ",B)
    
    # P = np.array([0.68, 0.954])
    # Ndf = 10
    # chi2 = P_chi2(P, Ndf)
    # B2 = chi2_bias(chi2, Ndf)
    # B = P_bias(P, Ndf)
    # print('P: ',P)
    # print('Bias: ',B)
    # print('Bias2: ',B2)
    # chi2 = [0.16489847, 15.57786862,  3.72270754, 1.72995148]
    # P = chi2_P(chi2, 1)
    # print('significance level: ',P)
    # print('confidence level: ',1-P)
    
    Pc = np.array([  0.99945618    ,  0.99829719   ,   0.99999901   ,   0.99749943])
    B = P_bias( Pc, 9)
    print('Bias: ',B)