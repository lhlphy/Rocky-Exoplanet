# chi2 to P-value
import numpy as np
from scipy import stats

def chi2_P(chi2, df):
    # degrees of freedom, denoted df in the implementation
    # describes the significance level
    return 1 - stats.chi2.cdf(chi2, df)


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

    # chi2 = [0.16489847, 15.57786862,  3.72270754, 1.72995148]
    # chi2 = np.array([11.95907534 , 9.84513952 , 23.94641758, 9.14017561 ])
    # dof = 5
    # P = chi2_P(chi2, dof)
    # print('significance level: ',P)
    # print('confidence level: ',1-P)
    
    # # Pc = np.array([  0.99945618    ,  0.99829719   ,   0.99999901   ,   0.99749943])
    # B = P_bias(1 - P, dof)
    # print('Bias: ',B)
    
    chi2_set = np.loadtxt('chi2_data.txt')
    chi2 = np.mean(chi2_set, axis = 0)  + np.array([11.95907534 , 9.84513952, 9.14017561, 23.94641758 ])
    print('chi2:', chi2)
    chi2 = np.abs(chi2  - chi2[0])
    print('delta chi2:', chi2)
    
    dof = 1
    P = chi2_P(chi2, dof)
    print('significance level: ',P)
    print('confidence level: ',1-P)
    
    B = P_bias(1 - P, dof)
    print('Bias: ',B)
    
    