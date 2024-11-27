import numpy as np
import matplotlib.pyplot as plt
from exo_rank_function import Rank_exo

class Rank_plot:
    def __init__(self, lam_range, Nlam=5, T_liq = 1600, Nname = 10, Nhead = 30):
        '''    
        lam_range: range of observed wavelengths
        Nlam: number of wavelengths to be plotted, divide the range into Nlam parts
        T_liq: liquid temperature, 作为完全融化区域的温度阈值
        # Nname: number of names of planets
        Nhead: number of planets in the head'''
        self.lam_range = lam_range
        self.Nlam = Nlam
        self.T_liq = T_liq
        self.lam_list = np.linspace(lam_range[0], lam_range[1], Nlam)
        self.Nname = Nname
        self.Nhead = Nhead
        
    def Exo_rank(self, lam_spec):
        return Rank_exo(lam_spec=lam_spec, T_liq=self.T_liq, Nhead=self.Nhead, rank_standard='Specular')
    
    # def Plot_Rank(self):
    #     plt.figure()
    #     for lam_spec in self.lam_list:
    #         top_exo = self.Exo_rank(lam_spec)
    #         for i in range(self.Nhead-1, -1, -1): # 倒序绘制，避免Rank 1的点被遮挡
    #             if i == 0:
    #                 plt.plot(lam_spec, top_exo['Specular'].iloc[i], 'o', color='red')
    #                 plt.text(lam_spec, top_exo['Specular'].iloc[i], top_exo['pl_name'].iloc[i].split('.txt')[0], fontsize=8)
    #                 continue
    #             plt.plot(lam_spec, top_exo['Specular'].iloc[i], 'o', color='gray')
                
    #     plt.xlabel('Wavelength (nm)')
    #     plt.ylabel('Specular_corr')
    #     plt.show()
    def Plot_Rank(self):
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()  # Create a second y-axis

        # Calculate bar width
        if len(self.lam_list) > 1:
            bar_width = 0.8 * (self.lam_list[1] - self.lam_list[0])
        else:
            bar_width = 0.1  # Default width if only one wavelength

        for lam_spec in self.lam_list:
            top_exo = self.Exo_rank(lam_spec)

            # Plot the 'Specular' data
            for i in range(self.Nhead - 1, -1, -1):  # Reverse order plotting
                if i == 0:
                    ax1.plot(lam_spec, top_exo['Specular'].iloc[i], 'o', color='red', zorder=2)
                    ax1.text(lam_spec, top_exo['Specular'].iloc[i], top_exo['pl_name'].iloc[i].split('.txt')[0], fontsize=8)
                else:
                    ax1.plot(lam_spec, top_exo['Specular'].iloc[i], 'o', color='black', zorder=2)

            # Plot the bar chart for 'Spl_CR' of the top-ranked exoplanet
            ax2.bar(lam_spec, top_exo['Spl_CR'].iloc[0], width=bar_width, color='blue', alpha=0.4, zorder=1)

        # Set labels and styles
        # ax1.set_yscale('log')
        ax1.set_xlabel('Wavelength (nm)', fontsize=15)
        ax1.set_ylabel('Specular_corr', fontsize=15)
        ax1.tick_params(axis='both', labelsize=12)

        ax2.set_ylabel('contrast ratio of specular reflection (ppm)', color='blue', fontsize=15)
        ax2.tick_params(axis='y', labelcolor='blue', labelsize=12)
        ax2.spines['right'].set_color('blue')
        # 设置边框的宽度
        ax2.spines['right'].set_linewidth(2)
        ax1.spines['left'].set_linewidth(2)
        ax1.spines['top'].set_linewidth(2)
        ax1.spines['bottom'].set_linewidth(2)

        plt.tight_layout()  # 自动调整布局
        plt.show()
            
if __name__ == '__main__':
    rank_plot = Rank_plot([1e-6, 5e-6], Nlam=20, T_liq=1600, Nname=10, Nhead=30)
    rank_plot.Plot_Rank()
        

    
