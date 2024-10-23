import numpy as np
from lib.main_function import *
from lib.parameter_list import  PPs
import matplotlib.pyplot as plt
# 理论推导模型和我的模型的交叉验证
ID = 4
vertify_radiation(np.array([4000,5000])* 1e-9, Temperature = PPs.Stellar_T , id = ID, Ntheta = 120, NWavelength = 2)

# plot the results
DIF = np.load(f'temp/V{ID}/variables/DIF.npy')
RAT = np.load(f'temp/V{ID}/variables/RAT.npy')
Theta = np.load(f'temp/V{ID}/variables/Theta.npy')
Wave = np.load(f'temp/V{ID}/variables/Wave.npy')
DIF = DIF[:,-1]
RAT = RAT[:,-1]

fname = 'R3'
I_diffuse = np.load(f'temp/{fname}/variables/I_diffuse.npy')
Thermal = np.load(f'temp/{fname}/variables/Thermal.npy')
Wave_list = np.load(f'temp/{fname}/variables/wave_list.npy')
Theta_list = np.load(f'temp/{fname}/variables/Theta.npy')
Thermal = Thermal[-1,:]
I_diffuse = I_diffuse[-1,:]

def template_plot(xname,x1 = [], y1=[], y1name='', y1label ='', x2=[], y2=[], y2name ='', y2label ='',x3=[], y3=[], y3name ='',y3label ='',\
                  X1=[] ,d1=[], d1name='',d1label ='',X2=[], d2=[], d2name ='',d2label ='',X3=[], d3=[], d3name ='',d3label ='', id = 0, pic_name='template'):

    mpl.rcParams['axes.unicode_minus']=False       #显示负号

    a = 0.45
    fig, ax = plt.subplots(figsize=(26*a,13*a))
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.3)
    plt.rcParams['ytick.direction'] = 'in'# 刻度线显示在内部
    plt.rcParams['xtick.direction'] = 'in'# 刻度线显示在内部

    if np.size(y1) != 0:
        axpos = [0.1, 0.15, 0.7, 0.7]
        bwith = 2
        ax.spines['bottom'].set_linewidth(bwith)
        ax.spines['left'].set_linewidth(bwith)
        ax.spines['top'].set_linewidth(bwith)
        ax.spines['right'].set_linewidth(bwith)
        ax.set_position(axpos)
        #ax.axhline(y=np.average(Tc[3:]), color='gray', ls='-', )
        ax.plot(x1, y1, 'k-', label=y1label)
        if np.size(d1) != 0:
            ax.plot(X1, d1, 'ok', label=d1label)
        ax.set_ylim(ymin=0, ymax= np.max(d1)*1.1)
        ax.set_xlabel(xname, fontsize=18)
        ax.set_ylabel(y1name, fontsize=18)
        ax.tick_params(length=6, width=2)
        ax.spines['right'].set_visible(False)

    if np.size(y2) != 0:
        lambda_color = 'blue'
        labmda_ax = ax.twinx()
        labmda_ax.set_position(axpos)
        labmda_ax.plot(x2, y2, color=lambda_color, label=y2label)
        if np.size(d2) != 0:
            labmda_ax.plot(X2, d2,'ob', color=lambda_color, label=d2label)
        labmda_ax.set_ylim(ymin=0, ymax= np.max(d2)*1.1)
        labmda_ax.set_xlabel(xname, fontsize=18)
        labmda_ax.tick_params(length=6, width=2, color=lambda_color, labelcolor=lambda_color)
        labmda_ax.set_ylabel(y2name, fontsize=18, color=lambda_color)
        labmda_ax.spines['right'].set(color=lambda_color, linewidth=2.0, linestyle=':')

    if np.size(y3) != 0:
        omglog_color = 'red'
        omglog_ax = ax.twinx()
        # 使用科学计数法的刻度
        omglog_ax.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
        # 获取 y 轴 OffsetText 对象
        offset_text = omglog_ax.yaxis.get_offset_text()
        # 调整位置示例，偏移 (dx, dy) 单位是像素 (points)
        offset_text.set_position((1.12, 0))
        # 调整字体大小
        offset_text.set_size(18)  # 或者使用 offset_text.set_fontsize(12)
        omglog_ax.spines['right'].set_position(('data', np.max(x3)*1.1))
        omglog_ax.set_ylim(0, np.max(d3)*1.1)
        omglog_ax.set_position(axpos)
        omglog_ax.plot(x3, y3, '-', color=omglog_color, label=y3label)
        if np.size(d3) != 0:
            omglog_ax.plot(X3, d3, color=omglog_color, label=d3label, marker='o')
        omglog_ax.set_ylabel(y3name, fontsize=18, color=omglog_color)
        omglog_ax.tick_params(length=6, width=2, color=omglog_color, labelcolor=omglog_color)
        omglog_ax.spines['right'].set(color=omglog_color, linewidth=2.0, linestyle='-.')


    plt.title(pic_name,fontsize=18)
    fig.legend(loc='upper left', fontsize=10, bbox_to_anchor=(0.095, 0.86))
    # save the plot to temp/ folder
    os.makedirs(f'temp/V{id}/Results', exist_ok=True)
    name = f'temp/V{id}/Results/{pic_name}.png'
    plt.savefig(name)
    plt.close()


yname = 'Contrast ratio (ppm)'
xname = 'Orbital phase (rad)'

print(DIF)
print(I_diffuse)

template_plot( xname,x1 = Theta, y1 = DIF* 1e6, y1name = yname, y1label = 'Diffuse: Theoretical calculation', \
              X1 = Theta_list, d1 = I_diffuse* 1e6, d1name = yname, d1label = 'Diffuse: My model',\
                 x2 = Theta, y2 =RAT* 1e6, y2name=yname, y2label='Thermal radiation: Theoretical calculation',\
                     X2 = Theta_list, d2 = Thermal* 1e6, d2name=yname, d2label='Thermal radiation: My model', \
                        id = ID, pic_name = 'Cross-validation of the model, LHS 3844 b, a = 0.006 AU')

