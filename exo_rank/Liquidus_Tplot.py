import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# 数据
# sample = [
#     'B-glass', 
#     'B-Nat', 
#     'Fasnia', 
#     'MORB', 
#     'B´arðarbunga\nbasalts',  # 添加换行符
#     'Quartz', 
#     'Hawaiian\nbasalt',        # 添加换行符
#     'Teide', 
#     'Erebus'
# ]
sample = [f'S{i}' for i in range(1, 10)]
T_liq_celsius = [1176.56, 1174.02, 1201.16, 1169.92, 1193.16, 1722.85, 1172.66, 1072.85, 1085.94]
T_liq = [temp + 273.15 for temp in T_liq_celsius]  # 转换为开尔文

# 额外数据
Planet = ['K2-141 b', 'Kepler-607 b',  'Kepler-1037 b', 'Kepler-78 b', 'Kepler-10 b', '55 Cancri e']
Teq_celsius = [2103, 1876, 1652, 2188, 2223, 1958]  # 假设这些温度以摄氏度为单位
Teq = [temp for temp in Teq_celsius]  # 转换为开尔文

# 设置图形大小
plt.figure(figsize=(7.5, 6))  # 增加宽度以适应标签

# 绘制柱状图
bars = plt.bar(sample, T_liq, color='skyblue')

# 添加标题和标签
# plt.title('Liquidus Temperature of Various Samples', fontsize=16)
plt.xlabel('Sample index', fontsize=18)
plt.ylabel('Liquidus Temperature (K)', fontsize=19)  # 修改单位为开尔文
ytick_values = [0, 500, 1000, 1500, 2000, 2234]
plt.yticks(ytick_values, fontsize=16)
for tick, value in zip(plt.gca().get_yticklabels(), ytick_values):
    if value == 1500:
        tick.set_color('blue')
    elif value == 2234:
        tick.set_color('red')
plt.ylim(0, 2500)

# 倾斜 x 轴标签
# plt.xticks(rotation=60, ha='center', fontsize=14)
plt.xticks(ha='center', fontsize=16)
# 自动调整 x 轴范围
plt.xlim(-0.5, len(sample) - 1.5)

# # 添加数据标签在每个柱子顶部
# for bar in bars:
#     height = bar.get_height()
#     plt.text(
#         bar.get_x() + bar.get_width() / 2., 
#         height + 10,
#         f'{height:.2f} K', 
#         ha='center', 
#         va='bottom', 
#         fontsize=10
#     )

# 绘制平行于x轴的虚线，并在附近标注Planet名称
ax = plt.gca()  # 获取当前轴
x_min, x_max = ax.get_xlim()  # 获取x轴范围

# 获取最后一个柱子的x坐标
last_bar = bars[-1]
last_bar_x = last_bar.get_x() + last_bar.get_width() / 2.

# 定义虚线颜色
line_color = 'gray'
label_x_positions = [
    x_min + offset * (x_max - x_min)
    for offset in [0.03, 0.03, 0.03, 0.31, 0.55, 0.31]
]
plt.xlim(x_min, x_max + 1)  # 增加x轴上限

plt.axhspan(2234, 2500, color='red', alpha=0.3, zorder=0)
plt.axhline(y=2234, color='red', linewidth=1.5)
plt.axhline(y=1500, color='blue', linewidth=1.5)
plt.text(
    (x_min + x_max + 1) / 2,
    (2234 + 2500) / 2,
    'Atmospheres become optically thick',
    color='red',
    fontsize=15,
    ha='center',
    va='center',
)

planet_labels = []
for planet, teq, label_x in zip(Planet, Teq, label_x_positions):
    # 在最后一个柱子的上方添加标签，稍微向左移动
    if planet == 'Kepler-808 b':
        VA = 'top'
    # elif planet == 'Kepler-1320 b':
    #     VA ='bottom'
    else:
        VA = 'center'
        
    label = plt.text(
        label_x,
        teq, 
        planet,
        color='gray',
        fontsize=15,
        va= VA,
        ha='left',  # 标签左对齐
        
    )
    planet_labels.append((teq, label))

plt.gcf().canvas.draw()
for teq, label in planet_labels:
    text_bbox = label.get_window_extent()
    text_x_min = ax.transData.inverted().transform((text_bbox.x0, text_bbox.y0))[0]
    text_x_max = ax.transData.inverted().transform((text_bbox.x1, text_bbox.y1))[0]
    plt.hlines(teq, x_min, text_x_min - 0.03, color=line_color, linestyle='--', linewidth=1)
    plt.hlines(teq, text_x_max - 0.11, x_max + 1, color=line_color, linestyle='--', linewidth=1)

# 调整布局以防止标签被截断
plt.tight_layout()
plt.savefig('Liquidus_Tplot.png')
plt.savefig('Liquidus_Tplot.pdf')
# 显示图形
plt.show()
