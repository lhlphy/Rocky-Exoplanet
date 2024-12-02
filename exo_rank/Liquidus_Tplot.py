import numpy as np
import matplotlib.pyplot as plt

# 数据
sample = [
    'B-glass', 
    'B-Nat', 
    'Fasnia', 
    'MORB', 
    'B´arðarbunga\nbasalts',  # 添加换行符
    'Quartz', 
    'Hawaiian\nbasalt',        # 添加换行符
    'Teide', 
    'Erebus'
]
T_liq_celsius = [1176.56, 1174.02, 1201.16, 1169.92, 1193.16, 1722.85, 1172.66, 1072.85, 1085.94]
T_liq = [temp + 273.15 for temp in T_liq_celsius]  # 转换为开尔文

# 额外数据
Planet = ['K2-141 b', 'Kepler-1320 b', 'Kepler-808 b', 'Kepler-607 b',  'Kepler-1037 b']
Teq_celsius = [2103, 1409, 1324, 1876, 1652]  # 假设这些温度以摄氏度为单位
Teq = [temp + 273.15 for temp in Teq_celsius]  # 转换为开尔文

# 设置图形大小
plt.figure(figsize=(6, 6))  # 增加宽度以适应标签

# 绘制柱状图
bars = plt.bar(sample, T_liq, color='skyblue')

# 添加标题和标签
plt.title('Liquidus Temperature of Various Samples', fontsize=16)
plt.xlabel('Sample Type', fontsize=15)
plt.ylabel('Liquidus Temperature (K)', fontsize=15)  # 修改单位为开尔文
plt.yticks(fontsize=13)

# 倾斜 x 轴标签
plt.xticks(rotation=60, ha='center', fontsize=12)
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
line_color = 'blue'

for planet, teq in zip(Planet, Teq):
    # 绘制虚线
    plt.axhline(y=teq, color=line_color, linestyle='--', linewidth=1)
    
    # 在最后一个柱子的上方添加标签，稍微向左移动
    if planet == 'Kepler-808 b':
        VA = 'top'
    elif planet == 'Kepler-1320 b':
        VA ='bottom'
    else:
        VA = 'center'
        
    plt.text(
        last_bar_x + 0.03 * (x_max - x_min),  # 将x位置向左移动10%
        teq, 
        planet,
        color=line_color,
        fontsize=13,
        va= VA,
        ha='right',  # 标签右对齐
        backgroundcolor='white',  # 背景为白色，避免覆盖图像
        bbox=dict(facecolor='white', edgecolor='none', pad=1, alpha = 1)
        
    )

# 调整x轴范围以给标签留出空间
plt.xlim(x_min, x_max + 1)  # 增加x轴上限

# 调整布局以防止标签被截断
plt.tight_layout()
plt.savefig('lava_Tplot.png')
plt.savefig('lava_Tplot.pdf')
# 显示图形
plt.show()
