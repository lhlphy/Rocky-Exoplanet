import os
from PIL import Image, ImageDraw, ImageFont, ImageEnhance
import imageio
import numpy as np
import lib.main_function as mf
from parameter_list import *

# 图片文件夹路径
folder_path = r'./temp/6/plots'  # 修改为你的图片文件夹路径
output_path = os.path.join(folder_path, 'output.gif')

# 动画参数
duration = 0.3  # 每帧展示时间 (秒)

# 加载图片
img_files = [os.path.join(folder_path, f"plot_0_{i}.png") for i in range(0,360, 10)]

frames = []
IT = mf.blackbody_radiation(Temperature, Wavelength) * np.pi *R2**2 /(4*np.pi*a**2) *np.pi *R1**2
id = 5
Star_flux = np.load(f'./temp/R{id}/variables/star_flux.npy')
Intensity = np.load(f'./temp/R{id}/variables/TOT_Intensity.npy')
NI = Intensity/IT

for i, img_file in enumerate(img_files, start=0):
    # 打开图片
    img = Image.open(img_file)

    # 创建可编辑对象
    draw = ImageDraw.Draw(img)

    # 设置字体
    try:
        font = ImageFont.truetype("arial.ttf", 36)
    except IOError:
        font = ImageFont.load_default()

    # 计算文本尺寸（改为用textbbox来获取）
    text = f'Orbit angle: {i*10}'
    bbox = draw.textbbox((0, 0), text, font=font)
    textwidth = bbox[2] - bbox[0]
    textheight = bbox[3] - bbox[1]

    # 指定文本位置
    width, height = img.size
    x, y = width - textwidth - 10, height - textheight - 10

    # 绘制文本
    draw.text((x, y), text, font=font, fill=(0,0,0))

    enhancer = ImageEnhance.Brightness(img)
    img = enhancer.enhance(0.8 + 10*NI[i])

    # 添加带编号的帧到帧列表
    frames.append(img)

# 保存为GIF
imageio.mimsave(output_path, frames, 'GIF', duration=duration)

print(f'动画保存为 {output_path}')