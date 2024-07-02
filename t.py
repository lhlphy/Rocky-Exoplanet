import numpy as np

def rotate_vector(vector, axis, angle):
    """
    旋转向量。
    
    Parameters:
    vector (array-like): 待旋转的向量，例如 [x, y, z]
    axis (str): 'x', 'y', 或 'z' 指定绕哪个轴旋转
    angle (float): 旋转角度（以弧度为单位）
    
    Returns:
    np.ndarray: 旋转后的向量
    """
    # 定义旋转矩阵
    if axis == 'x':
        R = np.array([[1, 0, 0],
                      [0, np.cos(angle), -np.sin(angle)],
                      [0, np.sin(angle), np.cos(angle)]])
    elif axis == 'y':
        R = np.array([[np.cos(angle), 0, np.sin(angle)],
                      [0, 1, 0],
                      [-np.sin(angle), 0, np.cos(angle)]])
    elif axis == 'z':
        R = np.array([[np.cos(angle), -np.sin(angle), 0],
                      [np.sin(angle), np.cos(angle), 0],
                      [0, 0, 1]])
    else:
        raise ValueError("轴必须是 'x', 'y' 或 'z'")

    # 旋转向量
    rotated_vector = R @ np.array(vector)
    return rotated_vector

# 示例使用
if __name__ == "__main__":
    vector = [1, 0, 0]  # 待旋转的向量
    axis = 'z'  # 旋转轴
    angle = np.pi / 2  # 旋转90度
    rotated_vector = rotate_vector(vector, axis, angle)
    print("原向量:", vector)
    print("旋转后的向量:", rotated_vector)

