import matplotlib.pyplot  as plt
import numpy as np

name =  'My 0.006 AU'
Tmap = np.load(f'temp/{name}/plots/Tmap0.npy')

SIZE = [180, 360]  # Size of the meshgrid
# Create meshgrid for the planet
phiP_list = np.linspace(-np.pi / 2, np.pi / 2, SIZE[0])
thetaP_list = np.linspace(0, 2 * np.pi, SIZE[1])

x, y = np.meshgrid(thetaP_list, phiP_list)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(x, y, Tmap, cmap='rainbow',linewidth=0, antialiased=False)
plt.show()
plt.savefig(f'temp/{name}/plots/Tmap_2D.png')
plt.close()
