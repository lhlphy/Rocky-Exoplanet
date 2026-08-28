import numpy as np
from save_txt import save_txt
sim_obs = np.loadtxt("telescope_measure/sim_obs_10_16.txt", delimiter=' ')
x = sim_obs[:,0]
y = sim_obs[:,3] *1e6
Nbin = 6 # the number of bin: 6
Bin = np.zeros(Nbin)  
xloc = np.zeros(Nbin)
NGroup = sim_obs.shape[0] // Nbin # the number of pixels in each bin
for i in range(Nbin):  # calculate the errorbar of each bins
    xq = x[NGroup*i:NGroup*(i+1)]  # group pixel points to each bin
    yq = y[NGroup*i:NGroup*(i+1)]
    Bin[i] = np.sqrt(1/np.sum(1/yq**2)) # calculate the errorbar of each bins
    xloc[i] = np.mean(xq) # calculate the mean wavelength of each bin
    
data = np.zeros((Nbin, 2))
data[:,0] = xloc
data[:,1] = Bin

save_txt(data, "Wavelength(micron)\tErrorbar(ppm)", "PandExo_6bin_G395M_18hrs.txt")
print(data)
