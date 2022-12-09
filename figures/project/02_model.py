#%%
import numpy as np
import matplotlib.pyplot as plt
import os

def perm2poro(K):
    """
    Convert the permeability to porosity.

    Parameters
    ----------
        K : array
            permeability.
            
    Returns
    -------
        phi : array
            porosity.        
    """    

    nz,nx = K.shape

    poro_model = np.zeros((nz,nx))
    for i in range(0, nz):
        for j in range(0, nx):
            # p = np.polynomial.Polynomial([-0.0314**2*K[i,j], 2*0.0314**2*K[i,j], -0.0314**2*K[i,j], 1.527**2])
            p = [1.527**2,-0.0314**2*K[i,j], 2*0.0314**2*K[i,j],-0.0314**2*K[i,j]]
            poro_model[i,j] = np.min(np.real(np.roots(p)[np.where(np.int64(np.real(np.roots(p)) == np.roots(p)) == 1)[0]]))

    dd = 0
    for i in range(10, nz):
        if i%3 == 0:
            dd += 1
        if (i>=38 and i<=40):
            poro_model[i, 0:70-dd] = 0.055

    # layer-4
    dd = 0
    for i in range(10, nz):
        if i%3 == 0:
            dd += 1
        if (i>=35 and i<=37):
            poro_model[i, 75-dd:nx] = 0.055

    return poro_model
    
#%% 1. Create Model
nt=100
dt=20
nz = 70
nx = 100
dx = 6
dz = 6
inject_z = 43
inject_x = 65
product_z = 50
product_x = 80
a = 1.03 # convert velocity to permeability. Kh=a*v
paras = np.array([nz,nx,dz,dx,inject_z,inject_x,product_z,product_x,nt,dt])

v_model = 2*np.ones((nz,nx))  # km/s
density_model = 2.3e3*np.ones((nz,nx))
Kh_model = a*v_model+200

# layer-1: water [0m - 50m]
v_model[0:10,:] = 1.5 # km/s
density_model[0:10,:] = 1.0e3 # kg/m3
Kh_model[0:10,:] = 1e-10

# layer-2: fault
dd = 0
for i in range(10, nz):
    if i%3 == 0:
        dd += 1
    v_model[i, 70-dd:75-dd] = (1.5e3 + i * 20)/1e3
    density_model[i, 70-dd:75-dd] = 2.1e3 + i * 5
    Kh_model[i, 70-dd:75-dd] = a*1.5e3/1e3 +300 - i * 2

# layer-3
dd = 0
for i in range(10, nz):
    if i%3 == 0:
        dd += 1
    if (i>=38 and i<=40):
        Kh_model[i, 0:70-dd] = 1e-3
    if (i > 40):
        v_model[i, 0:70-dd] = (2e3 + i * 30)/1e3
        density_model[i, 0:70-dd] = 2.5e3 + i * 5
        Kh_model[i, 0:70-dd] = a*(2e3)/1e3 +80- i * 0.5

# layer-4
dd = 0
for i in range(10, nz):
    if i%3 == 0:
        dd += 1
    if (i>=35 and i<=37):
        Kh_model[i, 75-dd:nx] = 1e-3
    if (i > 37):
        v_model[i, 75-dd:nx] = (2e3 + i * 30)/1e3
        density_model[i, 75-dd:nx] = 2.5e3 + i * 5
        Kh_model[i, 75-dd:nx] = a*(2e3)/1e3 +80-  i * 0.5

poro_model = perm2poro(Kh_model)


#%% 2. Plot Model
fig1, axs = plt.subplots(1,2, figsize=(10, 10))

neg = axs[0].imshow(v_model, interpolation='gaussian', cmap='viridis')
fig1.colorbar(neg, ax=axs[0], location='right', anchor=(0, 0.47), shrink=0.2)
axs[0].plot(inject_x, inject_z, marker="X", markersize=10, markeredgecolor="red", markerfacecolor="red")
axs[0].set_xlabel("X/(m)")
axs[0].set_ylabel("Z/(m)")
axs[0].set_title("Velocity Model (km/s)")
axs[0].set_xticks(np.linspace(0,nx,num=5), (np.linspace(0,nx*dx,num=5)).astype(int))
axs[0].set_yticks(np.linspace(0,nz,num=5), (np.linspace(0,nz*dz,num=5)).astype(int))

neg = axs[1].imshow(density_model, interpolation='gaussian', cmap='viridis')
fig1.colorbar(neg, ax=axs[1], location='right', anchor=(0, 0.47), shrink=0.2)
axs[1].set_xlabel("X/(m)")
axs[1].set_ylabel("Z/(m)")
axs[1].set_title("Density Model (kg/m3)")
axs[1].set_xticks(np.linspace(0,nx,num=5), (np.linspace(0,nx*dx,num=5)).astype(int))
axs[1].set_yticks(np.linspace(0,nz,num=5), (np.linspace(0,nz*dz,num=5)).astype(int))

fig1.savefig('./02_result/vel_density.pdf', dpi=800, format='pdf')


fig2, axs = plt.subplots(1,2, figsize=(10, 10))

neg = axs[0].imshow(Kh_model, interpolation='gaussian', cmap='viridis')
fig2.colorbar(neg, ax=axs[0], location='right', anchor=(0, 0.47), shrink=0.2)
axs[0].plot(inject_x, inject_z, marker="X", markersize=10, markeredgecolor="red", markerfacecolor="red")
axs[0].set_xlabel("X/(m)")
axs[0].set_ylabel("Z/(m)")
axs[0].set_title("Permeability Model")
axs[0].set_xticks(np.linspace(0,nx,num=5), (np.linspace(0,nx*dx,num=5)).astype(int))
axs[0].set_yticks(np.linspace(0,nz,num=5), (np.linspace(0,nz*dz,num=5)).astype(int))

neg = axs[1].imshow(poro_model, interpolation='gaussian', cmap='viridis')
fig2.colorbar(neg, ax=axs[1], location='right', anchor=(0, 0.47), shrink=0.2)
axs[1].set_xlabel("X/(m)")
axs[1].set_ylabel("Z/(m)")
axs[1].set_title("Porosity Model")
axs[1].set_xticks(np.linspace(0,nx,num=5), (np.linspace(0,nx*dx,num=5)).astype(int))
axs[1].set_yticks(np.linspace(0,nz,num=5), (np.linspace(0,nz*dz,num=5)).astype(int))

fig2.savefig('./02_result/pere_poro.pdf', dpi=800, format='pdf')


#%% 3. Save file
os.makedirs('./fwiflow_jobsfile', exist_ok=True)
np.save(file="./fwiflow_jobsfile/v_model.npy", arr=v_model.T)
np.save(file="./fwiflow_jobsfile/poro_model.npy", arr=poro_model.T)
np.save(file="./fwiflow_jobsfile/density_model.npy", arr=density_model.T)
np.save(file="./fwiflow_jobsfile/Kh_model.npy", arr=Kh_model.T)
np.save(file="./fwiflow_jobsfile/paras.npy", arr=paras)


# %%
