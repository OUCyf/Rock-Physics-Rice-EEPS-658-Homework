#%%
import numpy as np
import h5py
import matplotlib.pyplot as plt
from PatchyWhite import PatchySatWhiteDutta
import matplotlib.animation as animation
# from julia.api import Julia
# jl = Julia(compiled_modules=False)
# from julia import Main
# Main.include("utilize.jl")


#%% 1. read data
fp = h5py.File("./01_result_good/Output_good.jld2",'r')
S = np.transpose(fp["S"][()], (2,0,1))  # saturation
p = np.transpose(fp["p"][()], (2,0,1))  # pressure
Volume_inject = fp["Volume_inject"][()]
Volume_inow = fp["Volume_inow"][()]
nt = fp["nt"][()]
dt = fp["dt"][()]
Kh_model= fp["K"][()]
poro_model =fp["ϕ"][()]
inject_z=fp["inject_z"][()]
inject_x=fp["inject_x"][()]
fp.close()

fp = h5py.File("./Compass2km.jld2",'r')
n = fp["n"][()]
d = fp["d"][()]
o = fp["o"][()]
v_model = np.transpose(fp["v"][()]).T
rho = np.transpose(fp["rho"][()]).T
nz=n[1]
nx=n[0]
dz=d[1]
dx=d[0]
fp.close()

#%%
fig, axs = plt.subplots(1, 1, figsize=(10, 10))
neg = axs.imshow(Kh_model, interpolation='gaussian', cmap='viridis')
fig.colorbar(neg, ax=axs, location='right', anchor=(0, 0.4), shrink=0.6)
axs.set_xlabel("X/(m)")
axs.set_ylabel("Z/(m)")
axs.set_title("Permeability Model (mD)")
axs.set_xticks(np.linspace(0,nx,num=5), (np.linspace(0,nx*dx,num=5)).astype(int))
axs.set_yticks(np.linspace(0,nz,num=5), (np.linspace(0,nz*dz,num=5)).astype(int))
fig.savefig('./01_result/Permeability_Model.pdf', dpi=800, format='pdf')

fig, axs = plt.subplots(1, 1, figsize=(10, 10))
neg = axs.imshow(poro_model, interpolation='gaussian', cmap='viridis')
fig.colorbar(neg, ax=axs, location='right', anchor=(0, 0.4), shrink=0.6)
axs.set_xlabel("X/(m)")
axs.set_ylabel("Z/(m)")
axs.set_title("Porosity Model")
axs.set_xticks(np.linspace(0,nx,num=5), (np.linspace(0,nx*dx,num=5)).astype(int))
axs.set_yticks(np.linspace(0,nz,num=5), (np.linspace(0,nz*dz,num=5)).astype(int))
fig.savefig('./01_result/Porosity_Model.pdf', dpi=800, format='pdf')

fig, axs = plt.subplots(1, 1, figsize=(10, 10))
neg = axs.imshow(v_model, interpolation='gaussian', cmap='viridis')
fig.colorbar(neg, ax=axs, location='right', anchor=(0, 0.4), shrink=0.6)
axs.set_xlabel("X/(m)")
axs.set_ylabel("Z/(m)")
axs.set_title("Velocity Model (km/s)")
axs.set_xticks(np.linspace(0,nx,num=5), (np.linspace(0,nx*dx,num=5)).astype(int))
axs.set_yticks(np.linspace(0,nz,num=5), (np.linspace(0,nz*dz,num=5)).astype(int))
fig.savefig('./01_result/Velocity_Model.pdf', dpi=800, format='pdf')



#%% 2. GIF CO2 Saturation 
fig1= plt.figure()
ims = []
for i in range(0,nt):
    neg = plt.imshow(S[i,:,:], interpolation='gaussian', cmap='viridis',vmin = -0, vmax = 1).findobj()
    # fig1.colorbar(neg, ax=axs[0], location='right', anchor=(0, 0.47), shrink=0.2).findobj()
    plt.plot(inject_x, inject_z, marker="X", markersize=10, markeredgecolor="red", markerfacecolor="red")
    plt.xlabel("X/(m)")
    plt.ylabel("Z/(m)")
    plt.title("CO2 Saturation")
    plt.xticks(np.linspace(0,nx,num=5), (np.linspace(0,nx*dx,num=5)).astype(int))
    plt.yticks(np.linspace(0,nz,num=5), (np.linspace(0,nz*dz,num=5)).astype(int))
    ims.append(neg)
ani = animation.ArtistAnimation(fig1, ims, interval=200, repeat_delay=1000)
ani.save("./01_result/S_change.gif",writer='pillow',dpi=200)


obs = np.linspace(1, stop=nt-1, num=9).astype(int)
fig1, axs = plt.subplots(3, 3, figsize=(10, 10))
for i in range(0,9):
    num_x = i//3
    num_y = i%3
    neg = axs[num_x,num_y].imshow(S[obs[i],:,:], interpolation='gaussian', cmap='viridis')
    fig1.colorbar(neg, ax=axs[num_x,num_y], location='right', anchor=(0, 0.47), shrink=0.4)
    axs[num_x,num_y].plot(inject_x, inject_z, marker="X", markersize=10, markeredgecolor="red", markerfacecolor="red")
    axs[num_x,num_y].set_xlabel("X/(m)")
    axs[num_x,num_y].set_ylabel("Z/(m)")
    axs[num_x,num_y].set_title("CO2 Saturation")
    axs[num_x,num_y].set_xticks(np.linspace(0,nx,num=5), (np.linspace(0,nx*dx,num=5)).astype(int))
    axs[num_x,num_y].set_yticks(np.linspace(0,nz,num=5), (np.linspace(0,nz*dz,num=5)).astype(int))
fig1.tight_layout()
fig1.savefig('./01_result/CO2_Saturation.pdf', dpi=800, format='pdf')

# fig2, axs = plt.subplots(3, 3, figsize=(10, 10))
# for i in range(0,9):
#     num_x = i//3
#     num_y = i%3
#     neg = axs[num_x,num_y].imshow(p[obs[i],:,:], interpolation='gaussian', cmap='viridis')
#     fig2.colorbar(neg, ax=axs[num_x,num_y], location='right', anchor=(0, 0.47), shrink=0.2)
#     axs[num_x,num_y].plot(inject_x, inject_z, marker="X", markersize=10, markeredgecolor="red", markerfacecolor="red")
#     axs[num_x,num_y].set_xlabel("X/(m)")
#     axs[num_x,num_y].set_ylabel("Z/(m)")
#     axs[num_x,num_y].set_title("CO2 Pressure")
#     axs[num_x,num_y].set_xticks(np.linspace(0,nx,num=5), (np.linspace(0,nx*dx,num=5)).astype(int))
#     axs[num_x,num_y].set_yticks(np.linspace(0,nz,num=5), (np.linspace(0,nz*dz,num=5)).astype(int))
# fig2.tight_layout()
# fig2.savefig('./vel_density.pdf', dpi=800, format='pdf')




# %% 3. compute patchy model
# v_stack, rho_stack = Main.patchy(S, v_model, density_model, poro_model, float(dx), float(dz))
# Kdry= 5.4e8   #2.735e9 # Dry bulk modulus (Pa)
# MUdry= 7e8   #3.5e9  # Dry shear modulus                   (Pa)
MUdry= 5.0e9
Kdry= 3e9

Kmin=36.6e9  # Grain bulk modulus                  (Pa)
RHOmin=2650  # Grain Bulk density                  (kg/m^3)
RHOfl2=1e3   # Density of fluid-2                   (kg/m^3)
RHOfl1=0.353e3  # Density of fluid-1                (kg/m^3)
Kfl2=2100e6   # Bulk moduli of fluid-2              (Pa)
Kfl1=0.0258e6   # Bulk moduli of fluid-1            (Pa)
VISCfl2=89e-5 # [water]  Viscosity of fluid-2          (Pa.s)
VISCfl1=1.372e-5 # [CO2]  Viscosity of fluid-1         (Pa.s)
r=0.02 # patch radius                             unit/m
r_new=0.05
fq=100 # frequencies to calculate properties       hz

Bulk=np.zeros((nt+1,nz,nx))
Vp=np.zeros((nt+1,nz,nx))
Attn=np.zeros((nt+1,nz,nx))
RHO=np.zeros((nt+1,nz,nx))
# evolution time
for tt in range(0, nt+1):
    for i in range(0,nz):
        # for j in range(0,nx):
        MUdry= rho[i,:]*1e3*(v_model[i,:]*1e3/1.73)**2
        Kdry= rho[i,:]*1e3*(v_model[i,:]*1e3)**2 - 4/3*MUdry  #2.735e9 # Dry bulk modulus (Pa)
        
        POR=poro_model[i,:]  # porosity (fraction)
        PERM=Kh_model[i,:]*9.869233e-16  #  absolute  permeability  (m^2)
        SATfl2=1-S[tt,i,:] # Saturation of fluid 1  (fraction)
        SATfl1=S[tt,i,:] # Saturation of fluid 2  (fraction)
        bulk,bp,attn,rho1 = PatchySatWhiteDutta(Kdry,MUdry,Kmin,RHOmin,
            POR,PERM,SATfl1,SATfl2,RHOfl1,RHOfl2,Kfl1,Kfl2,VISCfl1,VISCfl2,r,fq)
        Bulk[tt,i,:] = bulk
        Vp[tt,i,:] = bp
        Attn[tt,i,:] = attn
        RHO[tt,i,:] = rho1

# compute new `patch radius``
rBulk=np.zeros((nt+1,nz,nx))
rVp=np.zeros((nt+1,nz,nx))
rAttn=np.zeros((nt+1,nz,nx))
rRHO=np.zeros((nt+1,nz,nx))
for tt in range(0, nt+1):
    for i in range(0,nz):
        # for j in range(0,nx):
        MUdry= rho[i,:]*1e3*(v_model[i,:]*1e3/1.73)**2
        Kdry= rho[i,:]*1e3*(v_model[i,:]*1e3)**2 - 4/3*MUdry  #2.735e9 # Dry bulk modulus (Pa)
        
        POR=poro_model[i,:]  # porosity (fraction)
        PERM=Kh_model[i,:]*9.869233e-16  #  absolute  permeability  (m^2)
        SATfl2=1-S[tt,i,:] # Saturation of fluid 1  (fraction)
        SATfl1=S[tt,i,:] # Saturation of fluid 2  (fraction)
        bulk,bp,attn,rho1 = PatchySatWhiteDutta(Kdry,MUdry,Kmin,RHOmin,
            POR,PERM,SATfl1,SATfl2,RHOfl1,RHOfl2,Kfl1,Kfl2,VISCfl1,VISCfl2,r_new,fq)
        rBulk[tt,i,:] = bulk
        rVp[tt,i,:] = bp
        rAttn[tt,i,:] = attn
        rRHO[tt,i,:] = rho1




#%% 4. plot curve
# point_x = 250
# point_z = 230

point_x = 200  # 340  200
point_z = 220  # 320  250


fig1, axs = plt.subplots(1,2, figsize=(10, 6))
neg = axs[0].scatter(S[:,point_z,point_x],rVp[:,point_z,point_x],label='r='+str(r_new))
neg = axs[0].scatter(S[:,point_z,point_x],Vp[:,point_z,point_x],label='r='+str(r))
axs[0].legend(loc='upper right',fontsize=10, shadow=False)
axs[0].set_xlim(0, 1)
axs[0].set_xlabel('CO2 Saturation')
axs[0].set_ylabel('Vp (m/s)')
neg = axs[1].scatter(S[:,point_z,point_x], rAttn[:,point_z,point_x],label='r='+str(r_new))
neg = axs[1].scatter(S[:,point_z,point_x], Attn[:,point_z,point_x],label='r='+str(r))
axs[1].legend(loc='upper right',fontsize=10, shadow=False)
axs[1].set_xlim(0, 1)
axs[1].set_xlabel('CO2 Saturation')
axs[1].set_ylabel(r'$Q^{-1}$')
fig1.savefig('./01_result/V_curve.pdf', dpi=800, format='pdf')





#%% 5. plot velocity change and Q
tt=nt
fig1, axs = plt.subplots(1,2, figsize=(10, 10))

neg = axs[0].imshow(v_model*1e3-Vp[tt,:,:], interpolation='gaussian', cmap='viridis')
fig1.colorbar(neg, ax=axs[0], location='right', anchor=(0, 0.47), shrink=0.2)
axs[0].plot(inject_x, inject_z, marker="X", markersize=10, markeredgecolor="red", markerfacecolor="red")
axs[0].set_xlabel("X/(m)")
axs[0].set_ylabel("Z/(m)")
axs[0].set_title("Velocity Change(m/s)")
axs[0].set_xticks(np.linspace(0,nx,num=5), (np.linspace(0,nx*dx,num=5)).astype(int))
axs[0].set_yticks(np.linspace(0,nz,num=5), (np.linspace(0,nz*dz,num=5)).astype(int))

neg = axs[1].imshow(Attn[tt,:,:], interpolation='gaussian', cmap='viridis')
fig1.colorbar(neg, ax=axs[1], location='right', anchor=(0, 0.47), shrink=0.2)
axs[1].set_xlabel("X/(m)")
axs[1].set_ylabel("Z/(m)")
axs[1].set_title(r'$Q^{-1}$')
axs[1].set_xticks(np.linspace(0,nx,num=5), (np.linspace(0,nx*dx,num=5)).astype(int))
axs[1].set_yticks(np.linspace(0,nz,num=5), (np.linspace(0,nz*dz,num=5)).astype(int))
fig1.savefig('./01_result/V_change2Q.pdf', dpi=800, format='pdf')



#%% 6. gif Velocity Change
fig1= plt.figure()
ims = []
for i in range(0,nt+1):
    neg = plt.imshow(v_model*1e3-Vp[i,:,:], interpolation='gaussian', cmap='viridis').findobj()
    # fig1.colorbar(neg, ax=axs[0], location='right', anchor=(0, 0.47), shrink=0.2)
    plt.plot(inject_x, inject_z, marker="X", markersize=10, markeredgecolor="red", markerfacecolor="red")
    plt.xlabel("X/(m)")
    plt.ylabel("Z/(m)")
    plt.title("Velocity Change(m/s)")
    plt.xticks(np.linspace(0,nx,num=5), (np.linspace(0,nx*dx,num=5)).astype(int))
    plt.yticks(np.linspace(0,nz,num=5), (np.linspace(0,nz*dz,num=5)).astype(int))
    ims.append(neg)

ani = animation.ArtistAnimation(fig1, ims, interval=200, repeat_delay=1000)
ani.save("./01_result/v_change.gif",writer='pillow',dpi=200)




