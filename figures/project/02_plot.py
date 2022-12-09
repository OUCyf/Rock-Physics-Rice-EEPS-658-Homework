#%%
import numpy as np
import h5py
import matplotlib.pyplot as plt
from PatchyWhite import PatchySatWhiteDutta
import matplotlib.animation as animation


#%% 1.read data
fp = h5py.File("./02_result/Output_good.jld2", 'r')
S = np.transpose(fp["S"][()], (2,0,1))  # saturation
p = np.transpose(fp["p"][()], (2,0,1))  # pressure
Volume_inject = fp["Volume_inject"][()]
Volume_inow = fp["Volume_inow"][()]
fp.close()

v_model = np.load("./fwiflow_jobsfile/v_model.npy").T
poro_model = np.load("./fwiflow_jobsfile/poro_model.npy").T
density_model = np.load("./fwiflow_jobsfile/density_model.npy").T
Kh_model = np.load("./fwiflow_jobsfile/Kh_model.npy").T
paras = np.load("./fwiflow_jobsfile/paras.npy").T
nz=paras[0]
nx=paras[1]
dz=paras[2]
dx=paras[3]
inject_z=paras[4]
inject_x=paras[5]
product_z=paras[6]
product_x=paras[7]
nt=paras[8]
dt=paras[9]


# %% 2. plot CO2 Saturation
obs = np.linspace(1, stop=nt-1, num=9).astype(int)
fig1, axs = plt.subplots(3, 3, figsize=(10, 10))
for i in range(0,9):
    num_x = i//3
    num_y = i%3
    neg = axs[num_x,num_y].imshow(S[obs[i],:,:], interpolation='gaussian', cmap='viridis')
    fig1.colorbar(neg, ax=axs[num_x,num_y], location='right', anchor=(0, 0.47), shrink=0.2)
    axs[num_x,num_y].plot(inject_x, inject_z, marker="X", markersize=10, markeredgecolor="red", markerfacecolor="red")
    axs[num_x,num_y].set_xlabel("X/(m)")
    axs[num_x,num_y].set_ylabel("Z/(m)")
    axs[num_x,num_y].set_title("CO2 Saturation")
    axs[num_x,num_y].set_xticks(np.linspace(0,nx,num=5), (np.linspace(0,nx*dx,num=5)).astype(int))
    axs[num_x,num_y].set_yticks(np.linspace(0,nz,num=5), (np.linspace(0,nz*dz,num=5)).astype(int))
fig1.tight_layout()
fig1.savefig('./02_result/CO2_Saturation.pdf', dpi=800, format='pdf')

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
ani.save("02_result/S_change.gif",writer='pillow',dpi=200)

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
Kdry=2.735e9 # Dry bulk modulus (Pa)
MUdry=3.5e9  # Dry shear modulus                   (Pa)
Kmin=36.6e9  # Grain bulk modulus                  (Pa)
RHOmin=2650  # Grain Bulk density                  (kg/m^3)
RHOfl2=1e3   # Density of fluid-2                   (kg/m^3)
RHOfl1=0.353e3  # Density of fluid-1                (kg/m^3)
Kfl2=2100e6   # Bulk moduli of fluid-2              (Pa)
Kfl1=0.0258e6   # Bulk moduli of fluid-1            (Pa)
VISCfl2=89 # [water]  Viscosity of fluid-2          (Pa.s)
VISCfl1=1.372 # [CO2]  Viscosity of fluid-1         (Pa.s)
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
        POR=poro_model[i,:]  # porosity (fraction)
        PERM=Kh_model[i,:]*9.869233e-16  #  absolute  permeability  (m^2)
        SATfl2=1-S[tt,i,:] # Saturation of fluid 1  (fraction)
        SATfl1=S[tt,i,:] # Saturation of fluid 2  (fraction)
        bulk,bp,attn,rho = PatchySatWhiteDutta(Kdry,MUdry,Kmin,RHOmin,
            POR,PERM,SATfl1,SATfl2,RHOfl1,RHOfl2,Kfl1,Kfl2,VISCfl1,VISCfl2,r,fq)
        Bulk[tt,i,:] = bulk
        Vp[tt,i,:] = bp
        Attn[tt,i,:] = attn
        RHO[tt,i,:] = rho

# compute new `patch radius``
rBulk=np.zeros((nt+1,nz,nx))
rVp=np.zeros((nt+1,nz,nx))
rAttn=np.zeros((nt+1,nz,nx))
rRHO=np.zeros((nt+1,nz,nx))
for tt in range(0, nt+1):
    for i in range(0,nz):
        # for j in range(0,nx):
        POR=poro_model[i,:]  # porosity (fraction)
        PERM=Kh_model[i,:]*9.869233e-16  #  absolute  permeability  (m^2)
        SATfl2=1-S[tt,i,:] # Saturation of fluid 1  (fraction)
        SATfl1=S[tt,i,:] # Saturation of fluid 2  (fraction)
        bulk,bp,attn,rho = PatchySatWhiteDutta(Kdry,MUdry,Kmin,RHOmin,
            POR,PERM,SATfl1,SATfl2,RHOfl1,RHOfl2,Kfl1,Kfl2,VISCfl1,VISCfl2,r_new,fq)
        rBulk[tt,i,:] = bulk
        rVp[tt,i,:] = bp
        rAttn[tt,i,:] = attn
        rRHO[tt,i,:] = rho

#%% 4. plot curve
point_x = 60
point_z = 42

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
axs[1].set_ylabel('Q')
fig1.savefig('./02_result/V_curve.pdf', dpi=800, format='pdf')




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
axs[1].set_title("Q")
axs[1].set_xticks(np.linspace(0,nx,num=5), (np.linspace(0,nx*dx,num=5)).astype(int))
axs[1].set_yticks(np.linspace(0,nz,num=5), (np.linspace(0,nz*dz,num=5)).astype(int))
fig1.savefig('./02_result/V_change2Q.pdf', dpi=800, format='pdf')



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
ani.save("./02_result/v_change.gif",writer='pillow',dpi=200)


# %%
