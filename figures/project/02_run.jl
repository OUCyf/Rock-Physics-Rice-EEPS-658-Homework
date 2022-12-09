using NPZ, JLD2
using Seis4CCS.RockPhysics 
using Seis4CCS.FlowSimulation

v = npzread("./fwiflow_jobsfile/v_model.npy");
rho = npzread("./fwiflow_jobsfile/density_model.npy");
rho = rho/1f3
K = npzread("./fwiflow_jobsfile/Kh_model.npy");
ϕ = npzread("./fwiflow_jobsfile/poro_model.npy");
paras = npzread("./fwiflow_jobsfile/paras.npy");

nz=paras[1]
nx=paras[2]
dz=paras[3]
dx=paras[4]
inject_z=paras[5]
inject_x=paras[6]
product_z=paras[7]
product_x=paras[8]
nt =paras[9]        # number of time steps
dt =paras[10]    # time interval between 2 adjacent time steps (in day), you should NOT set this large otherwise the solver will be unstable


#%% Hyperparameter for flow simulation
hy = (nx-1)*dx        # width of the cell in y direction (now it's 2D code but CO2 lives in 3D space!)

grid_ = comp_grid((nx, nz), (dx, dz), hy, nt, dt);


#%% injection point
qw = zeros(nt, nx, nz);
qw[:,Int(round(inject_x)), Int(round(inject_z))] .= 0.3;   # in [m^3/s]

qo = zeros(nt, nx, nz);
qo[:,Int(round(product_x)),Int(round(product_z))] .= -0.3; # also in [m^3/s]


#%% Forward
@time S, p = flow(K, ϕ, qw, qo, grid_);

Volume_inject = 0.3 * nt * dt * 24 * 60 * 60;
Volume_inow = sum(S[end,:,:].*ϕ*dx*dz*hy);

JLD2.@save "./02_result/Output.jld2" S p Volume_inject Volume_inow




#%%
# using PyPlot
# using SlimPlotting
# figure(figsize=(20,12));
# obs = Int.(round.(range(1, stop=nt+1, length=9)));   # 9 observed time samples
# for i = 1:9
#     subplot(3,3,i)
#     plot_velocity(S[obs[i],:,:]', (dx,dz); name="", vmax=1, new_fig=false);
#     scatter(inject_x*dx, inject_z*dz, marker=".", s=100, label="injection")
# end
# suptitle("CO2 saturation [0-1]");

# figure(figsize=(20,12));
# obs = Int.(round.(range(1, stop=9+1, length=9)));   # 9 observed time samples
# for i = 1:9
#     subplot(3,3,i)
#     plot_velocity(S[obs[i],:,:]', (dx,dz); name="", vmax=1, new_fig=false);
#     scatter(inject_x*dx, inject_z*dz, marker=".", s=100, label="injection")
# end
# suptitle("CO2 saturation [0-1]");