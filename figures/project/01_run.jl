using JLD2
using Seis4CCS.RockPhysics 
using Seis4CCS.FlowSimulation

#%%
JLD2.@load "./Compass2km.jld2";
K = VtoK(v, d);
ϕ = Ktoϕ(K, d);

#%% Hyperparameter for flow simulation
hy = (n[1]-1)*d[1]        # width of the cell in y direction (now it's 2D code but CO2 lives in 3D space!)
nt = 100         # number of time steps
dt = 20         # time interval between 2 adjacent time steps (in day), you should NOT set this large otherwise the solver will be unstable
grid_ = comp_grid(n, d, hy, nt, dt);

#%% injection point
qw = zeros(nt, n[1], n[2]);
inject_x = 200
inject_z = 250
inj_loc = (inject_x*d[1], inject_z*d[2]);    # injection location (at (100,300) cell in this case)
qw[:,Int(round(inj_loc[1]/d[1])), Int(round(inj_loc[2]/d[2]))] .= 0.3;   # in [m^3/s]

qo = zeros(nt, n[1], n[2]);
prod_loc = (300*d[1], 300*d[2]);  # injection location (at (300,300) cell in this case)
qo[:,Int(round(prod_loc[1]/d[1])),Int(round(prod_loc[2]/d[2]))] .= -0.3; # also in [m^3/s]

#%% Forward
@time S, p = flow(K, ϕ, qw, qo, grid_)
Volume_inject = 0.3 * nt * dt * 24 * 60 * 60;
Volume_inow = sum(S[end,:,:].*ϕ*d[1]*d[2]*hy);

JLD2.@save "./01_result/Output.jld2" S p Volume_inject Volume_inow K ϕ nt dt inject_x inject_z