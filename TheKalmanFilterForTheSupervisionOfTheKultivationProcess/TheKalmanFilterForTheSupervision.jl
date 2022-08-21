using Printf, LinearAlgebra, Symbolics, DifferentialEquations

include("./data/BC2_eth_pred.jl")
include("./data/BC2.jl")
include("./data/ME.jl")
include("./data/time.jl")
include("./data/timeE.jl")
# include("./data/P.jl")
using Plots

# Symbols for symbolic math calculations
@variables G E X P t
@variables Y_gx Y_ge Y_ex mu1 mu2 K_M_G K_M_E

#Variables and Parameters
initX = [2.5; 6.0; 0.2; 0.15; 0.08] # Initial State (Biomass, Glucose, Ethanol,
  # MaxGrowthRate on Glucose, Max Growth Rate on Ethanol) = (X, G, E, K_M_G, K_M_E)
initP = Diagonal([0.1,0.02,0.02,0.02,0.02]) # Initial Process Estimation / Initial Estimation Error Covariance Matrix
(initPRows, initPColumns) = size(initP)
init = [initX;reshape(initP,initPRows * initPColumns)] # Combined initial value vector

# For the odesolver
H = [0,0,1,0,0] # Observation Matrix
Q = Diagonal([0.001,0.001,0.001,0.001,0.001]) # Process Noise Covariance Matrix

R = 0.05 # Measurement noise covariance matrix
K1 = 0.1 # Monod Constant Glucose
K2 = 0.1 # Monod Constant Ethanol


#Estimated Parameter values
Ygx=0.15; # Yield glucose -> biomass
Yge=0.34; # Yield glucose -> ethanol
Yex=0.43; # Yield ethanol -> biomass

#Monod terms
mue1 = mu1*G /(G+K_M_G);
mue2 = mu2*E /(E+K_M_E) * (1 - mue1/mu1);

#Model ODE
dS = X * [ ( mue1  + mue2) ; # Biomass
 -mue1/Y_gx  ; # Glucose
 ( mue1/Y_gx*Y_ge - mue2/Y_ex) ; # Ethanol
 0;   # Maximum Specific Growth Rate on Glucose
 0;   # Maximum Specific Growth Rate on Ethanol
 ];

# Jacobian of Model with respect to state variables
F = Symbolics.jacobian(dS,[X,G,E,mu1,mu2]) # jacobian(matriz com as funções, lista com as variáveis a respeito das quais estamos tomando o jacobiano)
P = initP
dP = F * P + P * F' + Q # Defines the function dP

# Simulation / State prediction and filtering
# Replace all symbolic parameters with their respective numeric values where needed
F = substitute(F, Dict(Y_gx=>Ygx, Y_ge=>Yge, Y_ex=>Yex, K_M_G=>K1, K_M_E=>K2));
dS = substitute(dS, Dict(Y_gx=>Ygx, Y_ge=>Yge, Y_ex=>Yex, K_M_G=>K1, K_M_E=>K2));
dP = substitute(dP, Dict(Y_gx=>Ygx, Y_ge=>Yge, Y_ex=>Yex, K_M_G=>K1, K_M_E=>K2));

# function P_Update(oldP, )

println("Printing dS")
println(dS)

# Build the ODE function to feed the solver
function tkftcp!(du,u,p,t) # The Kalman Filter For the Cultivation Process's equation
  du[1] = substitute(dS[1], Dict([ X=>u[1], G=>u[2], E=>u[3], mu1=>u[4], mu2=>u[5] ]) ).val
  du[2] = substitute(dS[2], Dict([ X=>u[1], G=>u[2], E=>u[3], mu1=>u[4], mu2=>u[5] ]) ).val
  du[3] = substitute(dS[3], Dict([ X=>u[1], G=>u[2], E=>u[3], mu1=>u[4], mu2=>u[5] ]) ).val
  du[4] = substitute(dS[4], Dict([ X=>u[1], G=>u[2], E=>u[3], mu1=>u[4], mu2=>u[5] ]) ).val
  du[5] = substitute(dS[5], Dict([ X=>u[1], G=>u[2], E=>u[3], mu1=>u[4], mu2=>u[5] ]) ).val
end

tspan = (0.0,time[end])
prob = ODEProblem(tkftcp!,initX,tspan)
solutions = []

sol = solve(prob,save_everystep=true)
plot(sol)
savefig("./charts/result.png")

# Simulate the process from one ethanol gas measurement time to the next:
t0 = 0;
# MC = zeros(0,5); # store filtered states in these variables
# SimState = zeros(0,5);
# SimTime = [];

#for time in timeE
#  tspan = [t0, time]
#  
#end

# STARTING LOOP ################################################
#for i = 1:length(timeE)
# tspan = [t0 timeE(i)];
# [T,state] = ode45(OdeSys, tspan, init); #  simulate / solve model
# PS = state(end,1:5)';  #  predicted state
# MS = ME(i);   #  measured state
# P = reshape(state(end,6:end),5,5); #  process covariance matrix
# K = P*H'/(H*P*H'+R);  #  kalman gain matrix
# disp(P*H')
# FS = PS + K * (MS-PS(3));  #  filtered state
# Pfilt = P-K*H*P ;   #  filtered process covariance matrix
# init = [FS; Pfilt(:)];   #  new initial condition
# t0 = timeE(i);   #  new starting time for next iteration
#
# #  Save intermediate states for plotting
# MC  = [MC;FS'];
# state(end,1:3) = NaN;
# SimState = [SimState; state(:,1:5)];
# SimTime = [SimTime; T];
#end
println("Fin")