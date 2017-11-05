%model parameters
M = 1; 
m = 3/4;
b = 0.1; 
k = 0.15;
g = 9.81; 
L = 0.9; %uncertain parameter range [0.9-1.9]


%% Initial values
x_0=0;
x_dot0=0;
x_dotdot0=0;
theta_0=pi/2;
theta_dot0=0;
theta_dotdot0=0;
MODEL;