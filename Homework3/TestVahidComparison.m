%% Input data
L=1.4;      % [m]
M=1;        % [kg]
m=3/4;      % [kg]
b=0.1;      % [sec/m]
k=0.15;     % [N/m]
g=9.81;     % [m/s^2]
h=0.001;     % Step size
T_max=100;  % Simuation time
cd Local_controllers
load('K_A2.mat')
cd ..
Ts=0.001;
h=Ts;
%% Initial values
x_0 = 0;
x_dot0 = 0;
x_dotdot0 = 0;
theta_0 = 0.1;%deg2rad(10);
theta_dot0 =0 ;
theta_dotdot0 = 0;

%% State space
A = [0 0 1 0 0 0 ; 
   0 0 0 1 0 0 ;
   -k/M -m*g/M -b/M 0 0 0 ;
   k/(M*L) (m+M)/(M*L)*g b/(M*L) 0 0 0;
   0 0 0 0 -0.8 0;
   0 0 0 0 0 -0.8];
B=[0 0 1/M -1/(L*M) 0 0]';
H=[1 0 0 0 0 0;
   0 1 0 0 0 0];
E=[0 0;
    0 0 ;
    -1/M 0;
    1/(L*M) -1/(m*L);
    0.64 0;
    0 0.64];

% x, theta, x_dot, theta_dot
A_kalman= [0 0 1 0 0 0 ; 
   0 0 0 1 0 0 ;
   -k/M -m*g/M -b/M 0 -1 / M 0 ;
   k/(M*L) (m+M)/(M*L)*g b/(M*L) 0 1 / (L*M) -1 / (L*m);
   0 0 0 0 -0.8 0;
   0 0 0 0 0 -0.8];
E_kalman=[0 0;
    0 0 ;
    0 0;
    0 0;
    0.64 0;
    0 0.64];

