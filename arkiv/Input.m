%% Input data
L=1.4;      % [m]
M=1;        % [kg]
m=3/4;      % [kg]
b=0.1;      % [sec/m]
k=0.15;     % [N/m]
g=9.81;     % [m/s^2]
h=0.01;     % Step size
T_max=100;  % Simuation time
%% Initial values
x_0=0;
x_dot0=0;
x_dotdot0=0;
theta_0=0;
theta_dot0=0;
theta_dotdot0=0;

%% Input to Kalman
A=[0 0 1 0; 
   0 0 0 1;
   -k/M -m*g/M -b/M 0;
   k/(M*L) (m+M)/(M*L)*g b/(M*L) 0];
B=[0 0 1/M -1/(L*M)]';
H=[1 0 0 0;
   0 1 0 0];
E=[0 0;
    0 0 ;
    -1/M 0;
    1/(L*M) -1/(m*L)];

[Phi,delta]=c2d(A,B,h); 
[phi,gamma]=c2d(A,E,h);
X_0=[0 0 0 0]';
y_0=[0 0];
P_0=eye(4)*10^-6;
Q=eye(2);           % process white noise
R=eye(2)*10^-6;
