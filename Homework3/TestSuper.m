%% Input data
L= 1.025;      % [m]
M=1;        % [kg]
m=3/4;      % [kg]
b=0.1;      % [sec/m]
k=0.15;     % [N/m]
g=9.81;     % [m/s^2]
h=0.001;     % Step size
T_max=100;  % Simuation time

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

%% Secondary system
L_2 = 1.525;

A_2 = [0 0 1 0 0 0 ;
    0 0 0 1 0 0 ;
    -k/M -m*g/M -b/M 0 0 0 ;
    k/(M*L_2) (m+M)/(M*L_2)*g b/(M*L_2) 0 0 0;
    0 0 0 0 -0.8 0;
    0 0 0 0 0 -0.8];

E_2 = [0 0;
    0 0 ;
    -1/M 0;
    1/(L_2*M) -1/(m*L_2);
    0.64 0;
    0 0.64];

%% LQR stuff
Q = diag([1, 1, 1, 1,0, 0]);
R = 1;

K = lqr(A, B, Q, R);

%% Kalman stuff
[Phi_lin,delta_lin]=c2d(A,B,h);
[phi_lin,gamma_lin]=c2d(A,E,h);
X_0=[x_0 theta_0 x_dot0 theta_dot0 0 0]';
y_0=[0 0];
P_0=diag([1,1,1,1,100,100])*0.001;

Q_kalman=diag([(10), (10e-4)]);              % Process white noise

R_kalman=diag([(10e-7), (10e-7)]);           % Measurement noise
%R_kalman=diag([(0), (0)]);


%% Supervisor
L_sup = [1.025 1.275 1.525 1.775];
A_sup = cell(1,length(L_sup));
B_sup = cell(1,length(L_sup));
E_sup = cell(1,length(L_sup));
Phi = cell(1,length(L_sup));
phi = cell(1,length(L_sup));
delta = cell(1,length(L_sup));
gamma = cell(1,length(L_sup));

for i = 1:length(L_sup)
    A_sup{i} = [0 0 1 0 0 0 ;
        0 0 0 1 0 0 ;
        -k/M -m*g/M -b/M 0 0 0 ;
        k/(M*L_sup(i)) (m+M)/(M*L_sup(i))*g b/(M*L_sup(i)) 0 0 0;
        0 0 0 0 -0.8 0;
        0 0 0 0 0 -0.8];
    
    B_sup{i}=[0 0 1/M -1/(L_sup(i)*M) 0 0]';
    
    E_sup{i}=[0 0;
        0 0 ;
        -1/M 0;
        1/(L_sup(i)*M) -1/(m*L_sup(i));
        0.64 0;
        0 0.64];
    
    [Phi{i}, delta{i}] = c2d(A_sup{i} ,B_sup{i}, h);
    [phi{i}, gamma{i}] = c2d(A_sup{i}, E_sup{i}, h);
end

C_sup = [1 0 0 0 0 0;
         0 1 0 0 0 0];
R_sup = R_kalman;

load K_A1
load K_A2
load K_A3
load K_A4

sys1 = c2d(ss(A_kA1, B_kA1, C_kA1, D_kA1), h);
sys2 = c2d(ss(A_kA2, B_kA2, C_kA2, D_kA2), h);
sys3 = c2d(ss(A_kA3, B_kA3, C_kA3, D_kA3), h);
sys4 = c2d(ss(A_kA4, B_kA4, C_kA4, D_kA4), h);

