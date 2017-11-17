%% Input data
L=1.4;      % [m]
M=1;        % [kg]
m=3/4;      % [kg]
b=0.1;      % [sec/m]
k=0.15;     % [N/m]
g=9.81;     % [m/s^2]
h=0.001;     % Step size
T_max=300;  % Simuation time

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
%% LQR stuff
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
Q = diag([1, 1, 1, 1,0, 0]);
R = 1;

K = lqr(A, B, Q, R);

%% Kalman stuff
[Phi,delta]=c2d(A_kalman,B,h); 
[phi,gamma]=c2d(A_kalman,E_kalman,h);
X_0=[x_0 theta_0 x_dot0 theta_dot0 0 0]';

y_0=[0 0];
P_0=diag([1,1,1,1,100,100])*0.001;

Q_kalman=diag([(10), (10e-4)]);              % Process white noise

R_kalman=diag([(10e-7), (10e-7)]);           % Measurement noise
%R_kalman=diag([(0), (0)]);
%% Kalmanfilter 2
Ts=h;
D=0;
B_kalman=[B,E_kalman];
sys=ss(A_kalman,B_kalman,H,D);

[kest,L_kalman,P_kalman,M_kalman,Z_kalman] = kalmd(sys,Q_kalman,R_kalman,Ts);
A_disk=kest.a;
B_disk=kest.b;
C_disk=kest.c;
C_disk=C_disk(1:2,:);
D_disk=kest.d(3:8,:);
[A_ss, B_ss, C_ss, D_ss, T_ss]=ssdata(kest(1:2,:));