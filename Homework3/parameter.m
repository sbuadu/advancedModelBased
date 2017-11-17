%%%PARAMETER%%%%
m=3/4;            %mass of ball
M=1;              %mass of cart
g=9.8;            %geranesh zamin
b=0.1;            %friction of cart   
%L_max=1.2;       %Minimum of the length of the pendulum  
%L_min=.6;        %Maxiimum of the length of the pendulum  
Rx_SIMULATION=10^-10;         %Power of the measurment noise
Rth_SIMULATION=10^-11;       %Power of the measurment noise
Rx=10^-6;         %Power of the measurment noise to design control
Rth=10^-7;       %Power of the measurment noise to design control
Qf=1;              %Power of the disturbance f
Qw=10;             %Power of the disturbance w
alpha_f=.8;        %bandwidth of filter which makes disturbance f (effecting on ball) from white noise
gain_f=.008;       %gain of filter which makes disturbance f (effecting on ball) from white noise
alpha_w=.8;        %bandwidth of filter which makes disturbance w (effecting on cart) from white noise
gain_w=.8;         %gain of filter which makes disturbance w (effecting on cart) from white noise
%alpha_theta=1;   %bandwidth of filter which  penalizes the output (angel of the pendulum or theta)
%gain_Atheta=1;   %gain of filter which penalizes the output (angel of the pendulum or theta)
alpha_x=.8;        %bandwidth of filter which  penalizes the output (position of the cart or x)
%gain_Ax=1;       %gain of filter which penalizes the output (position of the cart or x)
Ts=.001;
%L_nom=1.8;
N=4;              %Number of models
k=.15;            %stifess coeficient of spring which connects cart o wall

% A_L=[0               0                1          0
%      0               0                0          1
%      -k/m        -m*g/M            -b/M          0
%  k/(M*L_nom)  (m+M)*g/(M*L_nom)   b/(L_nom*M)    0];
% 
% B_L=[ 0                0           0
%       0                0           0
%      1/M             -1/M          0
%    -1/(L_nom*M)  1/(L_nom*M)  -1/(L_nom*m)  ];
% 
% C_L =[0 1 0 0
%       1 0 0 0];
% 
% D_L =[0  0  0
%       0  0  0];
% A_L=[0             0                1           0
%      0             0                0           1
%   -k/M          -m*g/M            -b/M          0
% k/(M*L_nom)  (m+M)*g/(M*L_nom)   b/(L_nom*M)    0];
% 
% B_L=[ 0
%       0
%      1/M
%    -1/(L_nom*M)];
% 
% G_L=[ 0
%       0
%     -1/M
%     1/(L_nom*M)];
% 
% 
% 
% C_L =[0 1 0 0
%       1 0 0 0];
% 
% D_L =[0
%       0];
cd Local_controllers
load('K_A2.mat')
cd ..
Ts=0.01;
h=Ts;