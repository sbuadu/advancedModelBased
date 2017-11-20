clear
close all
%% Input data
L = 1.4;      % [m]
M = 1;        % [kg]
m = 3/4;      % [kg]
b = 0.1;      % [sec/m]
k = 0.15;     % [N/m]
g = 9.81;     % [m/s^2]
h = 0.001;     % Step size
T_max = 200;  % Simuation time


%% Input to Kalman
A = [0 0 1 0;
    0 0 0 1;
    -k/M -m*g/M -b/M 0;
    k/(M*L) (m+M)/(M*L)*g b/(M*L) 0];
B = [0 0 1/M -1/(L*M)]';
H = [1 0 0 0;
    0 1 0 0];
E = [0 0;
    0 0;
    -1/M 0;
    1/(L*M) -1/(m*L)];

%% Initial values
x_0 = 0;
x_dot0 = 0;
x_dotdot0 = 0;
theta_0 = 0;
theta_dot0 = 0;
theta_dotdot0 = 0;

Q = diag([1,1,1,1]);
R = 1;


K = lqr(A, B, Q, R ); %controller optimized for L=1.4



% %% Checking and plotting control action for cheap and expensive control
% 
% % cheap
% Q = diag([1,1,1,1]);
% R = 0.1;
% K = lqr(A, B, Q, R ); 
% sim('Model.slx');
% load('expensive.mat');
% u_cheap = data(2,:);
% 
% %expensive
% Q = diag([1,1,1,1]);
% R = 100;
% K = lqr(A, B, Q, R ); 
% sim('Model.slx');
% load('expensive.mat');
% u_exp=data(2,:);
% 
% figure
% plot(data(1,:), u_cheap);
% hold on
% plot(data(1,:), u_exp);
% grid on;
% legend({'Cheap control (Q/R= 10)','Expensive control (Q/R=0.01)'}, 'FontSize', 14)
% title('Cheap vs expensive control: Control input u', 'Interpreter', 'latex',  'FontSize', 18);
% xlabel('Time [s]','Interpreter', 'latex',  'FontSize', 18);
% ylabel('Force [N]', 'Interpreter', 'latex',  'FontSize', 18);


% %% checking controller gains for cheap and expensive control
% a = linspace(0,1000,500);
% k = [];
% ratio_1 = a./R;
% ratio_2 = a./R;
% ratio_3 = a./R;
% ratio_4 = a./R;
% for i = 1 :length(a)
% Q = diag([a(i), a(i), a(i), a(i)]);
% R = 1;
%
% K = lqr(A, B, Q, R );
% k(i,:) = K;
% end


% %% Plotting the controller gains.
% figure
% plot(ratio_1,k(:,1));
% title('$x$ controller gain', 'Interpreter', 'latex', 'FontSize', 18);
% grid('on');
% xlabel('$Q_x/R_u $', 'Interpreter', 'latex', 'FontSize', 18);
% ylabel('$K_x$', 'Interpreter', 'latex', 'FontSize', 18);
% set(get(gca,'ylabel'),'rotation',0)
%
% figure
% plot(ratio_2,k(:,2));
% title('$\theta$ controller gain', 'Interpreter', 'latex', 'FontSize', 18);
% grid('on');
% xlabel('$Q_\theta/R_u $', 'Interpreter', 'latex', 'FontSize', 18);
% ylabel('$K_{\theta}$', 'Interpreter', 'latex', 'FontSize', 18);
% set(get(gca,'ylabel'),'rotation',0)
%
% figure
% plot(ratio_3,k(:,3));
% title('$\dot{x}$ controller gain', 'Interpreter', 'latex', 'FontSize', 18);
% grid('on');
% xlabel('$Q_{\dot{x}}/R_u $', 'Interpreter', 'latex', 'FontSize', 18);
% ylabel('$K_{\dot{x}}$', 'Interpreter', 'latex', 'FontSize', 18);
% set(get(gca,'ylabel'),'rotation',0)
%
% figure
% plot(ratio_4,k(:,4));
% title('$\dot{\theta}$ controller gain', 'Interpreter', 'latex', 'FontSize', 18);
% grid('on');
% xlabel('$Q_{\dot{\theta}}/R_u $', 'Interpreter', 'latex', 'FontSize', 18);
% ylabel('$K_{\dot\theta}$', 'Interpreter', 'latex', 'FontSize', 18);
% set(get(gca,'ylabel'),'rotation',0)
%
%

% %% checking stability for nonlinear system
% L_vec=3.587:0.001:3.59;
% 
% unstable= zeros(1,length(L_vec));
% result =[];
% for i = 1:length(L_vec)
%     X = ['Round ', num2str(i)];
%     disp(X)
%     L = L_vec(i);
% 
%     A = [0 0 1 0;
%         0 0 0 1;
%         -k/M -m*g/M -b/M 0;
%         k/(M*L) (m+M)/(M*L)*g b/(M*L) 0];
%     B = [0 0 1/M -1/(L*M)]';
%     H = [1 0 0 0;
%     0 1 0 0];
% E = [0 0;
%     0 0;
%     -1/M 0;
%     1/(L*M) -1/(m*L)];
% 
% 
%     sim('Model.slx');
% 
% 
%  load('stability.mat');
%  figure
%  plot(data(1,:), data(4,:)*180/pi)
% 
%   title([ ' State $\theta$ for L = ' num2str(L_vec(i)) ' m'] , 'Interpreter', 'latex', 'Fontsize', 18);
%   grid on
%   xlabel('Time [s]','Interpreter', 'latex', 'Fontsize', 18);
% ylabel('[rad]', 'Interpreter', 'latex', 'Fontsize', 18);
% 
% 
% 
% end
% %
% 
% 





%% checking stability for linear system.

L_vec=4.1778;
unstable= zeros(1,length(L_vec));
for i = 1:length(L_vec)
    L = L_vec(i);

    A = [0 0 1 0;
        0 0 0 1;
        -k/M -m*g/M -b/M 0;
        k/(M*L) (m+M)/(M*L)*g b/(M*L) 0];
    B = [0 0 1/M -1/(L*M)]';

    res = eig(A-B*K);
    

    for j=1:length(res)
        if real(res(1)) >0
            unstable(i) = 1;

        end
        if real(res(2)) >0
            unstable(i) = 1;

        end
        if real(res(3)) >0
            unstable(i) = 1;

        end
        if real(res(4)) >0
            unstable(i) = 1;
        end
    end

end


for k = 1:length(unstable)
   if(unstable(k) == 1)
    disp(k);
    break
   end
end



