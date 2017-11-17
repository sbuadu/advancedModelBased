close all;

%%Importing data 
load('impulse noise/expensive.mat');
%load('expensive.mat');
t_expensive = data(1,:);
u_expensive = data(2,:);
x_expensive = data(3,:);
theta_expensive = data(4,:);
clear data;

load('impulse noise/cheap.mat');
%load('cheap.mat');
t_cheap = data(1,:);
u_cheap = data(2,:);
x_cheap = data(3,:);
theta_cheap = data(4,:);


% 
% %Creating figures for expensive control
% figure; 
% hold on;
% plot(t_expensive,u_expensive);
% hold off;
% grid on;
% title('Expensive Control: Control input u', 'Interpreter', 'latex');
% xlabel('Time [s]', 'Interpreter', 'latex');
% ylabel('Force [N]', 'Interpreter', 'latex');
% 
% 
% figure; 
% hold on;
% plot( t_expensive, x_expensive);
% hold off;
% grid on;
% title('Expensive Control: State x', 'Interpreter', 'latex');
% xlabel('Time [s]','Interpreter', 'latex');
% ylabel('[m]', 'Interpreter', 'latex');
% 
% figure; 
% hold on;
% plot( t_expensive, theta_expensive);
% hold off;
% grid on;
% title('Expensive Control: State $\theta$', 'Interpreter', 'latex');
% xlabel('Time [s]','Interpreter', 'latex');
% ylabel('[rad]', 'Interpreter', 'latex');
% 
% 
% %Creating figures for cheap control
% 
% figure; 
% hold on;
% plot(t_cheap,u_cheap);
% hold off;
% grid on;
% title('Cheap Control: Control input u', 'Interpreter', 'latex');
% xlabel('Time [s]', 'Interpreter', 'latex');
% ylabel('Force [N]', 'Interpreter', 'latex');
% 
% 
% figure; 
% hold on;
% plot( t_cheap, x_cheap);
% hold off;
% grid on;
% title('Cheap Control: State x', 'Interpreter', 'latex');
% xlabel('Time [s]','Interpreter', 'latex');
% ylabel('[m]', 'Interpreter', 'latex');
% 
% figure; 
% hold on;
% plot( t_cheap, theta_cheap);
% hold off;
% grid on;
% title('Cheap Control: State $\theta$', 'Interpreter', 'latex');
% xlabel('Time [s]','Interpreter', 'latex');
% ylabel('[rad]', 'Interpreter', 'latex');

%%relevant plots

figure; 
hold on; 
plot( t_expensive, x_expensive);
plot( t_cheap, x_cheap);
hold off; 
grid on;
legend('Cheap control (Q/R= 10)','Expensive control (Q/R=0.01)');
title('Cheap vs expensive control:  State x', 'Interpreter', 'latex',  'FontSize', 14);
xlabel('Time [s]','Interpreter', 'latex', 'FontSize', 14);
ylabel('[m]', 'Interpreter', 'latex', 'FontSize', 14);
set(get(gca,'ylabel'),'rotation',0)



figure; 
hold on; 
plot( t_cheap, theta_cheap);
plot( t_expensive, theta_expensive);
hold off; 
grid on;
legend('Cheap control (Q/R= 10)','Expensive control (Q/R=0.01)')
title('Cheap vs expensive control: State $\theta$', 'Interpreter', 'latex',  'FontSize', 14);
xlabel('Time [s]','Interpreter', 'latex',  'FontSize', 14);
ylabel('[rad]', 'Interpreter', 'latex',  'FontSize', 14);
set(get(gca,'ylabel'),'rotation',0)


figure; 
hold on; 
plot(t_cheap,u_cheap);
 plot(t_expensive,u_expensive);
hold off; 
grid on;
legend('Cheap control (Q/R= 10)','Expensive control (Q/R=0.01)')
title('Cheap vs expensive control: Control input u', 'Interpreter', 'latex',  'FontSize', 14);
xlabel('Time [s]','Interpreter', 'latex',  'FontSize', 14);
ylabel('Force [N]', 'Interpreter', 'latex',  'FontSize', 14);



