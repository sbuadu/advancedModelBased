%% Input
TestKalman
%% Simulate
sim('Model_kalmanTest4.slx')
%% Plot

figure(1)
plot(Time,KalmanFilterOutLin(:,1));
hold on
plot(Time,KalmanFIlterOutNon(:,1));
grid on
xlabel('Time','interpreter','latex','FontSize',18)
ylabel('$x_{linear}$ and $x_{nonlinear}$ [m]','interpreter','latex','FontSize',18);
title('x - from linear and nonlinear system','interpreter','latex','FontSize',18);
lgd =legend({'$x_{linear}$',' $x_{nonlinear}$'},'interpreter','latex','FontSize',18);

figure(2)
plot(Time,KalmanFilterOutLin(:,1));
hold on
plot(Time,KalmanFIlterOutNon(:,1));
grid on
xlabel('Time','interpreter','latex','FontSize',18)
ylabel('$x_{linear}$ and $x_{nonlinear}$ [m]','interpreter','latex','FontSize',18);
title('x - from linear and nonlinear system','interpreter','latex','FontSize',18);
lgd =legend({'$x_{linear}$',' $x_{nonlinear}$'},'interpreter','latex','FontSize',18);

plot(Time,KalmanFilterOutLin(:,1));
hold on
plot(Time,KalmanFIlterOutNon(:,1));
grid on
xlabel('Time','interpreter','latex','FontSize',18)
ylabel('$x_{linear}$ and $x_{nonlinear}$ [m]','interpreter','latex','FontSize',18);
title('x - from linear and nonlinear system','interpreter','latex','FontSize',18);
lgd =legend({'$x_{linear}$',' $x_{nonlinear}$'},'interpreter','latex','FontSize',18);