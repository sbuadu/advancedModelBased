close all

%comparing x state before and after kalman Linear system
figure;
hold on;
title('Comparing state $x$ with and without filtering','interpreter','latex','FontSize',18);

plot(Time, NoiseSigLin(:,1));
plot(Time,KalmanFilterOutLin(:,1));

legend({'$x$',' filtered $x$'},'interpreter','latex','FontSize',18);
grid on
xlabel('Time','interpreter','latex','FontSize',18)
ylabel('$x$ [m]','interpreter','latex','FontSize',18);

axes('position',[.65 .175 .25 .25])
box on
plot(Time(30001:30300), NoiseSigLin(30001:30300,1));
hold on
plot(Time(30001:30300), KalmanFilterOutLin(30001:30300,1));
axis tight 



%comparing theta state before and after kalman Linear system
figure;
hold on;
title('Comparing state $theta$ with and without filtering','interpreter','latex','FontSize',18);

plot(Time, NoiseSigLin(:,2));
plot(Time,KalmanFilterOutLin(:,2));
legend({'$\theta$',' filtered $\theta$'},'interpreter','latex','FontSize',18);
grid on
xlabel('Time','interpreter','latex','FontSize',18)
ylabel('$\theta$ [rad]','interpreter','latex','FontSize',18);

axes('position',[.65 .175 .25 .25])
box on
plot(Time(30001:30300), NoiseSigLin(30001:30300,2));
hold on
plot(Time(30001:30300), KalmanFilterOutLin(30001:30300,2));
axis tight 




