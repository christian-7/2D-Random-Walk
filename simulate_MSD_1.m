D=0.001;                % diffusion constant ?m2/s             
t=1:1:1000;    % time
v=0.0001; % speed of active transport

confA=0.1;
A1=1.01;
A2=0.71;

alpha=0.7;


%% Simulate

normal=4*D*t;

active=4*D*t+v*t.^2;

confined=confA*(1-A1*exp((-4*A2*D*t)/confA));

anomalous=4*D*t.^alpha;

figure('Position',[200 400 1000 200])

subplot(1,4,1)
plot(t,normal,'b');hold on;
title('normal');
% title('MSD=4Dt','FontSize',14');
xlabel('lag time t');
ylabel('MSD r^2');
% legend('normal')
% text(50,3,'MSD=4*D*t','FontSize',12)


subplot(1,4,2)
plot(t,anomalous,'black');hold on;
% title('MSD=4*D*t^\alpha','FontSize',14');
title('anomalous');
xlabel('lag time t');
ylabel('MSD r^2');
% legend('anomalous')
% text(50,0.4,'MSD=4*D*t^\alpha','FontSize',12)



subplot(1,4,3)
plot(t,active,'r');hold on;
% title('MSD=4Dt+vt^2','FontSize',14);
title('active');
xlabel('lag time t');
ylabel('MSD r^2');
% legend('active')
% text(50,90,'MSD=4*D*t+v*t.^2','FontSize',12)


subplot(1,4,4)
plot(t,confined,'g');hold on;
% title('MSD=r_c(1-A1exp((-4A2Dt)/confA))','FontSize',14);
title('confined');
xlabel('lag time t');
ylabel('MSD r^2');
% text(50,9.91e-3,'MSD=confA*(1-A1*exp((-4*A2*D*t)/confA))','FontSize',8);

figure('Position',[50 100 300 300])
plot(t,normal/max(normal),'b');hold on;
plot(t,anomalous/(max(anomalous)*1.5),'black');hold on;
plot(t,active/(max(active)*0.5),'r');hold on;
plot(t,confined/(max(confined)*1.5),'g');hold on;
xlabel('lag time (s)','FontSize',14);
ylabel('MSD (r^2)','FontSize',14);
legend('normal','anomalous','active','confined','FontSize',14)


