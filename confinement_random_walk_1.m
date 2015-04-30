clear all, clc, close all
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D=0.000605085473694187;                     % diffusion constant, ?m2/s             
t=1:1:1000;                                 % number of time steps, seconds        
R=1;                                        % radius of the confinement zone, ?m

start=[1 1];                % starting coordinates
num_steps=250;              % number of steps
                    
dt=0.5;                     % time step, seconds
dx=1;                     % pixel size, ?m
step_size=sqrt(4*D*dt);     % step size from D, ?m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate from Simson et al, 1995

omega=0.2048-2.5117*((D*t)./(R^2));
% 
% figure('Position',[20 400 700 300])
% 
% h=gcf;
% set(h,'PaperOrientation','landscape');

subplot(1,2,1)
plot(t,10.^omega); hold on;
% title('Probability to stay in confined region with radius R (Simson et al, 1995)')
xlabel('time (s)','FontSize',12);
ylabel('\psi (R,t)','FontSize',12);
box on;

subplot(1,2,2)
plot(D*t/R^2,omega);hold on;
xlabel('Dt/R^2','FontSize',12);
ylabel('log(\psi)','FontSize',12);
box on;

%% Generate random walk Trajectory

pos = zeros(num_steps,2);
pos(1,1)=start(1);
pos(1,2)=start(2);
pos(1,3)=0;

a = 1;      
b = 360;  
angle = (b-a).*rand(num_steps,1) + a;

for k=2:num_steps;
    
    pos(k,1)=(pos(k-1,1))+sin(angle(k))*step_size;
    pos(k,2)=(pos(k-1,2))+cos(angle(k))*step_size;
    pos(k,3)=k*dt; % time in seconds
    
end

% subplot(2,3,1)
figure('Position',[730 400 400 300])
line(pos(:,1)*dx,pos(:,2)*dx);hold on;
scatter(pos(:,1)*dx,pos(:,2)*dx,3,pos(:,3));hold on;
plot(pos(1,1)*dx,pos(1,2)*dx,'*b','MarkerSize',12);hold on;
text(pos(1,1)*dx,pos(1,2)*dx, 'Start');
plot(pos(length(pos),1)*dx,pos(length(pos),2)*dx,'+b','MarkerSize',12);hold on;
text(pos(length(pos),1)*dx,pos(length(pos),2)*dx,'End');hold on;
title('XY scatter trajectory');
xlabel('x (\mu m)');
ylabel('y (\mu m)');

%% Calculate confinement

% generate variable frame

frame=pos(:,3);             % time step in seconds
frame=frame/dt;             % time step in frames
frame=frame-min(frame);     % starting from 0
frame(1,1)=1;               % starting from 1

% i = frame --> Reihe
% j = gap; --> Spalte

prob=[];    %zeros(max(frame), 5);
prob2=[];   %zeros(5, 2);
d=[];
vx=[];
vy=[];


for i=1:max(frame)-4;    % for all frames
    vx=find(frame == i);
    
    if isempty(vx)==1;                   % if frame does not exist, skip   
    else
        
    c=1; 
    for j=1:4;                                         % segment length
          
        vy=find(frame <= (i+j) & frame >= i );         % select segment
        subset(:,1)=pos(vy);                           % define segment as subset
        subset(:,2)=pos(vy,2);
        
        if length(vy)==1;   % if subset is only 1 frame --> distance is 0
                     R=0;
        else    
        
            for k=2:length(subset);
            
                 d(k,1)=sqrt(((subset(k,1)-subset(1,1))^2)+((subset(k,2)-subset(1,2))^2));   
                 R=max(d);
        
            end
        
                                                                % maximum distance within subset
        prob(i:(i+j),c)=0.2048-2.5117*((D*j)./(R^2));               % probability within subset

        clear subset
        end
    
    c=c+1;    
        
    end
    clear vx vy R d c;
    
    end
   
end
clear subset

for l=1:length(prob)
    
prob2(l,1)=l;
prob2(l,2)=mean(nonzeros(prob(l,:))); % this is omega

end


L=[];

for i=1:length(prob2)
    
    if 10.^(prob2(i,2))>0.1
       L(i,1)=0;
       
    else
        
        L(i,1)=((prob2(i,2))*(-1)-1);
        
    end
L(i,2)=i;
    
end

% s1=smooth(D*prob2(:,1)*dt,prob2(:,2),0.05,'loess');
% s2=smooth(prob2(:,1)*dt,10.^(prob2(:,2)),0.05,'loess');
% s3=smooth(L(:,2)*dt, L(:,1),0.05,'loess');


figure('Position',[200 20 900 300])

subplot(1,3,1)
% plot(D*prob2(:,1)*dt,prob2(:,2));hold on;
% plot(D*prob2(:,1)*dt,s1,'-r');
plot(prob2(:,1)*dt,prob2(:,2));hold on;
% plot(prob2(:,1)*dt,s1,'-r');
xlabel('time (s)','FontSize',12);
ylabel('mean log(\psi)','FontSize',12);

subplot(1,3,2)
plot(prob2(:,1)*dt,10.^(prob2(:,2)));hold on;
% plot(prob2(:,1)*dt,s2,'-r');
xlabel('time (s)','FontSize',12);
ylabel('mean \psi','FontSize',12);

subplot(1,3,3)
plot(L(:,2)*dt, L(:,1)); hold on;
% plot(L(:,2)*dt, s3,'-r');       
xlabel('time (s)','FontSize',12);
ylabel('probability level L','FontSize',12);


%% Smoothed curves

% figure('Position',[200 20 900 300])
% subplot(1,3,1)
% % plot(D*prob2(:,1)*dt,prob2(:,2));hold on;
% plot(D*prob2(:,1)*dt,s1,'-r');
% xlabel('Dt/R^2','FontSize',12);
% ylabel('mean log(\psi)','FontSize',12);
% 
% subplot(1,3,2)
% % plot(prob2(:,1)*dt,10.^(prob2(:,2)));hold on;
% plot(prob2(:,1)*dt,s2,'-r');
% xlabel('time (s)','FontSize',12);
% ylabel('mean \psi','FontSize',12);
% 
% subplot(1,3,3)
% % plot(L(:,2)*dt, L(:,1)); hold on;
% plot(L(:,2)*dt, s3,'-r');       
% xlabel('time (s)','FontSize',12);
% ylabel('probability level L','FontSize',12);