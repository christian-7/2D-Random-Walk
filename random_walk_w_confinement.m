% 1. simulates a random walk according to the initialized parameters and adds confined areas with D_conf
% 2. calculates the confinement index according to Simson, et al. (1995), Biophysical Journal, 69(September), 989�993. 
% 
% 
% D-diffusion coefficient
% D_conf-diffusion coefficient in confined area
% dt-time step
% 
% Date: 30/04/15
% Author: Christian Sieben

clear all, close all, clc

%% Initialize parameters for random walk and add circles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D=0.000605085473694187;                     % diffusion constant, ?m2/s 
D_conf=D/5;                                 % diffusion constant in confined areas, ?m2/s

start=[0 0];                                % starting coordinates
num_steps=300;                              % number of steps
                    
dt=0.3;                                     % time step, seconds
dx=1;                                       % pixel size, ?m
step_size=sqrt(4*D*dt);                     % step size according to D
step_size_conf=sqrt(4*D_conf*dt);           % step size according to D_conf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define confined areas as circles with radius and center position

radius=0.05;
center=[0.1 0.1];
center2=[-0.2 0.2];
center3=[-0.2 -0.1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate and plot random walk trajectory

pos = zeros(num_steps,2);
pos(1,1)=start(1);
pos(1,2)=start(2);
pos(1,3)=0;

a = 1;      
b = 360;  
angle = (b-a).*rand(num_steps,1) + a;

for k=2:num_steps;
    
    if  sqrt(((pos(k-1,1)-center(1,1))^2)+((pos(k-1,2)-center(1,2))^2)) <= radius | ...
        sqrt(((pos(k-1,1)-center2(1,1))^2)+((pos(k-1,2)-center2(1,2))^2)) <= radius | ...
        sqrt(((pos(k-1,1)-center3(1,1))^2)+((pos(k-1,2)-center3(1,2))^2)) <= radius;
        
    pos(k,1)=(pos(k-1,1))+sin(angle(k))*step_size_conf;
    pos(k,2)=(pos(k-1,2))+cos(angle(k))*step_size_conf;
    pos(k,3)=k*dt;
    pos(k,4)=step_size_conf;
    
    else
        
    pos(k,1)=(pos(k-1,1))+sin(angle(k))*step_size;
    pos(k,2)=(pos(k-1,2))+cos(angle(k))*step_size;
    pos(k,3)=k*dt; % time in seconds
    pos(k,4)=step_size;

    end
    
end


figure('Position',[0 400 400 300],'name','XY scatter trajectory')
line(pos(:,1)*dx,pos(:,2)*dx);hold on;
scatter(pos(:,1)*dx,pos(:,2)*dx,3,pos(:,3));hold on;
plot(pos(1,1)*dx,pos(1,2)*dx,'*b','MarkerSize',12);hold on;
text(pos(1,1)*dx,pos(1,2)*dx, 'Start');
plot(pos(length(pos),1)*dx,pos(length(pos),2)*dx,'+b','MarkerSize',12);hold on;
text(pos(length(pos),1)*dx,pos(length(pos),2)*dx,'End');hold on;
well1 = viscircles(center, radius);
well2 = viscircles(center2, radius);
well3 = viscircles(center3, radius);
axis([-0.5 0.5 -0.5 0.5])
title('XY scatter trajectory');
xlabel('x (\mum)','FontSize',12);
ylabel('y (\mum)','FontSize',12);
box on;

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

c=1;
segment=19; % Sm, segment length in frames


for i=1:max(frame)-segment;                       % for all frames
    vx=find(frame == i);                     % find frame i
    
    if isempty(vx)==1;                       % if frame does not exist, skip   
    else
        
     
    for j=4:segment;                                           % segment length
          
        vy=find(frame <= (i+j) & frame >= i );          % select segment
        subset(:,1)=pos(vy);                            % define segment as subset
        subset(:,2)=pos(vy,2);
        
        if length(vy)==1;   % if subset is only 1 frame --> distance is 0
                     R=0;
        else    
        
            for  k=2:length(subset);
                 d(k,1)=sqrt(((subset(k,1)-subset(1,1))^2)+((subset(k,2)-subset(1,2))^2));    % calculate the distance to each point in subset from point i  
            end
        R=max(d);                                                      % maximum distance within subset
        prob(i:(i+j),c)=0.2048-2.5117*((D*j)./(R^2));                  % probability within subset
%         prob(c,i:(i+j))=((D*j)./(R^2));

%       prob(c,i:(i+j))=horzcat(prob(c,i:(i+j)),((D*j)./(R^2)));
        c=c+1;  
        clear subset
        end
    
%     c=c+1;    
        
    end
    clear vx vy R d;
    
    end
   
end
clear subset

for l=1:length(frame)
    
prob2(l,1)=l;                           % frame
prob2(l,2)=mean(nonzeros(prob(l,:)));   % this is psi

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


figure('Position',[400 400 900 300], 'name','Probability psi, log(psi), confinement index L')

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


figure('Position',[200 800 300 300], 'name','confinement index L')
plot(L(:,2)*dt, L(:,1)); hold on;
% plot(L(:,2)*dt, s3,'-r');       
xlabel('time (s)','FontSize',12);
ylabel('confinement index L','FontSize',12);


%% Figures to save

figure('Position',[200 800 300 300], 'name','confinement index L')
plot(L(:,2)*dt, L(:,1)); hold on;
% plot(L(:,2)*dt, s3,'-r');       
xlabel('time (s)','FontSize',12);
ylabel('confinement index L','FontSize',12);


figure('Position',[600 800 600 300], 'name','XY scatter and confinement index L')
subplot(1,2,1)
line(pos(:,1)*dx,pos(:,2)*dx);hold on;
scatter(pos(:,1)*dx,pos(:,2)*dx,10,pos(:,3),'filled');hold on;
plot(pos(1,1)*dx,pos(1,2)*dx,'*b','MarkerSize',12);hold on;
text(pos(1,1)*dx,pos(1,2)*dx, 'Start');
plot(pos(length(pos),1)*dx,pos(length(pos),2)*dx,'+b','MarkerSize',12);hold on;
text(pos(length(pos),1)*dx,pos(length(pos),2)*dx,'End');hold on;
well1 = viscircles(center, radius);
well2 = viscircles(center2, radius);
well3 = viscircles(center3, radius);
axis([-0.5 0.5 -0.5 0.5])
title('XY scatter trajectory');
xlabel('x (\mum)','FontSize',12);
ylabel('y (\mum)','FontSize',12);
box on;



subplot(1,2,2)
plot(L(:,2)*dt, L(:,1)); hold on;
scatter(L(:,2)*dt, L(:,1),15,pos(:,3),'filled')
% plot(L(:,2)*dt, s3,'-r');       
xlabel('time (s)','FontSize',12);
ylabel('confinement index L','FontSize',12);
title('Confinement index L');