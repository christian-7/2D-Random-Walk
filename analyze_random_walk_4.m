% 1. simulates a random walk according to the initialized parameters
% 2. calculates trajectory parameters: scatter, MSD, velocity, distance from origin
% 3. calculates the confinement index according to Simson, et al. (1995), Biophysical Journal, 69(September), 989?993. 
% 
% 
% D-diffusion coefficient
% dt-time step
% 
% Date:     04/05/15
% Author:   Christian Sieben

clear all, clc, close all

%% Initiate parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start=[1 1];                                % starting coordinates
num_steps=300;                              % number of steps
D=0.046;                                     % diffusion constant ?m2/s             
segment=20;                                 % Sm, segment length in frames


dt=0.3;                     % time step
dx=1;                       % pixel size
step_size=sqrt(4*D*dt);     % D=(dx^2)/dt --> in ?m  

minAxis=-2.5;               % min axis 
maxAxis=2.5;                % max axis


figure('Position',[200 400 900 650],'name','Overview Figure: scatter, velocity, MSD, cum. displacement')
h=gcf;
set(h,'PaperOrientation','landscape');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate random walk

pos = zeros(num_steps,2);
pos(1,1)=start(1);
pos(1,2)=start(2);
pos(1,3)=0;

a = 1;       % randomize angle
b = 360;  
angle = (b-a).*rand(num_steps,1) + a;

for k=2:num_steps;
    
    pos(k,1)=(pos(k-1,1))+sin(angle(k))*step_size;      % pos in ?m
    pos(k,2)=(pos(k-1,2))+cos(angle(k))*step_size;      % pos in ?m
    pos(k,3)=k*dt;                                      % time in seconds
    
end

%% Plot Data

subplot(2,3,1)
line(pos(:,1)*dx,pos(:,2)*dx);hold on;
scatter(pos(:,1)*dx,pos(:,2)*dx,3,pos(:,3));hold on;
plot(pos(1,1)*dx,pos(1,2)*dx,'*b','MarkerSize',12);hold on;
text(pos(1,1)*dx,pos(1,2)*dx, 'Start');
plot(pos(length(pos),1)*dx,pos(length(pos),2)*dx,'+b','MarkerSize',12);hold on;
text(pos(length(pos),1)*dx,pos(length(pos),2)*dx,'End');hold on;
title('XY scatter trajectory');
axis([minAxis maxAxis minAxis maxAxis])
xlabel('x (\mu m)', 'FontSize',12);
ylabel('y (\mu m)', 'FontSize',12);
box on;
% c=colorbar%('northoutside');


%% Plot displacement from origin

dcum = zeros(num_steps,1);
dcum(1,1)=0;

for k=2:num_steps;
    
    dist(k,1)=sqrt(((pos(k,1)-pos(1,1))^2)+((pos(k,2)-pos(1,2))^2));

    dcum(k)=dist(k)+dcum(k-1);
end

subplot(2,3,4)
scatter(pos(:,3),dist*dx,2,pos(:,3));
title('Distance from origin');
xlabel('time (s)', 'FontSize',12);
ylabel('distance from origin (\mu m)', 'FontSize',12);

subplot(2,3,5)
scatter(pos(:,3),dcum*dx,2,pos(:,3));
set(gca,'xscale','log')
set(gca,'yscale','log')
axis([1 100 0.001 1e5])
title('Cumulative Distance from Origin', 'FontSize',12);
xlabel('time (s)', 'FontSize',12);
ylabel('cumulative distance from origin (\mu m)', 'FontSize',12);

%% My MSD

frame=pos(:,3);
frame=frame/dt; % random walk generate in s --> transform to frames
frame=frame-min(frame);

% i = frame --> Reihe
% j = gap; --> Spalte

msd=zeros(max(frame), 50);
msd2=zeros(50, 3);

for i=1:1:max(frame);    % for all frames
    vx=find(frame == i);
    
    if isempty(vx)==1; % if frame does not exist, skip   
    else
        
    
    for j=1:300; % time gap
    
    if vx+j>length(pos) || i+j~=frame(vx+j) % if point plus gap exeeds lentgh, skip
    
        msd(i,j)=0;
    
    else
        
        msd(i,j)=((pos(vx,1)-pos(vx+j,1))^2)+((pos(vx,2)-pos(vx+j,2))^2);
    
%         msd2(j,1)=j;
%         msd2(j,2)=mean(nonzeros(msd(:,j)));
%         msd2(j,3)=std(nonzeros(msd(:,j)));
    end
    
    end
        
    end
    
end

for m=1:300;
    
        msd2(m,1)=m;
        msd2(m,2)=mean(nonzeros(msd(:,m)));
        msd2(m,3)=std(nonzeros(msd(:,m)));
        
end


msd2(:,4)=msd2(:,2);        % mean msd per second
msd2(:,5)=msd2(:,1)*dt;     % lag time in second
msd2(:,6)=msd2(:,3);        % error msd per second

subplot(2,3,2)
% scatter(msd2(:,5),msd2(:,4)); hold on;
% errorbar(msd2(:,5),msd2(:,4),msd2(:,6)); hold on;
plot(msd2(:,5),msd2(:,4)); hold on;
plot(msd2(:,5),4*D*msd2(:,5),'-g')

%% Velocity

vel=zeros(length(pos),1);

for k=2:num_steps;
    
    vel(k,1)=sqrt(((pos(k,1)-pos(k-1,1))^2)+((pos(k,2)-pos(k-1,2))^2));

    
end

vel=vel*dx/dt; % vel in ?m/s
subplot(2,3,3)
scatter(pos(:,3),vel,3,pos(:,3));
title('Velocity');
xlabel('time (s)');
ylabel('velocity (\mu m/s)');


subplot(2,3,6)
hist(vel,20)
title('Histogram of Velocity');
xlabel('velocity (\mu m/s)', 'FontSize',12);
ylabel('count', 'FontSize',12);

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

figure('Position',[200 20 900 300],'name','Confinement parameters')
h=gcf;
set(h,'PaperOrientation','landscape');

subplot(1,3,1)
plot(D*prob2(:,1)*dt,prob2(:,2));hold on;
xlabel('Dt/R^2','FontSize',12);
ylabel('mean log(\psi)','FontSize',12);

subplot(1,3,2)
plot(prob2(:,1)*dt,10.^(prob2(:,2)));hold on;
xlabel('time (s)','FontSize',12);
ylabel('mean \psi','FontSize',12);

subplot(1,3,3)
plot(L(:,2)*dt, L(:,1)); hold on;       
xlabel('time (s)','FontSize',12);
ylabel('probability level L','FontSize',12);

%% Figures to save

figure('Position',[200 800 300 300], 'name','confinement index L')

plot(L(:,2)*dt, L(:,1)); hold on;
% plot(L(:,2)*dt, s3,'-r');       
xlabel('time (s)','FontSize',12);
ylabel('confinement index L','FontSize',12);


figure('Position',[600 800 600 300], 'name','XY scatter and confinement index L')
h=gcf;
set(h,'PaperOrientation','landscape');

subplot(1,2,1)
line(pos(:,1)*dx,pos(:,2)*dx);hold on;
scatter(pos(:,1)*dx,pos(:,2)*dx,10,pos(:,3),'filled');hold on;
plot(pos(1,1)*dx,pos(1,2)*dx,'*b','MarkerSize',12);hold on;
text(pos(1,1)*dx,pos(1,2)*dx, 'Start');
plot(pos(length(pos),1)*dx,pos(length(pos),2)*dx,'+b','MarkerSize',12);hold on;
text(pos(length(pos),1)*dx,pos(length(pos),2)*dx,'End');hold on;
axis([minAxis maxAxis minAxis maxAxis])
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
%% 

figure('Position',[500 800 500 600], 'name','Compare time with displacement and confinement')
h=gcf;
set(h,'PaperOrientation','portrait');

subplot(4,1,1)
line(pos(:,3)*dt,pos(:,1));hold on;
scatter(pos(:,3)*dt,pos(:,1),5,pos(:,3)*dt);hold on;
xlabel('time (s)','FontSize',12);
ylabel('x (\mu m)','FontSize',12);
box on;



subplot(4,1,2)
line(frame*dt,pos(:,2));hold on;
scatter(frame*dt,pos(:,2),5,pos(:,3));hold on;
xlabel('time (s)','FontSize',12);
ylabel('y (\mu m)','FontSize',12);
box on;

subplot(4,1,3)
plot(L(:,2)*dt, L(:,1)); hold on;
scatter(L(:,2)*dt, L(:,1),15,pos(:,3),'filled')     
xlabel('time (s)','FontSize',12);
ylabel('confinement index L','FontSize',12);
box on;

dcum = zeros(num_steps,1);
dcum(1,1)=0;

for k=2:num_steps;
    
    dist(k,1)=sqrt(((pos(k,1)-pos(1,1))^2)+((pos(k,2)-pos(1,2))^2));

    dcum(k)=dist(k)+dcum(k-1);
end


subplot(4,1,4)
plot(frame*dt,dist); hold on;
scatter(frame*dt,dist,2,pos(:,3)*dt);
xlabel('time (s)','FontSize',12);
ylabel('distance from origin (\mu m)','FontSize',12);
box on;


