clear all, clc, close all

%% Initiate parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start=[1 1];            % starting coordinates
num_steps=800;          % number of steps
D=1e-3;                 % diffusion constant ?m2/s             


dt=0.5;                     % time step
dx=0.1;                     % pixel size
step_size=sqrt(4*D*dt);     % D=(dx^2)/dt --> in ?m  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate random walk in loop

%%%
n=50;
random_walks=zeros(num_steps,n*2); % matrix 
%%%
 j=1;

for i=1:n;

   
    
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

random_walks(:,j)=pos(:,1);
random_walks(:,j+1)=pos(:,2);

j=j+2;

clear pos a b angle k
end

%% Calculate My MSD

l=1;
k=1;
% n=10;
msd_all=zeros(100,n*3); % matrix 


for l=1:2:2*n;
    
    pos(:,1)=random_walks(:,l);
    pos(:,2)=random_walks(:,l+1);
    pos(:,3)=1:num_steps;

    frame=pos(:,3);
    
% frame=frame/dt; % random walk generate in s --> transform to frames
% frame=frame-min(frame);

% i = frame --> Reihe
% j = gap; --> Spalte

msd=zeros(max(frame), 50);
msd2=zeros(50, 3);

for i=1:1:max(frame);    % for all frames
    vx=find(frame == i);
    
    if isempty(vx)==1; % if frame does not exist, skip   
    else
        
    
    for j=1:100; % time gap
    
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

for m=1:100;
    
        msd2(m,1)=m;
        msd2(m,2)=mean(nonzeros(msd(:,m)));
        msd2(m,3)=std(nonzeros(msd(:,m)));        
end



msd_all(:,k)=msd2(:,2);          % mean msd per second
msd_all(:,k+1)=msd2(:,1)*dt;     % lag time in second
msd_all(:,k+2)=msd2(:,3);        % error msd per second

k=k+3;
clear msd msd2 pos frame

end

w=1;
figure;
for v=1:n;
    
    plot(msd_all(:,w+1),msd_all(:,w)); hold on;
    w=w+3;
end

plot(msd_all(:,2),4*D*msd_all(:,2),'og')
title('Mean squared displacement')
xlabel('lag time (s)','FontSize',12);
ylabel('MSD (\mu m^2)','FontSize',12);

%% Calculate confinement

l=1;
u=1;
% n=100;
% conf_all=zeros(num_steps,n); % matrix 


for l=1:2:(2*n)-1;
    
    pos(:,1)=random_walks(:,l);
    pos(:,2)=random_walks(:,l+1);
    pos(:,3)=1:num_steps;

    frame=pos(:,3);
    
%     l=l+2

% generate variable frame

% frame=pos(:,3);             % time step in seconds
% frame=frame/dt;             % time step in frames
% frame=frame-min(frame);     % starting from 0
% frame(1,1)=1;               % starting from 1

% i = frame --> Reihe
% j = gap; --> Spalte

prob=[];    %zeros(max(frame), 5);
prob2=[];   %zeros(5, 2);
d=[];
vx=[];
vy=[];


for i=1:max(frame)-30;    % for all frames
    vx=find(frame == i);
    
    if isempty(vx)==1;                   % if frame does not exist, skip   
    else
        
    c=1; 
    for j=1:30;                                         % segment length
          
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
prob2(l,2)=mean(nonzeros(prob(l,:))); % this is log omega

end
  
conf_all(:,u)=prob2(:,2); % this is omega


u=u+1;

clear prob prob2

end

%% Plot probability level

L=[];
frames=1:length(conf_all);
figure('Position',[200 20 900 300])
h=gcf;
set(h,'PaperOrientation','landscape');



for j=1:n;

for i=1:length(conf_all)
    
    if 10.^(conf_all(i,j))>0.1
       L(i,j)=0;
       
    else
        
        L(i,j)=((conf_all(i,j))*(-1)-1);
        
    end
        
end

s1(:,j)=smooth(D*frames*dt,conf_all(:,j),0.05,'loess');
s2(:,j)=smooth(frames*dt,10.^(conf_all(:,j)),0.05,'loess');
s3(:,j)=smooth(frames*dt, L(:,j),0.05,'loess');


subplot(1,3,1)
% plot(D*prob2(:,1)*dt,prob2(:,2));hold on;
% plot(D*frames*dt,s1(:,j),'-r');hold on;
plot(frames*dt,s1(:,j));hold on;
xlabel('time(s)','FontSize',12);
ylabel('mean log(\psi)','FontSize',12);

subplot(1,3,2)
% plot(prob2(:,1)*dt,10.^(prob2(:,2)));hold on;
plot(frames*dt,s2(:,j));hold on;
xlabel('time (s)','FontSize',12);
ylabel('mean \psi','FontSize',12);

subplot(1,3,3)
% plot(L(:,2)*dt, L(:,1)); hold on;
plot(frames*dt, s3(:,j));hold on;   
xlabel('time (s)','FontSize',12);
ylabel('probability level L','FontSize',12);

end