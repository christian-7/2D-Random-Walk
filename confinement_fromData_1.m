
clear all, clc, close all

%% Load data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('tracks2.mat')

pos(:,1)=nonzeros(tracks.Ch2(:,2)*0.1);                                      % x coord
pos(:,2)=nonzeros(tracks.Ch2(:,3)*0.1);                                      % y coord
% help(1,1)=1;
help=(tracks.Ch2(:,1))+1;
% pos(:,3)=nonzeros(tracks.Ch2(:,1))-(min(nonzeros(tracks.C2(:,1))));          % set time in sec (from data)
pos(:,3)=help(1:length(pos),1);  

% pos(:,1)=(tracks.Ch2(:,2)*0.1);                                      % x coord
% pos(:,2)=(tracks.Ch2(:,3)*0.1);                                      % y coord       
% pos(:,3)=(tracks.Ch2(:,1)-min((tracks.C2(:,1))));  

num_steps=length(pos);
D=0.001;

dt=0.5; % time step, seconds
dx=0.1; % pixel size, ?m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Data

figure('Position',[200 400 400 300])
line(pos(:,1)*dx,pos(:,2)*dx);hold on;
scatter(pos(:,1)*dx,pos(:,2)*dx,3,pos(:,3));hold on;
plot(pos(1,1)*dx,pos(1,2)*dx,'*b','MarkerSize',12);hold on;
text(pos(1,1)*dx,pos(1,2)*dx, 'Start');
plot(pos(length(pos),1)*dx,pos(length(pos),2)*dx,'+b','MarkerSize',12);hold on;
text(pos(length(pos),1)*dx,pos(length(pos),2)*dx,'End');hold on;
title('XY scatter trajectory','FontSize',12);
xlabel('x (\mu m)','FontSize',12);
ylabel('y (\mu m)','FontSize',12);                  % pixel size
colorbar



%% Calculate confinement

% generate variable frame

frame=pos(:,3);             % time step in seconds
frame=frame/dt;             % time step in frames
frame=frame-min(frame);     % starting from 0
frame(1,1)=1;               % starting from 1

% i = frame --> Reihe
% j = gap; --> Spalte

prob=[];%zeros(max(frame), 5);
prob2=[];%zeros(5, 2);
d=[];
vx=[];
vy=[];


for i=1:max(frame);    % for all frames
    vx=find(frame == i);
    
    if isempty(vx)==1;                   % if frame does not exist, skip   
    else
        
    c=1; 
    for j=1:30;                                         % segment length
        
        subset=[];
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
%         prob(i:(i+j),c)=(0.2048-2.5117*D*j./1*R^2);               % probability within subset
          prob(i:(i+j),c)=0.2048-2.5117*((D*j)./(R^2));
        
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

L=[];

for i=1:length(prob2)
    
    if 10.^(prob2(i,2))>0.1
       L(i,1)=0;
    else   
        L(i,1)=((prob2(i,2))*(-1)-1);   
    end
L(i,2)=i;
    
end

s1=smooth(D*prob2(:,1)*dt,prob2(:,2),0.05,'rloess');
s2=smooth(prob2(:,1)*dt,10.^(prob2(:,2)),0.05,'rloess');
s3=smooth(L(:,2)*dt, L(:,1),0.05,'moving');


figure('Position',[200 20 900 300])
subplot(1,3,1)
plot(D*prob2(:,1)*dt,prob2(:,2));hold on;
plot(D*prob2(:,1)*dt,s1,'-r');
xlabel('Dt/R^2','FontSize',12);
ylabel('mean log(\psi)','FontSize',12);

subplot(1,3,2)
plot(prob2(:,1)*dt,10.^(prob2(:,2)));hold on;
plot(prob2(:,1)*dt,s2,'-r');
xlabel('time (s)','FontSize',12);
ylabel('mean \psi','FontSize',12);

subplot(1,3,3)
plot(L(:,2)*dt, L(:,1)); hold on;
plot(L(:,2)*dt, s3,'-r');       
xlabel('time (s)','FontSize',12);
ylabel('probability level L','FontSize',12);

%% Plot smoothed curves

s1=smooth(D*prob2(:,1)*dt,prob2(:,2),0.05,'rloess');
s2=smooth(prob2(:,1)*dt,10.^(prob2(:,2)),0.05,'rloess');
s3=smooth(L(:,2)*dt, L(:,1),0.05,'moving');


figure('Position',[200 20 900 300])
subplot(1,3,1)
% plot(D*prob2(:,1)*dt,prob2(:,2));hold on;
plot(D*prob2(:,1)*dt,s1,'-r');
xlabel('Dt/R^2','FontSize',12);
ylabel('mean log(\psi)','FontSize',12);

subplot(1,3,2)
% plot(prob2(:,1)*dt,10.^(prob2(:,2)));hold on;
plot(prob2(:,1)*dt,s2,'-r');
xlabel('time (s)','FontSize',12);
ylabel('mean \psi','FontSize',12);

subplot(1,3,3)
% plot(L(:,2)*dt, L(:,1)); hold on;
plot(L(:,2)*dt, s3,'-r');       
xlabel('time (s)','FontSize',12);
ylabel('probability level L','FontSize',12);