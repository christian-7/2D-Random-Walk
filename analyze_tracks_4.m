clear all, clc, close all
%% Load data

load('tracks2.mat')

pos(:,1)=nonzeros(tracks.Ch4(:,2)*0.1);                                      % x coord
pos(:,2)=nonzeros(tracks.Ch4(:,3)*0.1);                                      % y coord
% help(1,1)=1;
help=(tracks.Ch4(:,1))+1;
% pos(:,3)=nonzeros(tracks.Ch2(:,1))-(min(nonzeros(tracks.C2(:,1))));          % set time in sec (from data)
pos(:,3)=help(1:length(pos),1);  

% pos(:,1)=(tracks.Ch2(:,2)*0.1);                                      % x coord
% pos(:,2)=(tracks.Ch2(:,3)*0.1);                                      % y coord       
% pos(:,3)=(tracks.Ch2(:,1)-min((tracks.C2(:,1))));  


dt=0.5; % time step
dx=0.1; % pixel size

num_steps=length(pos);
D=0.001;

figure('Position',[200 400 900 500])
h=gcf;
set(h,'PaperOrientation','landscape');

%% Plot Data

subplot(2,3,1)
line(pos(:,1)*dx,pos(:,2)*dx);hold on;
scatter(pos(:,1)*dx,pos(:,2)*dx,3,pos(:,3));hold on;
plot(pos(1,1)*dx,pos(1,2)*dx,'*b','MarkerSize',12);hold on;
text(pos(1,1)*dx,pos(1,2)*dx, 'Start');
plot(pos(length(pos),1)*dx,pos(length(pos),2)*dx,'+b','MarkerSize',12);hold on;
text(pos(length(pos),1)*dx,pos(length(pos),2)*dx,'End');hold on;
title('XY scatter trajectory');
xlabel('x (\mu m)','FontSize',12);
ylabel('y (\mu m)','FontSize',12);
colorbar


%% Plot displacement from origin

dcum = zeros(num_steps,1);
dcum(1,1)=0;

for k=2:num_steps;
    
    d(k,1)=sqrt(((pos(k,1)-pos(1,1))^2)+((pos(k,2)-pos(1,2))^2));

    dcum(k)=d(k)+dcum(k-1);
end

subplot(2,3,4)
scatter(pos(:,3),d*dx,2,pos(:,3));
title('Distance from origin');
xlabel('time (s)','FontSize',12);
ylabel('distance from origin (\mu m)','FontSize',12);

subplot(2,3,5)
scatter(pos(:,3),dcum*dx,2,pos(:,3));
set(gca,'xscale','log')
set(gca,'yscale','log')
% axis([1 100 0.001 1e5])
title('Cumulative Distance from Origin');
xlabel('time (s)','FontSize',12);
ylabel('cumulative distance from origin (\mu m)','FontSize',12);

%% Calculate MSD (from Nathanael)

tab(:,1)=pos(:,3); %1:1:length(pos);
tab(:,2)=pos(:,1);
tab(:,3)=pos(:,2);
X=tab(:,2:3);

% dx=0.1;     % pixel size
% dt=0.5;       % time step

% X=data(:,2:3);

frame=tab(:,1);
frame=frame-min(frame);           % if it does not start with 0
% frame=frame/dt;                 % frame in sec
N=(max(frame)-min(frame)+1);      % number of frames

msd=zeros(1,N);
Time = [0:N-1]*dt;
f=zeros(1,N);


% Find the frames that have been recorded

for j=0:N
    if ~isempty(find(frame==j))
        f(j+1)=find(frame==j,1);
    end
end



for i=1:N-1
    c=0;
    
    for j=0:N-i-1
        if f(i+j+1)>0 && f(j+1)>0
            c=c+1;
            msd(i+1)=msd(i+1)+ sum((X(f(i+j+1),:)-X(f(j+1),:)).^2 );
        end
    end
    if c>0
        msd(i+1)=msd(i+1)/c;
    end
end


msd=msd*dx^2/dt;



subplot(2,3,2)
plot(Time,msd,'-r');hold on;
title('Mean square displacement');
xlabel('time step','FontSize',12);
ylabel('MSD (\mu m^2)','FontSize',12);
% axis([0 50 0 5]);
% legend('Nathanael');
% set(leg2,'FontSize',12);
hold on;

% clear dt;

%% My MSD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%transform time in frames

frame=[];
frame=pos(:,3); 
% frame=frame/dt;
% frame=frame-min(frame);
% frame(1,1)=1;

% i = frame --> Reihe
% j = gap; --> Spalte

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate MSD

msd=zeros(max(frame), 50);
msd2=zeros(50, 3);

for i=1:max(frame);    % find frame, for all frames
    vx=find(frame == i);
    
    if isempty(vx)==1; % if frame does not exist, skip   
    else
        
    
    for j=1:floor(max(frame)/4); % time gap
    
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

for m=1:floor(max(frame)/4);
    
        msd2(m,1)=m;
        msd2(m,2)=mean(nonzeros(msd(:,m)));
        msd2(m,3)=std(nonzeros(msd(:,m)));
        
end

msd2(:,4)=msd2(:,2)*dx^2/dt;
msd2(:,5)=msd2(:,1)*dt;
msd2(:,6)=msd2(:,3)*dx^2/dt;

subplot(2,3,2)
% errorbar(msd2(:,5),msd2(:,4),msd2(:,6)); hold on;
plot(msd2(:,5),msd2(:,4),'-g'); hold on;
% plot(msd2(:,5), 4*D*msd2(:,5),'-g')



%% Velocity

vel=zeros(length(pos),1);

for k=2:num_steps;
    
    vel(k,1)=sqrt(((pos(k,1)-pos(k-1,1))^2)+((pos(k,2)-pos(k-1,2))^2));

    
end

vel=vel*dx/dt; % vel in ?m/s
subplot(2,3,3)
scatter(pos(:,3),vel,3,pos(:,3));
title('Velocity');
xlabel('time (s)','FontSize',12);
ylabel('velocity (\mu m/s)','FontSize',12);


subplot(2,3,6)
hist(vel,20)
title('Histogram of Velocity');
xlabel('time (s)','FontSize',12);
ylabel('velocity (\mu m/s)','FontSize',12);

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
h=gcf;
set(h,'PaperOrientation','landscape');

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
h=gcf;
set(h,'PaperOrientation','landscape');

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



