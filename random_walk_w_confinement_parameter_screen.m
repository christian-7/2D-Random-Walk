% Steps 1 and 2 for a matrix of combinations between radius and dwell time
%
% 1. simulates a random walk according to the initialized parameters and adds confined areas with D_conf
% 2. calculates the confinement index according to Simson, et al. (1995), Biophysical Journal, 69(September), 989?993. 
% 
% 
% dwell time - min radius - step size - max radius
% ratio - min ratio - step size - max ratio
% dt-time step
% 
% Date:     17/05/15
% Author:   Christian Sieben

clear all, close all, clc

%% Initialize parameters for random walk and add circles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D=0.04;                                          % diffusion constant, ?m2/s 
radius=[5:5:100];                                % number of steps the particle stays confined 
ratio=[5:5:100];                                  % ratio D/D_conf

[radius_loop, ratio_loop]=meshgrid(radius, ratio);

radiusdb=radius_loop(:);
ratiodb=ratio_loop(:);
D_conf_db=transpose(D./ratiodb);

start=[0 0];                                % starting coordinates
num_steps=300;                              % number of steps
n=1;                                        % number of cycles per parameter pair
                    
dt=0.5;                                     % time step, seconds
dx=1;                                       % pixel size, ?m
step_size=sqrt(4*D*dt);                     % step size according to D


segment=30;                                 % Sm, segment length in frames

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define confined areas as circles with radius and center position

center=[1 0.5];
center2=[-1 0.5];
center3=[0 -0.6];

minAxis=-2.5;
maxAxis=2.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate and plot random walk trajectory

for index=1:length(radiusdb);
    
    radius=radiusdb(index);
    D_conf=D_conf_db(index);
    step_size_conf=sqrt(4*D_conf*dt);           % step size according to D_conf

pos = zeros(num_steps,2);
pos(1,1)=start(1);
pos(1,2)=start(2);
pos(1,3)=0;

a = 1;      
b = 360;  
angle = (b-a).*rand(num_steps,1) + a;

for k=2:num_steps;
    
    if k>=150 & k<=150+radius %| k>=200 & k<=200+radius
        
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


% Calculate confinement

% generate variable "frame"

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

result(index,1)=max(L(:,1));
result(index,2)=sum(L(:,1));

end
%% 

intmat=reshape(result(:,2),length(radius_loop),length(ratio_loop));
maxmat=reshape(result(:,1),length(radius_loop),length(ratio_loop));

% intmat=cell2mat(intmat);
% maxmat=cell2mat(maxmat);

figure('Position',[100 200 800 600])
h=gcf;
set(h,'PaperOrientation','landscape');

subplot(2,2,1);
surf(radius_loop, ratio_loop, intmat)
% shading interp
view(-42,27)

title('Integral');
xlabel('dwell time (sec)');
ylabel('D/D_conf');
zlabel('Integral');
% colorbar;

subplot(2,2,2);
surf(radius_loop, ratio_loop, maxmat)
% shading interp
view(-42,27)

title('Max');
xlabel('dwell time (sec)');
ylabel('D/D_conf');
zlabel('max L');

subplot(2,2,3);
surf(radius_loop, ratio_loop, intmat)
% shading interp
view(0,90)

title('Integral');
xlabel('dwell time (sec)');
ylabel('D/D_conf');
zlabel('Integral');
colorbar;

subplot(2,2,4);
surf(radius_loop, ratio_loop, maxmat)
% shading interp
view(0,90)
colorbar;

title('Max');
xlabel('dwell time (sec)');
ylabel('D/D_conf');
zlabel('max L');