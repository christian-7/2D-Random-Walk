%% Perform least square fitting to determine D_conf%%%%%%%%%%%%%%%%%%%%%
%
% script performs fitting of confinement model to MSD vs time curves
% MSD curves are taken from structure allMSD
%
% info on fitting: http://www.walkingrandomly.com/?p=5196
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Position',[500 800 600 600], 'name','MSD per well analysis using conf model')

n=15; % number of MSD curves

for r=1:n;

    x=allMSD.MSD{r,1}(1:10,1);
    y=allMSD.MSD{r,1}(1:10,2);
    
xfit=0:0.1:5.5;

%%%%%% with all 4 parameters free %%%%%%

% fun = @(p,x) p(1)^2*(1-p(2)*exp((-4*p(3)*p(4)*x)/p(1)^2));
% pguess = [0.05,2,5,0.001]; % r, A1, A2, D

%%%%%% with only 2 free parameters %%%%%%

fun = @(p,x) p(1)^2*(1-1*exp((-4*2*p(2)*x)/p(1)^2));
pguess = [0.1,0.001]; % r, D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[p,fminres] = lsqcurvefit(fun,pguess,x,y)


allMSD.fit{1,1}(r,1)=p(1); % radius
allMSD.fit{1,1}(r,2)=p(2); % A1, D
% allMSD.fit{1,1}(r,3)=p(3); % A2
% allMSD.fit{1,1}(r,4)=p(4); % D
allMSD.fit{1,1}(r,5)=fminres; % minres

%%%%%% Plot MSD curves and overlay fitting %%%%%%%%%%%%%%% 

subplot(2,2,1)
% plot(xfit,p(1)^2*(1-p(2)*exp((-4*p(3)*p(4)*xfit)/p(1)^2)));hold on;
plot(xfit,p(1)^2*(1-1*exp((-4*2*p(2)*xfit)/p(1)^2)));hold on;
scatter(x,y);
xlabel('lag time (s)','FontSize',12);
ylabel('D_{conf} (\mum^2/s)','FontSize',12);

clear x y;
end

%%%%%%%%%%%%%%%%%%%%Plot the results%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,2)
scatter(cell2mat(allMSD.MSD(:,2)),allMSD.fit{1,1}(:,2))
xlabel('dwell time (s)','FontSize',12);
ylabel('D_{conf} (\mum^2/s)','FontSize',12);

subplot(2,2,3)
scatter(allMSD.fit{1,1}(:,1),allMSD.fit{1,1}(:,2))
xlabel('radius (\mum)','FontSize',12);
ylabel('D_{conf} (\mum^2/s)','FontSize',12);

subplot(2,2,4)
scatter(cell2mat(allMSD.MSD(:,2)),allMSD.fit{1,1}(:,1))
xlabel('dwell time (s)','FontSize',12);
ylabel('radius (\mum)','FontSize',12);

%% Perform linear fitting to determine D_conf%%%%%%%%%%%%%%%%%%%%%
%
% script performs fitting of linear model to MSD vs time curves
% MSD curves are taken from structure allMSD
%
% info on fitting: http://www.walkingrandomly.com/?p=5196
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

figure('Position',[500 800 600 300], 'name','MSD per well analysis using linear regression')

n=15;

for r=1:n;

a1 = polyfit(allMSD.MSD{r,1}(2:10,1),allMSD.MSD{r,1}(2:10,2),1);

allMSD.lin_fit{r,1}(:,1)=a1(:,1);
allMSD.lin_fit{r,1}(:,2)=a1(:,2);

xfit=0:0.1:2.5;

subplot(1,2,1);
scatter(allMSD.MSD{r,1}(1:10,1),allMSD.MSD{r,1}(1:10,2),'*'); hold on;
plot(xfit,a1(1)*xfit+a1(2));hold on;
legendInfo{r} = ['MSD well = ' num2str(r)]; 

subplot(1,2,2)
scatter(cell2mat(allMSD.MSD(r,2)),allMSD.lin_fit{r,1}(:,1));hold on;
xlabel('dwell time (s)','FontSize',12);
ylabel('D_{conf} (\mum^2/s)','FontSize',12)

end

%% %% Perform least square fitting using jacobian matrix to determine D_conf and radius %%%%%%%%%%%%%%%%%%%%%
%
% script performs fitting of confinement model to MSD vs time curves
% MSD curves are taken from structure allMSD
%
% fitting function in conf_model_jacobian_fit.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc
load('/Users/Christian/Documents/Arbeit/Manuscripts/Receptor_organization/trajectories/data/well_analysis/allMSD.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par0=[0.2,1,2,0.005]; % starting parameters
lowerBound=[0.1,0.01,0.01,0.001];
upperBound=[1,10,10,0.1];
xfit=0:0.1:5.5;

figure('Position',[500 800 600 300], 'name','MSD per well analysis using linear regression')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r=1:15;
    
    if r==10 | r==11;
        
    else
                
x=allMSD.MSD{r,1}(1:10,1);
y=allMSD.MSD{r,1}(1:10,2);

[fitParam,resnorm,residual,exitflag]=conf_model_jacobian_fit(par0,lowerBound,upperBound,x,y);

fit_result(:,r)=fitParam;
residual2(:,r)=residual;

% 1 - radius
% 2 - A1
% 3 - A2
% 4 - D


% subplot(2,2,1)
% plot(xfit,fit_result(1,r)^2*(1-fit_result(2,r)*exp((-4*fit_result(3,r)*fit_result(4,r)*xfit)/fit_result(1,r)^2)));hold on;
% scatter(x,y);hold on;
% xlabel('lag time (s)','FontSize',12);
% ylabel('MSD (\mum^2/s)','FontSize',12);
% % legendInfo{r} = ['MSD well = ' num2str(r)]; 
% % legend(legendInfo);
% box on;
% 
% subplot(2,2,3)
% plot(residual2(:,r));hold on;
% ylabel('residuals','FontSize',12);
% xlabel('lag time (s)','FontSize',12);
% axis([1 10 -0.2 0.2])
% box on;

% subplot(2,2,[2,4])
scatter(cell2mat(allMSD.MSD(r,2)),fit_result(4,r));hold on;
xlabel('dwell time (s)','FontSize',12);
ylabel('D_{conf} (\mum^2/s)','FontSize',12)
box on;

clear fitParam resnorm residual exitflag x y

    end
end

% figure
% scatter(cell2mat(allMSD.MSD(:,2)),fit_result(1,:))
% xlabel('dwell time (s)','FontSize',12);
% ylabel('radius (\mum)','FontSize',12)