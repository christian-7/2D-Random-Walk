clear all, close all, clc

dwell_time3=struct('dwell', []);


for y=4:20; % loop for ratios

    dwell_time2=[];
    
    
for qw=1:5; % number of iterations
    
    ratio=2*y;
    
    random_walk_w_confinement_fixed_wells_2

     
        if      isempty(dwell_time)==1;
       
        else
        
            dwell_time2=cat(1,dwell_time2,dwell_time);
    
        end
    
    dwell_time3.dwell{y,1}=dwell_time2;
    dwell_time3.dwell{y,2}=ratio;
    
    
    clear dwell_time 
    
end

clear dwell_time2 
end


figure('Position',[500 800 400 350], 'name','Dwell time vs. D_conf simulation')

for x=1:y;
    
        %errorbar(mean(dwell_time3.dwell{x, 1}),D/x,std(dwell_time3.dwell{x, 1})); hold on;
        %set(gca,'CameraUpVector',[1,0,0],'YDir','reverse','XAxisLocation','top');
        scatter(mean(dwell_time3.dwell{x, 1}),D/x); hold on;%std(dwell_time3.dwell{x, 1}));hold on;
        xlabel('Dwell time (s)','FontSize',12);
        ylabel('D_{conf} (\mum^2/s)','FontSize',12)
        box on;
  
end

  
  