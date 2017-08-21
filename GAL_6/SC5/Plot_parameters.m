clear all
close all
%figure
% Output file columns content: 
%
% Column 1  =  Number of single stars + binaries        
% Column 2  =  Number of binaries                       
% Column 3  =  Time (Myr)                               
% Column 4  =  Relaxation time (Myr)                    
% Column 5  =  Total cluster mass (Msun)                
% Column 6  =  Mass in core (Msun)                      
% Column 7  =  Mass outside the tidal radius (Msun)     
% Column 8  =  Maximum stellar distance from cluster centre of mass (pc)        
% Column 9  =  Half-mass radius (pc)
% Column 10 =  Radius containing inner 10% of cluster mass (pc)                 
% Column 11 =  Core radius - as determined by Nbody code (pc)                   
% Column 12 =  Number of systems (stars + binaries) inside the half-mass radius  
% Column 13 =  Number of systems within 1pc of the cluster centre               
% Column 14 =  Number of systems within the inner lagrangian radius (10%)       
% Column 15 =  Number of systems within the core radius
% Column 16 =  Velocity dispersion (km/s) 


% Ectract the data from the input file

data = load('extrct.dat'); % Ectract the data from the input file

N_sin_bin = data(:,1); 

N_bin = data(:,2);

t = data(:,3); 

t_relax = data(:,4);

M_tot = data(:,5);

M_core = data(:,6);

M_out_tidal = data(:,7);

D_max = data(:,8);

R_hm = data(:,9);

R_10 = data(:,10);

R_core = data(:,11);

N_hm = data(:,12);

N_1pc = data(:,13);

N_lagr = data(:,14);

N_core = data(:,15);

sigma_v = data(:,16);




% Plot of the results

subplot(1,3,1), plot(t, M_tot,'k')
xlabel('t (Myr)','fontsize',13)
ylabel('M_{tot} (M_{o})','fontsize',13)
set(gca,'fontsize',10)
axis([0 max(t) 0 8e3])
axis square



subplot(1,3,2), plot(t, R_hm,'k')
xlabel('t (Myr)','fontsize',13)
ylabel('R_{hm} (pc)','fontsize',13)
set(gca,'fontsize',10)
axis([0 max(t) 0 10])
axis square


subplot(1,3,3), plot(t, sigma_v,'k')
xlabel('t (Myr)','fontsize',13)
ylabel('\sigma_v (km s^{-1})','fontsize',13)
set(gca,'fontsize',10)
axis([0 max(t) 0 6])
axis square



