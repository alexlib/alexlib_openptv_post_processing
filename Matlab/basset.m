function [F,du,y,tmp1] = basset(t,uRel,T)
% basset(urel) estimates the basset force using
% 6 points along the trajectory
%


% T = 6; 

dt = 1; % 1/160;
tau = T:-1:1;
k_v = 1;
density_l = 1;
r = 1;

% t = 1:length(uRel);

F = zeros(length(uRel),1);


figure
hold on

% basset force loop
for ind = T:length(t)
    
    % urel = [data_velocity_p_f{kFiles,iEvents,ind-5:ind}]; % 6 x [U,V,W]
    urel = uRel(ind-T+1:ind);
    
    for k = 1:1 % 3
        % dt in seconds, urel in m/sec
        % du = gradient(urel(k:3:end),dt);
        du = gradient(urel(k:1:end),dt);
        
        y = du./ sqrt(pi * k_v * tau * dt);
        
        tmp1 =  6*pi * r^2 * k_v * density_l * dt * trapz(y);
        
        % basset_force1{kFiles,iEvents,ind}(k) = tmp1;
        F(ind) = tmp1; 
        
        % plot(t(ind-T+1:ind),urel,':',t(ind-T+1:ind),du,'-',t(ind-T+1:ind),y,'--',t(ind),F(ind),'ro');
        
    end
    
end % basset