% function res=RQ(first,last,tail,stepSize);
%RQ(1,6000,0,200)


% plot_cos_s_S



clear all

files = {'0104_water','0104_hompol'};
diss = {2.2e-5, 1.9e-5}; % m^2/s^3
nu = 1e-6; % m^2/s

urms = {0.0075,0.0069};
% L = {0.0194, 0.015};  % these are chosen to fit epsilon and u^3
L = {0.02,0.02};



styles = {'-','--'};


nBins = 101;


for ifile = 1:2
    
    
    
    load(files{ifile})
    
    [s11,s12,s13,...
        s22,s23,s33,ind] = select99percent(s11,s12,s13,s22,s23,s33);
    
    
    
    % -------------------------------------------------------------------------
    s = [s11,s12,s13,...
        s12,s22,s23,...
        s13,s23,s33];
    
    
    
    numPoints = length(s11);
    
    
    mean_s = (sum(s(:).^2)/numPoints).^(0.5);
    
    
    w1 = w1(ind);
    w2 = w2(ind);
    w3 = w3(ind);
    
    
    
    ux = s11;
    uy = s12-0.5*w3;
    uz = s13+0.5*w2;
    vx = s12+0.5*w3;
    vy = s22;
    vz = s23-0.5*w1;
    wx = s13-0.5*w2;
    wy = s23+0.5*w1;
    wz = s33;
    
    [R,Q] = deal(ux);
    
    for ii = 1:numPoints
        
        uM = [ux(ii), uy(ii), uz(ii);...
            vx(ii), vy(ii), vz(ii);...
            wx(ii) wy(ii) wz(ii)];
        
        Q(ii) = -(1/2)*(trace(uM^2))/mean_s^2;
        R(ii) = -(1/3)*(trace(uM^3))/mean_s^3;
    end
    
    
    jointPDF(R,Q ,nBins,nBins,'R','Q')
    
    % Quiver map
    
    ends = find(diff(age(ind)) < 0 );
    starts = [0;ends(1:end-1)+1];
    nTraj = length(starts);
    
    u = diff([R;0]);
    v = diff([Q;0]);
    
    minLength = 120;
    
    figure, hold on;
    
    for i = 1:nTraj-1 % for all trajectories
        
        ii = starts(i):ends(i);
        
        if length(ii) > minLength
            quiver(R(ii(3:end-2)),Q(ii(3:end-2)),u(ii(3:end-2)),v(ii(3:end-2)),1);
        end
    end
    
    hold off
end









