figure, hold on,

hfiles{1}  = ('/Volumes/WD Disk/My_Documents/MATLAB/Alex/HDF/0104_water.hdf');
hfiles{2}  = ('/Volumes/WD Disk/My_Documents/MATLAB/Alex/HDF/0104_hompol.hdf');
% hfiles{3}  = ('/Volumes/WD Disk/My_Documents/MATLAB/Alex/HDF/2604_water.hdf');
% hfiles{4}  = ('/Volumes/WD Disk/My_Documents/MATLAB/Alex/HDF/2604_pol.hdf');

for fileNum = 1:2
    
    hfile = hfiles{fileNum}
    
    [pathstr,namestr] = fileparts(hfile);
    
    [traj,attr] = readTrajHDF_v9(hfile,...
        'x',[],'y',[],'z',[],...
        'u',[],'v',[],'w',[],...
        's11',[],'s12',[],'s13',[],...
        's22',[],'s23',[],'s33',[],...
        'w1',[],'w2',[],'w3',[],...
        'minlength',20); % ,'frames',[1 100]);
    
    numTraj = length(traj)
    numPoints = length(cat(1,traj.x))
    trajLen = cat(1,traj.trajLen);
    
    for i = 1:numTraj
        traj(i).dudx = traj(i).s11;
        traj(i).dudy = traj(i).s12 - .5*traj(i).w3;
        traj(i).dudz = traj(i).s13 + .5*traj(i).w2;
        traj(i).dvdx = traj(i).s12 + .5*traj(i).w3;
        traj(i).dvdy = traj(i).s22;
        traj(i).dvdz = traj(i).s23 - .5*traj(i).w1;
        traj(i).dwdx = traj(i).s13 - .5*traj(i).w2;
        traj(i).dwdy = traj(i).s23 + .5*traj(i).w1;
        traj(i).dwdz = traj(i).s33;
    end
    
    k = 0;
    good = zeros(numTraj,1);
    fn = fieldnames(traj); fn1 =   find(strcmp(fn,'trajLen'));
    for i=1:numTraj
        absdiv = abs(traj(i).dudx +  traj(i).dvdy + traj(i).dwdz);
        reldiv = absdiv./(abs(traj(i).dudx) +  abs(traj(i).dvdy) + abs(traj(i).dwdz));
        %    ind = find(reldiv < .1); %  & absdiv < .1);
        ind = find(reldiv < .1 & absdiv < .1);
        ind1 = diff([ind;Inf]);
        ind2 = [0; find(ind1 > 1)];
        ind3 = diff(ind2);
        [r,c] = max(ind3);
        ind = ind(ind2(c)+1:ind2(c+1));
        if length(ind) > 15
            k = k + 1;
            good(k) = i;
            for j = [1:fn1-1,fn1+1:length(fn)]
                %             if length(traj(i).(fn{j}) > length(ind))
                traj(i).(fn{j}) = traj(i).(fn{j})(ind);
                %             end
            end
            traj(i).trajLen = length(ind);
        end
        
    end
    good = good(1:k);
    traj = traj(good);
    numTraj = length(traj)
    numPoints = length(cat(1,traj.dudx))
    trajLen = cat(1,traj.trajLen);
    trajnum = cat(1,traj.trajnum);
    
    
    dudx = cat(1,traj.dudx);
    dudy = cat(1,traj.dudy);
    dudz = cat(1,traj.dudz);
    dvdx = cat(1,traj.dvdx);
    dvdy = cat(1,traj.dvdy);
    dvdz = cat(1,traj.dvdz);
    dwdx = cat(1,traj.dwdx);
    dwdy = cat(1,traj.dwdy);
    dwdz = cat(1,traj.dwdz);
    
    maxLen    = max(trajLen);
    minLen    = min(trajLen);
    meanLen    = mean(trajLen);
    
    u = cat(1,traj.u);
    v = cat(1,traj.v);
    w = cat(1,traj.w);
    omega = [dwdy - dvdz, dudz - dwdx, dvdx - dudy];
    
    vel = [u,v,w];
    
    C = cross(omega,vel);
    
%     omegarms = mean((omega(:,1).^2 + omega(:,2).^2 + omega(:,3).^2)).^(.5);
%     velrms = mean((vel(:,1).^2 + vel(:,2).^2 + vel(:,3).^2)).^(.5);
    Crms = (C(:,1).^2 + C(:,2).^2 + C(:,3).^2).^(0.5);
%     nhist(Crms/omegarms/velrms);
nhist(Crms);
    
    
end
