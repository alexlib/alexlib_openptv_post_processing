
% Run version of the analysisHDF_ver0804.m
% 
hfile = '0211_20.hdf'

runtime = datestr(now)
diary(['analysisHDF_ver0804_',strrep(strrep(runtime,' ','_'),':','_'),'.txt'])

[pathstr,namestr] = fileparts(hfile);

    
[traj,attr] = readTrajAccHDF(hfile,...
    'x',[],'y',[],'z',[],...
    'u',[],'v',[],'w',[],...
    'dudx',[],'dudy',[],'dudz',[],...
    'dvdx',[],'dvdy',[],'dvdz',[],...
    'dwdx',[],'dwdy',[],'dwdz',[],...
    'minlength',20,'frames',[1 3000],'trajnum',[1 Inf]);

numTraj = length(traj)

numPoints = length(cat(1,traj.dudx))
trajLen = cat(1,traj.trajLen);

% 
% 
dudx = cat(1,traj.dudx);
dudy = cat(1,traj.dudy);
dudz = cat(1,traj.dudz);

dvdx = cat(1,traj.dvdx);
dvdy = cat(1,traj.dvdy);
dvdz = cat(1,traj.dvdz);

dwdx = cat(1,traj.dwdx);
dwdy = cat(1,traj.dwdy);
dwdz = cat(1,traj.dwdz);


reldiv = abs(dudx + dvdy + dwdz)./(abs(dudx)+abs(dvdy)+abs(dwdz));
absdiv = abs(dudx + dvdy + dwdz);
inddiv = absdiv < 0.1 & reldiv < 0.1;

k = 0;
good = zeros(numTraj,1);
fn = fieldnames(traj);
for i=1:numTraj
    absdiv = abs(traj(i).dudx +  traj(i).dvdy + traj(i).dwdz);
    reldiv = absdiv./(abs(traj(i).dudx) +  abs(traj(i).dvdy) + abs(traj(i).dwdz));
    ind = find(reldiv < .1 & absdiv < .1);
    ind1 = diff([ind;Inf]);
    ind2 = [0; find(ind1 > 1)];
    ind3 = diff(ind2);
    [r,c] = max(ind3);
    ind = ind(ind2(c)+1:ind2(c+1));
    if length(ind) > 15
        k = k + 1;
        good(k) = i;
        for j = 1:length(fn)- 1 
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


% return



dudx = cat(1,traj.dudx);
dudy = cat(1,traj.dudy);
dudz = cat(1,traj.dudz);

dvdx = cat(1,traj.dvdx);
dvdy = cat(1,traj.dvdy);
dvdz = cat(1,traj.dvdz);

dwdx = cat(1,traj.dwdx);
dwdy = cat(1,traj.dwdy);
dwdz = cat(1,traj.dwdz);



maxLen    = max(trajLen)
minLen    = min(trajLen)
meanLen    = mean(trajLen)

save([namestr,'_good'],'good','numTraj','numPoints','trajnum','trajLen','minLen','runtime');
diary off % diary(['analysisHDF_ver0804_',strrep(strrep(datestr(now),' ','_'),':','_'),'.txt'])