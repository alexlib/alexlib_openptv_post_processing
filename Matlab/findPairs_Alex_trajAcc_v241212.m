% 25-Dec-2012
% - File is modified to load the data from the MAT file
% - trajACC is now a separate function in the same directory
% - looking for pairs rather than for folding 
% - project with Anna Frishman and Grisha Falkovich
% 
% ----
% 08-Mar-2006 10:43:15
% version recompiled for reading trajAcc files directly, not HDF
% in order to compare with C++ from Beat
% the comparison is on the F:\PTV\Nov2004\res_02Novw\ data
% C++ is running from: D:\My Documents\ETH\People\Beat\folding\c++\
% debugged and fixed on my VC 6++
% - start using nested functions, in order to include all the codes in one
%
%
%
% 01-Mar-2006 22:43:49
% Version 1.0
% - change the filtered derivatives to the sphere around the center of the
% pair
% - added try/catch and fclose(fid)
%

function findPairs_Alex_trajAcc_v241212


% Some service functions
dist = @(ax,ay,az,bx,by,bz)( (ax - bx).^2 + (ay - by).^2 + (az -bz ).^2 ).^.5;

% formatString = ['''',repmat('%6.3f ',1,34),'\n'''];

% Globals
fromR = 0.0; % mm, the inner "pipe" radius
toR =   2.2;   % mm, the outer pipe radius
maxR = 100; % mm

minCloseEnough = 10; % minimum 45 frames together, 3 \tau_eta
% coarseR = 2; % size of the coarse-grained filter, maybe it's adaptive to the local size, check.

tau_eta = 0.23; % sec, Kolmogorov time, mai10
dT = 1/60; % sec, time separation between two frames

outfile = sprintf('mai10_pairs_%dto%dmm_1taueta.txt',fromR,toR);
outmatfile = sprintf('mai10_pairs_%dto%dmm_1taueta.mat',fromR,toR);

MINLENGTH = 10; % this is just a flag, go into readTrajAcc and change THERE
% [traj] = readTrajAcc('../datasets/april26_trajAcc_Beat',10000,16035);
% load('april26_trajAcc_Beat'); 


% for coarse-graining
% t = cat(2,traj.t).';
% 
% x = cat(2,traj.x).';
% y = cat(2,traj.y).';
% z = cat(2,traj.z).';

% ux = cat(2,traj.ux).';
% uy = cat(2,traj.uy).';
% uz = cat(2,traj.uz).';
% vx = cat(2,traj.vx).';
% vy = cat(2,traj.vy).';
% vz = cat(2,traj.vz).';
% wx = cat(2,traj.wx).';
% wy = cat(2,traj.wy).';
% wz = cat(2,traj.wz).';

numTraj = length(traj);
% maxLen    = max(cat(2,traj.trajLen)).';
maxLen = 200;


% Main loop
fid = fopen(outfile,'w');
pairCounter = 0;
pair = zeros(100000,2);

try
    for i = 1:numTraj-1 % through all spagetties
        for j = i+1:min(i+maxLen,numTraj) %through all OTHER spagetties
            
            % Remove the pair from the list
            %         if ismember([i,j],pair(1:pairCounter),'rows') || ismember([j,i],pair(1:pairCounter),'rows')
            %             continue
            %         end
            
            %                 fprintf(1,'\b\b\b\b\b\b\b\b\b');
            %                 fprintf(1,'%9d',pairCounter);
            
            [pairT,iA,iB] = intersect(traj(i).t,traj(j).t);
            
            if ~isempty(pairT)
                
                rPair = dist(...
                    traj(i).x(iA),traj(i).y(iA),traj(i).z(iA),...
                    traj(j).x(iB),traj(j).y(iB),traj(j).z(iB));
                
                %             iniClose    = find(rPair <= iniR,1,'first');
                iniClose    = find(rPair <= toR & rPair > fromR,1,'first'); % OR LAST? ????
                stillClose  = min(iA(end),iB(end)); % find(rPair <= maxR,1,'last');
                
                if ~isempty(iniClose)  && (stillClose - iniClose) >= minCloseEnough
                    
                    pairCounter = pairCounter + 1;
                    pair(pairCounter,1:2) = [i,j];
                    
                    k = iniClose:stillClose;
                    
                    disp(sprintf('%d : %d + %d, %d frames \n',pairCounter,i,j,length(k)))
                    
                    buffer = zeros(length(k),7);
                    
                    buffer(:,1) = (traj(i).t(iA(k)) - traj(i).t(iA(k(1))))/tau_eta*dT; % first column, time
                    
                    %                 buffer(:,2) = traj(i).x(iA(k)); % x,y,z
                    %                 buffer(:,3) = traj(i).y(iA(k));
                    %                 buffer(:,4) = traj(i).z(iA(k));
                    
                    buffer(:,2) = (- traj(i).x(iA(k)) + traj(j).x(iB(k))); % l_1, l_2, l_3 = l_x,l_y, l_z
                    buffer(:,3) = (- traj(i).y(iA(k)) + traj(j).y(iB(k)));
                    buffer(:,4) = (- traj(i).z(iA(k)) + traj(j).z(iB(k)));
                    
                    buffer(:,5) = (- traj(i).u(iA(k)) + traj(j).u(iB(k))); % l_1, l_2, l_3 = l_x,l_y, l_z
                    buffer(:,6) = (- traj(i).y(iA(k)) + traj(j).v(iB(k)));
                    buffer(:,7) = (- traj(i).w(iA(k)) + traj(j).w(iB(k)));
                    
                    fprintf(fid,'%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f  \n',buffer');
                end % if ~isempty
            end % if ~isempty(iniClose)
        end % for j
    end % for i
    
catch
    lasterr
    fclose(fid)
end

fclose(fid)
pair = pair(1:pairCounter,:);
save(outmatfile,'pair');
end


% ----------------------------------------------------

function testResults
% Test results:
tmp = load('pairs.txt');
startPoints = find(tmp(:,1) == 0);
numPairs = length(startPoints);
lenPairs = diff([startPoints(1:end);length(tmp)]);

pairs = repmat(struct('t',[],'lx',[],'ly',[],'lz',[],'du',[],'dv',[],'dw',[]),[numPairs,1]);


for i = 1:numPairs
    k = startPoints(i);
    k = k:k+lenPairs(i)-1;
    pairs(i).t = tmp(k,1);
    pairs(i).lx = tmp(k,2);
    pairs(i).ly = tmp(k,3);
    pairs(i).lz = tmp(k,4);
    pairs(i).du = tmp(k,5);
    pairs(i).dv = tmp(k,6);
    pairs(i).dw = tmp(k,7);
end

r0 = zeros(length(pairs),1);
for i = 1:length(pairs); r0(i) = (pairs(i).lx(1).^2 + pairs(i).ly(1).^2 + pairs(i).lz(1).^2)^.5; end

figure, nhist(r0,100)



r = {};
du = {};
for i = 1:numPairs, 
    r{i} = (pairs(i).lx.^2 + pairs(i).ly.^2 + pairs(i).lz.^2).^.5; 
    du{i} = (pairs(i).du.^2 + pairs(i).dv.^2 + pairs(i).dw.^2).^.5; 
end

figure, hold on;
for i = 1:10:numPairs, 
    scatter(log10(du{i}),log10(r{i}),'.'); 
end
hold off


figure, hold on;
for i = 1:10:numPairs, 
    plot(pairs(i).t,r{i},'DisplayName',i); 
end
hold off

% set(gca,'yscale','log','xscale','log')


% ind = find(rand(numPairs,1) > .99);
% figure, hold on
% for k = 1:length(ind), i = ind(k);
%     quiver3(pair(i).x,pair(i).y,pair(i).z,pair(i).l1,pair(i).l2,pair(i).l3,'b')
% end
% 
% end

end




