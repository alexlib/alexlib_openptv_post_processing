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

function findPairs_Alex_1_trajAcc


% Some service functions
dist = @(ax,ay,az,bx,by,bz)( (ax - bx).^2 + (ay - by).^2 + (az -bz ).^2 ).^.5;

% formatString = ['''',repmat('%6.3f ',1,34),'\n'''];

% Globals
fromR = 1; % mm, the inner "pipe" radius
toR = 2;   % mm, the outer pipe radius
maxR = 30; % mm 

minCloseEnough = 45; % minimum 45 frames together, 3 \tau_eta
% coarseR = 2; % size of the coarse-grained filter, maybe it's adaptive to the local size, check.

tau_eta = 0.23; % sec, Kolmogorov time, mai10
dT = 1/60; % sec, time separation between two frames

outfile = sprintf('mai10_pairs_%dto%dmm_3taueta.txt',fromR,toR);


MINLENGTH = 20; % this is just a flag, go into readTrajAcc and change THERE
[traj] = readTrajAcc('test',10000,10500);


% for coarse-graining
t = cat(2,traj.t).';

x = cat(2,traj.x).';
y = cat(2,traj.y).';
z = cat(2,traj.z).';

ux = cat(2,traj.ux).';
uy = cat(2,traj.uy).';
uz = cat(2,traj.uz).';
vx = cat(2,traj.vx).';
vy = cat(2,traj.vy).';
vz = cat(2,traj.vz).';
wx = cat(2,traj.wx).';
wy = cat(2,traj.wy).';
wz = cat(2,traj.wz).';

numTraj = length(traj);
% maxLen    = max(cat(2,traj.trajLen)).';
maxLen = 200;


% Main loop
fid = fopen(outfile,'w');
pairCounter = 0;
pair = zeros(100,2);

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

            rPair = dist(...
                traj(i).x(iA),traj(i).y(iA),traj(i).z(iA),...
                traj(j).x(iB),traj(j).y(iB),traj(j).z(iB));

            %             iniClose    = find(rPair <= iniR,1,'first');
            iniClose    = find(rPair <= toR & rPair > fromR,1,'first');
            stillClose  = find(rPair <= maxR,1,'last');

            if ~isempty(iniClose) && (stillClose - iniClose) >= minCloseEnough

                pairCounter = pairCounter + 1;
                pair(pairCounter,1:2) = [i,j];

                k = iniClose:stillClose;

                disp(sprintf('%d : %d + %d, %d frames \n',pairCounter,i,j,length(k)))

                buffer = zeros(length(k),34);

                buffer(:,1) = (traj(i).t(iA(k)) - traj(i).t(iA(k(1))))/tau_eta*dT; % first column, time

                buffer(:,2) = traj(i).x(iA(k)); % x,y,z
                buffer(:,3) = traj(i).y(iA(k));
                buffer(:,4) = traj(i).z(iA(k));

                buffer(:,5) = (- traj(i).x(iA(k)) + traj(j).x(iB(k))); % l_1, l_2, l_3 = l_x,l_y, l_z
                buffer(:,6) = (- traj(i).y(iA(k)) + traj(j).y(iB(k)));
                buffer(:,7) = (- traj(i).z(iA(k)) + traj(j).z(iB(k)));

                buffer(:,8)  = traj(i).ux(iA(k)); % derivatives, dudx
                buffer(:,9)  = traj(i).uy(iA(k));     % dudy
                buffer(:,10) = traj(i).uz(iA(k));    % dudz
                buffer(:,11) = traj(i).vx(iA(k));
                buffer(:,12) = traj(i).vy(iA(k));
                buffer(:,13) = traj(i).vz(iA(k));
                buffer(:,14) = traj(i).wx(iA(k));
                buffer(:,15) = traj(i).wy(iA(k));
                buffer(:,16) = traj(i).wz(iA(k));

                buffer(:,17) = traj(j).ux(iB(k));                           % dudx
                buffer(:,18) = traj(j).uy(iB(k));   % dudy
                buffer(:,19) = traj(j).uz(iB(k));    % dudz
                buffer(:,20) = traj(j).vx(iB(k));    % dvdx
                buffer(:,21) = traj(j).vy(iB(k));                            % dvdy
                buffer(:,22) = traj(j).vz(iB(k));    % dvdz
                buffer(:,23) = traj(j).wx(iB(k));    % dwdx
                buffer(:,24) = traj(j).wy(iB(k));    % dwxy
                buffer(:,25) = traj(j).wz(iB(k));                            % dwdz


                % Additional loop, to find all the neighbours and get
                % "coarse-grained" derivatives

                for m = 1:length(k)

                    ind = find(t == traj(i).t(iA(k(m))));

                    cgx  = (traj(i).x(iA(k(m))) + traj(j).x(iB(k(m))))/2; % center of gravity
                    cgy =  (traj(i).y(iA(k(m))) + traj(j).y(iB(k(m))))/2;
                    cgz =  (traj(i).z(iA(k(m))) + traj(j).z(iB(k(m))))/2;



                    distAll = sqrt(...
                        (x(ind) - cgx).^2 + ...
                        (y(ind) - cgy).^2 + ...
                        (z(ind) - cgz).^2);

                    localR = sqrt(sum(buffer(m,5:7).^2,2));

                    n = find(distAll <= localR);
                    %                     disp(sprintf('%d neighbours',n))

                    buffer(m,26) = mean(ux(ind(n)));
                    buffer(m,27) = mean(uy(ind(n)));
                    buffer(m,28) = mean(uz(ind(n)));
                    buffer(m,29) = mean(vx(ind(n)));
                    buffer(m,30) = mean(vy(ind(n)));
                    buffer(m,31) = mean(vz(ind(n)));
                    buffer(m,32) = mean(wx(ind(n)));
                    buffer(m,33) = mean(wy(ind(n)));
                    buffer(m,34) = mean(wz(ind(n)));
                end % for m
                fprintf(fid,'%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f \n',buffer');
            end % if
        end % for j
    end % for i

catch
    lasterr
    fclose(fid)
end

fclose(fid)

end


% ----------------------------------------------------

function testResults
% Test results:
pairs = load('pairs.txt');
startPoints = find(pairs(:,1) == 0);
numPairs = length(startPoints);
lenPairs = diff([startPoints(1:end);length(pairs)]);

pair = repmat(struct('t',[],'x',[],'y',[],'z',[],'l1',[],'l2',[],'l3',[]),numPairs);


for i = 1:numPairs
    k = startPoints(i);
    k = k:k+lenPairs(i)-1;
    pair(i).t = pairs(k,1);
    pair(i).x = pairs(k,2);
    pair(i).y = pairs(k,3);
    pair(i).z = pairs(k,4);
    pair(i).l1 = pairs(k,5);
    pair(i).l2 = pairs(k,6);
    pair(i).l3 = pairs(k,7);
end


ind = find(rand(numPairs,1) > .99);
figure, hold on
for k = 1:length(ind), i = ind(k);
    quiver3(pair(i).x,pair(i).y,pair(i).z,pair(i).l1,pair(i).l2,pair(i).l3,'b')
end

end

% ----------------------------------------------------
function [traj] = readTrajAcc(varargin)

MINLENGTH = 20;

traj = struct('x',[],'y',[],'z',[],...
    ...'u',[],'v',[],'w',[],...
    'ux',[],'uy',[],'uz',[],...
    'vx',[],'vy',[],'vz',[],...
    'wx',[],'wy',[],'wz',[],...
    ...'dudt',[],'dvdt',[],'dwdt',[],...
    ...'ax',[],'ay',[],'az',[],...
    ...'daxdx',[],'daxdy',[],'daxdz',[],...
    ...'daydx',[],'daydy',[],'daydz',[],...
    ...'dazdx',[],'dazdy',[],'dazdz',[],...
    'trajnum',0,'t',[]);

fieldNames = fieldnames(traj);

if nargin == 0 % no inputs

[filename1, pathname] = uigetfile('trajAcc.*', 'Pick a FIRST trajectory file');

if isequal(filename1,0) || isequal(pathname,0)
    error('Wrong selection')
end

wd = cd;
cd(pathname);

[filename2, pathname] = uigetfile('trajAcc.*', 'Pick a LAST trajectory file');
if isequal(filename2,0) || isequal(pathname,0)
    cd(wd); error('Wrong selection')
end
    

% First and last index of the files
firstFileIndx =  eval(filename1(findstr(filename1,'.') + 1:end));
lastFileIndx =  eval(filename2(findstr(filename2,'.') + 1:end));

cd(wd);

else
    pathname = varargin{1};
    firstFileIndx = varargin{2};
    lastFileIndx = varargin{3};
end

fprintf(1,'Wait please ...  ')
disp('');
trajIndx = 1;

% Loop through all files
for n = firstFileIndx:lastFileIndx

    runIndx = n - firstFileIndx + 1; % running counter

    fprintf(1,'\b\b\b\b\b\b\b\b\b');
    fprintf(1,'%9d',n);

    % f = load(sprintf('%s%d',baseName,n),'ASCII');
    fid = fopen(sprintf('%s%s%s%d',pathname,filesep,'trajAcc.',n));
    if fid == -1, error('no file'), continue, end
    % f = fscanf(fid,'%g',[35 inf]); % It has two rows now.
    f = fscanf(fid,'%g',[32 inf]); % It has two rows now.
    if isempty(f), continue , end; % if the file is empty, skip to next file
    f = f(1:32,:)';
    fclose(fid);



    % prepare indices of the single trajectories
    [startTraj, endTraj, lenTraj] = singleTraj(f(:,32));

    % pick only long trajectories, longer than ' MINLENGTH '
    longTraj = find(lenTraj >= MINLENGTH);
    if isempty(longTraj), fprintf(1,'!'), continue, end; % if there are no long trajectories, move to the next file


    numLongTraj = length(longTraj);     % number of long trajectories
    startTraj = startTraj(longTraj);    % starting
    endTraj = endTraj(longTraj);        % ending points
    lenTraj = lenTraj(longTraj);        % lengthes
    %     maxLenTraj = max(lenTraj);          % maximum length

    for i = 1:numLongTraj
        
        for q = 1:3 % only first 3, x,y,z
            traj(trajIndx).(fieldNames{q})(1:lenTraj(i)) = ...
                1000 * f(startTraj(i):endTraj(i),q)'; % 1000 for mm
        end
        
        for q = 4:12 % only dudx
            traj(trajIndx).(fieldNames{q})(1:lenTraj(i)) = ...
                f(startTraj(i):endTraj(i),q)';
        end
        % trajIndex
        traj(trajIndx).(fieldNames{13})(1:lenTraj(i)) = ...
            f(startTraj(i):endTraj(i),31)';
        % running time
        traj(trajIndx).(fieldNames{14})(1:lenTraj(i)) = ...
            f(startTraj(i):endTraj(i),32)'+runIndx;
        trajIndx = trajIndx + 1;
    end

end % FILES
end % of function

% ---------------------------------------------------------------------------------------------------------
function [startTraj, endTraj, lenTraj] = singleTraj(age)
% Split the file into single trajectories, according to the
% AGE vector, column 29 of the trajPoint.* files
startTraj = find(age == 0);      % start points of the trajectories
endTraj = startTraj(2:end)-1;    % end points of the trajectories
endTraj = [endTraj;length(age)];
lenTraj = endTraj - startTraj + 1;
%    numTraj = length(startTraj);
end % of function




