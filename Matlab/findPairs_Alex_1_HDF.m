% 01-Mar-2006 22:43:49
% Version 1.0
% - change the filtered derivatives to the sphere around the center of the
% pair
% - added try/catch and fclose(fid)
%
% Initialization
hfile = '../HDF/mai10.hdf';
addpath ../

% Some service functions
dist = @(ax,ay,az,bx,by,bz)( (ax - bx).^2 + (ay - by).^2 + (az -bz ).^2 ).^.5;

% formatString = ['''',repmat('%6.3f ',1,34),'\n'''];

% Globals
fromR = 1; % mm, the inner "pipe" radius
toR = 2;   % mm, the outer pipe radius
maxR = 30;% mm, the maximum distance, cancelled
minCloseEnough = 45; % minimum 45 frames together, 3 \tau_eta
% coarseR = 2; % size of the coarse-grained filter, maybe it's adaptive to the local size, check.

tau_eta = 0.23; % sec, Kolmogorov time, mai10
dT = 1/60; % sec, time separation between two frames

outfile = sprintf('april26_pairs_%dto%dmm_3taueta.txt',fromR,toR);

[traj,attr] = readTrajHDF_v9(hfile,'x',[],'y',[],'z',[],'u',[],'v',[],'w',[],'t',[],...
    's11',[],'s12',[],'s13',[],'s22',[],'s23',[],'s33',[],'w1',[],'w2',[],'w3',[],...
    'frames',[1 3000]);
disp('');


% for coarse-graining
t = cat(1,traj.t);
x = cat(1,traj.x);
y = cat(1,traj.y);
z = cat(1,traj.z);



s12 = cat(1,traj.s12);
s13 = cat(1,traj.s13);
s23 = cat(1,traj.s23);

w1 = cat(1,traj.w1);
w2 = cat(1,traj.w2);
w3 = cat(1,traj.w3);


ux = cat(1,traj.s11);
uy = s12 - .5*w3;
uz = s13 + .5*w2;
vx = s12 + .5*w3;
vy = cat(1,traj.s22);
vz = s23 - .5*w1;
wx = s13 - .5*w2;
wy = s23 + .5*w1;
wz = cat(1,traj.s33);

clear s12 s13 s23 w1 w2 w3




numTraj = length(traj);
maxLen    = max(cat(1,traj.trajLen));
% ''''''''''''''''''''''''''''
% outfile = 'mai10_pairs_5to6mm_3taueta.txt';
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

                buffer(:,8)  = traj(i).s11(iA(k)); % derivatives, dudx
                buffer(:,9)  = traj(i).s12(iA(k)) - 0.5*traj(i).w3(iA(k));     % dudy
                buffer(:,10) = traj(i).s13(iA(k)) + 0.5*traj(i).w2(iA(k));    % dudz
                buffer(:,11) = traj(i).s12(iA(k)) + 0.5*traj(i).w3(iA(k));
                buffer(:,12) = traj(i).s22(iA(k));
                buffer(:,13) = traj(i).s23(iA(k)) - 0.5*traj(i).w1(iA(k));
                buffer(:,14) = traj(i).s13(iA(k)) - 0.5*traj(i).w2(iA(k));
                buffer(:,15) = traj(i).s23(iA(k)) + 0.5*traj(i).w1(iA(k));
                buffer(:,16) = traj(i).s33(iA(k));

                buffer(:,17)  = traj(j).s11(iB(k));                           % dudx
                buffer(:,18)  = traj(j).s12(iB(k)) - 0.5*traj(j).w3(iB(k));   % dudy
                buffer(:,19) = traj(j).s13(iB(k)) + 0.5*traj(j).w2(iB(k));    % dudz
                buffer(:,20) = traj(j).s12(iB(k)) + 0.5*traj(j).w3(iB(k));    % dvdx
                buffer(:,21) = traj(j).s22(iB(k));                            % dvdy
                buffer(:,22) = traj(j).s23(iB(k)) - 0.5*traj(j).w1(iB(k));    % dvdz
                buffer(:,23) = traj(j).s13(iB(k)) - 0.5*traj(j).w2(iB(k));    % dwdx
                buffer(:,24) = traj(j).s23(iB(k)) + 0.5*traj(j).w1(iB(k));    % dwxy
                buffer(:,25) = traj(j).s33(iB(k));                            % dwdz


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
    fclose(fid)
end

fclose(fid)


%{
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


ind = find(rand(numPairs,1) > .99)
figure, hold on
for k = 1:length(ind), i = ind(k);
quiver3(pair(i).x,pair(i).y,pair(i).z,pair(i).l1,pair(i).l2,pair(i).l3,'b')
end

%}







