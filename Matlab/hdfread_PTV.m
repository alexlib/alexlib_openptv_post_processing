% Global
            f = @(x1,y1,z1,x2,y2,z2)([x2(:)-x1(:),y2(:)-y1(:),z2(:)-z1(:)]);
            ar = @(r)(line([0 r(1)],[0 r(2)],[0 r(3)]));

% Data location
if findstr(computer,'PCWIN')
    cd('C:\Documents and Settings\user\My Documents\MATLAB\hdf');
else
    cd('/host/Documents and Settings/user/My Documents/MATLAB/HDF/')
end
% read the file info, later we'll add checks of contents
finfo = hdfinfo('0104_hompol.hdf')';
numFrames = length(finfo.Vdata);

[x,y,z,u,v,w,t] = deal([]);

for i = 1:1 % 5
    traj = sprintf('/%d',i);
    A1 = hdfread('0104_water.hdf',traj, 'Fields', 'x,y,z,u,v,w,t', 'FirstRecord',1);
    x = cat(2,x,A1{1});
    y = cat(2,y,A1{2});
    z = cat(2,z,A1{3});
    u = cat(2,u,A1{4});
    v = cat(2,v,A1{5});
    w = cat(2,w,A1{6});
    t = cat(2,t,A1{7});
    
    maxN = 1e6;
    [uplus,uminus] = deal(zeros(maxN,3));
    f1 = uplus(:,1);
    lr = f1;
    k = 0;
    
    
    uniT = unique(t);
    for tt = 1:12 % length(uniT)
        p = find(t == uniT(tt));
        for p1 = 1:length(p)-1
            for p2 = p1+1:length(p)
                r = f(x(p1),y(p1),z(p1),x(p2),y(p2),z(p2));
                k = k + 1

                lr(k) = sqrt(sum(r.^2,2));
                
                u1 = norm([u(p1),v(p1),w(p1)])*cosine([u(p1),v(p1),w(p1)],r);
                u2 = norm([u(p2),v(p2),w(p2)])*cosine([u(p2),v(p2),w(p2)],r);
                uplus(k) = (u1 + u2)/2;
                uminus(k,:) = (u2 - u1)/2;
                % $$\langle u_{+}^2 u_{-} \rangle = \langle \epsilon \rangle r/30$$
                f1(k) = uplus(k).^2*uminus(k)/lr(k)*30;
%                f2 = ??? next structure function
            end
        end
    end
end

conditional_average_f1
