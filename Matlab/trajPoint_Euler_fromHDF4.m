function res=trajPoint_Euler_fromHDF4(first,last,tail)

if nargin == 0
    first = 1;
    last = 2997;
    tail = 0;
end

    hfile = 'E:/Matlab/Alex/HDF/0104_hompol.hdf'
    finfo = hdfinfo(hfile);
    numFrames = length(finfo.Vdata);

value=zeros(300,30,2000);
pointInd=zeros(300,1);

for i=first:last
    clear f;
    i
    traj = sprintf('/%d',i);
    A1 = hdfread(hfile,traj, 'Fields', 'x,y,z,u,v,w,ax,ay,az,w1,w2,w3,s11,s12,s13,s22,s23,s33,ww1,ww2,ww3,wws,sss,R,Q,diss,div,trajnum,t,alx,aly,alz', 'FirstRecord',1);
%     name=['D:/data/Riso_micro_PTV/trajPoint.',num2Str(i)];
%     f=load(name);
%     s=size(f);
    s = length(A1{1});
    if s(1,1)>0
        u = A1{4}; % f(:,4);
        v = A1{5}; % f(:,5);
        w = A1{6}; % f(:,6);
        absu=(u.^2+v.^2+w.^2).^0.5;
        ax = A1{7}; % f(:,7);%%%%%%%%%%%
        ay = A1{8}; % f(:,8);%%%%%%%%%%%
        az = A1{9}; % f(:,9);%%%%%%%%%%%
        absa=(ax.^2+ay.^2+az.^2).^0.5;%%%%%%%%%%%%%
        curv=(((v.*az-w.*ay).^2+(w.*ax-u.*az).^2+(u.*ay-v.*ax).^2).^0.5)./(absu.^3);
        
        
        
        reldiv = A1{27}; % f(:,31);
        age = A1{29};    % f(:,32);
        age = age - age(1);
        

        ende=0;
        for iii=1:s(1,1)-1
            start=ende+1;
            ende=start;
            leng=1;
            searching=1;
            while ende<s(1,1)-1 & searching==1
                if age(ende+1)==age(ende)+1
                    ende=ende+1;
                    leng=leng+1;
                else
                    searching=0;
                end
            end
            if ende<s(1,1)
                for ii=start:ende
                    if age(ii)>tail & age(ii)<leng-tail+1 & reldiv(ii)<0.2 & absa(ii)^2/8e-4<100 & curv(ii)<3000
                        pointInd(age(ii))=pointInd(age(ii))+1;
                        for attr=1:30
                            tmp = A1{attr};
                            value(age(ii),attr,pointInd(age(ii)))=tmp(ii); % f(ii,attr);
                        end
                        
                    end
                end
            end
        end    
    end
    %write out
    name=['Euler.',num2Str(i)];
    fid=fopen(name,'w');
    pointInd(1)
    for ii=1:pointInd(1)
       for attr=1:29
           dummy=value(1,attr,ii);
           %[attr dummy]
           fprintf(fid, '%s\t', num2Str(value(1,attr,ii)));
       end
       fprintf(fid, '%s\n', num2Str(value(1,30,ii)));
    end
    fclose(fid);
    %shift upwards
    value(1:199,:,:)=value(2:200,:,:);
    pointInd(1:199,1)=pointInd(2:200,1);
    pointInd(200)=0;
end


