function checks_beat_hdf_v2(first,last,tail,stepSize) %#ok<*FNDEF>
%checks(10000,18600,0,100)

nPoints = 10000;
nBins = 101;


first = 1;
last = 10; % 2990;
stepSize = 1;

rmsax=0.01;%tail 7, step 1
rmsay=0.01;
rmsaz=0.01;


% Create bins
lbx = linspace(-2*rmsax,2*rmsax,nBins);
lby = linspace(-2*rmsax,2*rmsax,nBins);
lbz = linspace(-2*rmsax,2*rmsax,nBins);


tail = 1;

hfile = ['../HDF/0104_hompol.hdf'];


%Hi Alex, I outcomment everything that is not needed for the figures
% the code is a bit stupid as you have to iterate for rmsax, which is used
% to normalize the 'forPDFax' stuff

% jcurvGradcurvAcc=zeros(201,201);
% jacaLaBx=zeros(101,101);
% jacaLaBy=zeros(101,101);
% jacaLaBz=zeros(101,101);
jsxLi=zeros(nBins,nBins);
jsyLi=zeros(nBins,nBins);
jszLi=zeros(nBins,nBins);
forRMSax=zeros(nPoints,1);
forPDFax=zeros(nPoints,1);
forERRax=zeros(nPoints,1);
forRMSay=zeros(nPoints,1);
forPDFay=zeros(nPoints,1);
forERRay=zeros(nPoints,1);
forRMSaz=zeros(nPoints,1);
forPDFaz=zeros(nPoints,1);
forERRaz=zeros(nPoints,1);
% forMeanDiss=zeros(nPoints,1);%%%%1.9546e-005
%
% statistic=zeros(nPoints,1);

good=1;
allpoints=0;

% a0_x =1.0710
% a0_y =0.6892
% a0_z =1.3392


for i=first:stepSize:last
    clear f;
    i
    good/allpoints
    %     name=['D:/data/Riso_micro_PTV/g_trajPoint.',num2Str(i)];
    %     f=load(name);
    %     s=size(f);
    
    f1 = hdfread(hfile,int2str(i),...
        'Fields','x,y,z,u,v,w,ax,ay,az,w1,w2,w3,s11,s12,s13,s22,s23,s33,ww1,ww2,ww3,wws,sss,R,Q,diss,div,trajnum,t,alx,aly,alz',...
        'FirstRecord',1);
    s = length(f1{1});
    
    % if s(1,1) > 0
    
    % assigning columns of the file to the variables:
    [...
        xr,yr,zr,...                % 1-3
        u,v,w,...                   % 4-6
        ax,ay,az,...                % 7-9
        w1,w2,w3,...                % 10-12
        s11,s12,s13,s22,s23,s33,... % 13-18
        ww1,ww2,ww3,...             % 19-21
        wws,sss,R,Q,diss,...        % 22,23,24,25,26
        reldiv,trajnum,...             % 27,28
        age,...                       % 29 (!)
        alx,aly,alz] = ...          % 30-32
        deal(f1{1:32});
    
    
    
    
    %             xr=f(:,1);
    %             yr=f(:,2);
    %             zr=f(:,3);
    %         u=f(:,4);
    %         v=f(:,5);
    %         w=f(:,6);
    %         absu=(u.^2+v.^2+w.^2).^0.5;
    %         ax=f(:,7);%%%%%%%%%%%
    %         ay=f(:,8);%%%%%%%%%%%
    %         az=f(:,9);
    %         absa=(ax.^2+ay.^2+az.^2).^0.5;
    %         alx=f(:,19);%%%%%%%%%%%
    %         aly=f(:,20);%%%%%%%%%%%
    %         alz=f(:,21);
    %         absal=(alx.^2+aly.^2+alz.^2).^0.5;
    
    %         w1=f(:,10);
    %         w2=f(:,11);
    %         w3=f(:,12);
    %         s11=f(:,13);
    %         s12=f(:,14);
    %         s13=f(:,15);
    %         s22=f(:,16);
    %         s23=f(:,17);
    %         s33=f(:,18);
    ux=s11;
    uy=s12-0.5*w3;
    uz=s13+0.5*w2;
    vx=s12+0.5*w3;
    vy=s22;
    vz=s23-0.5*w1;
    wx=s13-0.5*w2;
    wy=s23+0.5*w1;
    wz=s33;
    %         enstrophy=w1.^2+w2.^2+w3.^2;
    %         absw=enstrophy.^0.5;
    
    
    wS1=w1.*s11+w2.*s12+w3.*s13;
    wS2=w1.*s12+w2.*s22+w3.*s23;
    wS3=w1.*s13+w2.*s23+w3.*s33;
    %         absws=(wS1.^2+wS2.^2+wS3.^2).^0.5;
    %         wws=w1.*wS1+w2.*wS2+w3.*wS3;
    %             leng=1;
    
    acx=u.*ux+v.*uy+w.*uz;%%%%%%%%%%%%%%%%%%%%
    acy=u.*vx+v.*vy+w.*vz;%%%%%%%%%%%%%%%%%%%%
    acz=u.*wx+v.*wy+w.*wz;%%%%%%%%%%%%%%%%%%%%
    %         absac=(acx.^2+acy.^2+acz.^2).^0.5;%%%%%%%%%%%%%%%%%%%%%%%%%
    %             aLamx=w2.*w-w3.*v;%%%%%%%%%%%%%%%%%%%%%%%
    %             aLamy=w3.*u-w1.*w;%%%%%%%%%%%%%%%%%%%%%%%
    %             aLamz=w1.*v-w2.*u;%%%%%%%%%%%%%%%%%%%%%%%%%
    %             absaLam=(aLamx.^2+aLamy.^2+aLamz.^2).^0.5;%%%%%%%%%%%%%%%%%
    %             unx=u./absu;
    %             uny=v./absu;
    %             unz=w./absu;
    %             aun=ax.*unx+ay.*uny+az.*unz;
    %             aPx=aun.*unx;%%%%%%%%%%%%%%%%%%%%%%%%
    %             aPy=aun.*uny;%%%%%%%%%%%%%%%%%%%%%%%%%
    %             aPz=aun.*unz;%%%%%%%%%%%%%%%%%%%%%%%%%%
    %             absaP=(aPx.^2+aPy.^2+aPz.^2).^0.5;%%%%%%%%%%%%%%%%%%%%%
    %             aOx=ax-aPx;%%%%%%%%%%%%%%%%%%%%
    %             aOy=ay-aPy;%%%%%%%%%%%%%%%%%%%
    %             aOz=az-aPz;%%%%%%%%%%%%%%%%%%%%
    %             absaO=(aOx.^2+aOy.^2+aOz.^2).^0.5;%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %             curvAcc=(((v.*az-w.*ay).^2+(w.*ax-u.*az).^2+(u.*ay-v.*ax).^2).^0.5)./(absu.^3);
    %
    %             reldiv=abs(ux+vy+wz)./(abs(ux)+abs(vy)+abs(wz));
    %         age=f(:,32);
    
    %             axx=f(:,22);
    %             axy=f(:,23);
    %             axz=f(:,24);
    %             ayx=f(:,25);
    %             ayy=f(:,26);
    %             ayz=f(:,27);
    %             azx=f(:,28);
    %             azy=f(:,29);
    %             azz=f(:,30);
    %
    %             aw1=azy-ayz;
    %             aw2=axz-azx;
    %             aw3=ayx-axy;
    %             aenstrophy=aw1.^2+aw2.^2+aw3.^2;
    %             absaw=aenstrophy.^0.5;
    %             as11=0.5*(axx+axx);
    %             as12=0.5*(axy+ayx);
    %             as13=0.5*(axz+azx);
    %             as22=0.5*(ayy+ayy);
    %             as23=0.5*(ayz+azy);
    %             as33=0.5*(azz+azz);
    %
    %             diva=as11+as22+as33;
    
    %         reldiv=f(:,31);
    
    
    %       so we have:
    %       absu,     1
    %       absa,     2
    %       absal,    3
    %       absac,    4
    %       absaLam,  5
    %       absaB,    6
    %       absaP,    7
    %       absaO,    8
    
    starts = find(age == i);
    ends = [starts(2:end)-1, s];
    nTraj = length(starts);
    
    for trajInd = 1:nTraj
        start = starts(trajInd);
        ende = ends(trajInd);
        %         leng = ende - start + 1;
        ii = start+tail:ende-tail+1;
        ii = reldiv(ii) < 0.2;
        ngood = sum(ii);
        
        ind = good:good+ngood-1;
        good = good + ngood;
        
        forRMSax(ind,1) = ax(ii).^2.';
        forRMSay(ind,1) = ay(ii).^2.';
        forRMSaz(ind,1) = az(ii).^2.';
        
        
        
        forPDFax(ind,1)=ax(ii)./rmsax;
        forPDFay(ind,1)=ay(ii)./rmsay;
        forPDFaz(ind,1)=az(ii)./rmsaz;
        forERRax(ind,1)=abs(ax(ii)-alx(ii)-acx(ii))./(abs(ax(ii))+abs(alx(ii))+abs(acx(ii)));
        forERRay(ind,1)=abs(ay(ii)-aly(ii)-acy(ii))./(abs(ay(ii))+abs(aly(ii))+abs(acy(ii)));
        forERRaz(ind,1)=abs(az(ii)-alz(ii)-acz(ii))./(abs(az(ii))+abs(alz(ii))+abs(acz(ii)));
        %                                     statistic(good,1)=(alx(ii)*acx(ii)+aly(ii)*acy(ii)+alz(ii)*acz(ii))/(absal(ii)*absac(ii));
        
        indJ = bindex(ax(ii),lbx,1);
        indI = bindex(alx(ii)+acx(ii),lbx,1);
        for j = 1:numel(indJ)
            jsxLi(indI(j),indJ(j))=jsxLi(indI(j),indJ(j))+1;
        end
        
        indJ = bindex(ay(ii),lby,1);
        indI = bindex(aly(ii)+acy(ii),lby,1);
        for j = 1:numel(indJ)
            jsyLi(indI(j),indJ(j))=jsxLi(indI(j),indJ(j))+1;
        end
        indJ = bindex(az(ii),lbz,1);
        indI = bindex(alz(ii)+acz(ii),lbz,1);
        for j = 1:numel(indJ)
            jszLi(indI(j),indJ(j))=jsxLi(indI(j),indJ(j))+1;
        end
        
        %                         indJ=floor((ax(ii)/rmsax+2)*25)+1;
        %                         indI=floor(((alx(ii)+acx(ii))/rmsax+2)*25)+1;
        %
        %
        %                         if indI>0 && indJ>0 && indI<102 && indJ<102
        %                             jsxLi(indI,indJ)=jsxLi(indI,indJ)+1;
        %                         end
        
        %                         indJ=floor((ay(ii)/rmsay+2)*25)+1;
        %                         indI=floor(((aly(ii)+acy(ii))/rmsay+2)*25)+1;
        %                         if indI>0 && indJ>0 && indI<102 && indJ<102
        %                             jsyLi(indI,indJ)=jsyLi(indI,indJ)+1;
        %                         end
        %                         indJ=floor((az(ii)/rmsaz+2)*25)+1;
        %                         indI=floor(((alz(ii)+acz(ii))/rmsaz+2)*25)+1;
        %                         if indI>0 && indJ>0 && indI<102 && indJ<102
        %                             jszLi(indI,indJ)=jszLi(indI,indJ)+1;
        %                         end
        %
    end
end



success=good/allpoints

if good>nPoints
    good=nPoints;
end
rmsax=mean(forRMSax(1:good,1))^0.5
rmsay=mean(forRMSay(1:good,1))^0.5
rmsaz=mean(forRMSaz(1:good,1))^0.5

% diss=2*1e-6*mean(forMeanDiss(1:good,1))
% a0_x=rmsax^2/(3*(diss)^1.5*1e-6^(-0.5))
% a0_y=rmsay^2/(3*(diss)^1.5*1e-6^(-0.5))
% a0_z=rmsaz^2/(3*(diss)^1.5*1e-6^(-0.5))

figure;
hold on;
grid on;
box on;
[N,X] = nhist(forERRax(1:good,1),100);
plot(X,N,'r');
[N,X] = nhist(forERRay(1:good,1),100);
plot(X,N,'g');
[N,X] = nhist(forERRaz(1:good,1),100);
plot(X,N,'b');
xlabel('rel. error a_i')
ylabel('pdf')

err_ax=mean(forERRax(1:good,1))
err_ay=mean(forERRaz(1:good,1))
err_az=mean(forERRay(1:good,1))

figure;
subplot(1,3,1)
hold on;
grid on;
box on;
[N,X] = nhist(forPDFax(1:good,1),1000);
plot(X,N,'r');
[N,X] = nhist(forPDFay(1:good,1),1000);
plot(X,N,'g');
[N,X] = nhist(forPDFaz(1:good,1),1000);
plot(X,N,'b');
xlabel('a_i/rms(a_i)')
ylabel('pdf')
set(gca,'YScale','log')

subplot(1,3,2)
hold on;
grid on;
box on;
[N,X] = nhist(forPDFax(1:good,1),1000);
plot(X,N.*X.^2,'r');
[N,X] = nhist(forPDFay(1:good,1),1000);
plot(X,N.*X.^2,'g');
[N,X] = nhist(forPDFaz(1:good,1),1000);
plot(X,N.*X.^2,'b');
xlabel('a_i/rms(a_i)')
ylabel('a^2pdf')
set(gca,'YScale','log')

subplot(1,3,3)
hold on;
grid on;
box on;
[N,X] = nhist(forPDFax(1:good,1),1000);
plot(X,N.*X.^4,'r');
[N,X] = nhist(forPDFay(1:good,1),1000);
plot(X,N.*X.^4,'g');
[N,X] = nhist(forPDFaz(1:good,1),1000);
plot(X,N.*X.^4,'b');
xlabel('a_i/rms(a_i)')
ylabel('a^4pdf')
set(gca,'YScale','log')

figure
x=-2:0.04:2;
contourf(x,x,jsxLi,10);
xlabel('a_{x}, rel. div<0.2','FontSize',12,'FontName','Times New Roman');
ylabel('a_{l,x}+a_{c,x}','FontSize',12,'FontName','Times New Roman');
axis equal
axis([-2 2 -2 2])
colorbar;
figure
contourf(x,x,jsyLi,10);
xlabel('a_{y}, rel. div<0.2','FontSize',12,'FontName','Times New Roman');
ylabel('a_{l,y}+a_{c,y}','FontSize',12,'FontName','Times New Roman');
axis equal
axis([-2 2 -2 2])
colorbar;
figure
contourf(x,x,jszLi,10);
xlabel('a_{z}, rel. div<0.2','FontSize',12,'FontName','Times New Roman');
ylabel('a_{l,z}+a_{c,z}','FontSize',12,'FontName','Times New Roman');
axis equal
axis([-2 2 -2 2])
colorbar;


save 0104_hompol_check

