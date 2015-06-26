% load 0104_water


nu = 1e-6; % kinematic viscosity of water/polymer solution

diss = mean(s11.^2 + s22.^2 + s33.^2 + ...
    2*s12.^2 + 2*s13.^2 + 2*s23.^2)*2*nu 

tau_eta = (nu/diss)^(1/2)
eta = ((nu)^3/diss)^(1/4)

u_rms = mean(1/3*(u.^2 + v.^2 + w.^2)).^(1/2)



nBins = 1000;

ux=s11;
uy=s12-0.5*w3;
uz=s13+0.5*w2;
vx=s12+0.5*w3;
vy=s22;
vz=s23-0.5*w1;
wx=s13-0.5*w2;
wy=s23+0.5*w1;
wz=s33;



acx=u.*ux+v.*uy+w.*uz;
acy=u.*vx+v.*vy+w.*vz;
acz=u.*wx+v.*wy+w.*wz;


ind = reldiv < 0.1; 
sum(ind)/length(ind)



al = [alx,aly,alz]; al = al(ind,:);
ac = [acx,acy,acz]; ac = ac(ind,:);
figure, nhist(cosine(al,ac),500,'r-')


a = [ax,ay,az]; a = a(ind,:); V = [u,v,w]; V = V(ind,:);
figure, nhist(cosine(a,V),500);








t = cumsum(age);

reldiv=abs(ux+vy+wz)./(abs(ux)+abs(vy)+abs(wz));
figure, plot(t,reldiv);
figure, nhist(reldiv,101);


figure, plot(t,ax,'r',t,alx+acx,'g'); %,t,acx,'b');



% cosalac = cosine([alx,aly,alz],[acx,acy,acz]);
% nhist(cosalac,101)


% cos(\omega, W)

wS1=w1.*s11+w2.*s12+w3.*s13;
wS2=w1.*s12+w2.*s22+w3.*s23;
wS3=w1.*s13+w2.*s23+w3.*s33;
absws=(wS1.^2+wS2.^2+wS3.^2).^0.5;
wws=w1.*wS1+w2.*wS2+w3.*wS3;

figure, nhist(cosine([w1,w2,w3],[wS1,wS2,wS3]),101);


% pdf of the acceleration components


% a_{Lamb}
aLamx=w2.*w-w3.*v;
aLamy=w3.*u-w1.*w;
aLamz=w1.*v-w2.*u;
absaLam=(aLamx.^2+aLamy.^2+aLamz.^2).^0.5;


% a_{B}
aBx=acx-aLamx;
aBy=acy-aLamy;
aBz=acz-aLamz;
absaB=(aBx.^2+aBy.^2+aBz.^2).^0.5;


% unit vector in the direction of velocity
absu=(u.^2+v.^2+w.^2).^0.5;
unx=u./absu;
uny=v./absu;
unz=w./absu;
aun=ax.*unx+ay.*uny+az.*unz;

% a_{\parallel}
aPx=aun.*unx;
aPy=aun.*uny;
aPz=aun.*unz;
absaP=(aPx.^2+aPy.^2+aPz.^2).^0.5;


% a_{\perp}
aOx=ax-aPx;
aOy=ay-aPy;
aOz=az-aPz;
absaO=(aOx.^2+aOy.^2+aOz.^2).^0.5;


            absa=(ax.^2+ay.^2+az.^2).^0.5;
            
 
            absal=(alx.^2+aly.^2+alz.^2).^0.5;
            
            absac = (acx.^2 + acy.^2 + acz.^2).^0.5;

            nBins = 1001;
figure,
nhist(absa,nBins,'r')
hold on
nhist(absal,nBins,'g')
nhist(absac,nBins,'b');
nhist(absaLam,nBins,'k')
nhist(absaB,nBins,'m')
nhist(absaP,nBins,'c');
set(gca,'yscale','log','xscale','log')


% a_{\ell}
alx1 = ax - acx;
aly1 = ay - acy;
alz1 = az - acz;
absal1=(alx1.^2+aly.^2+alz.^2).^0.5;

absal = (alx.^2 + aly.^2 + alz.^2).^0.5;

figure, nhist(absal1,nBins),set(gca,'yscale','log','xscale','log');
hold on
nhist(absal,nBins,'r');


jointPDF(alx,alx1,51,51,'a_{lx}','a_{x} - a_{cx}')





