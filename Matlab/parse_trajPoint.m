function [xr,yr,zr,...
    u,v,w,...
    ax,ay,az,...
    alx,aly,alz,...
    w1,w2,w3,...
    s11,s12,s13,s22,s23,s33,...
    axx,axy,axz,ayx,ayy,ayz,azx,azy,azz,...
    reldiv, age] = parse_trajPoint(f)

% from post_process.cpp, in macosx/ branch
%{
                     fprintf(fpp, "%lf\t", xp[ii]);//1
                     fprintf(fpp, "%lf\t", yp[ii]);//2
                     fprintf(fpp, "%lf\t", zp[ii]);//3
                     fprintf(fpp, "%lf\t", up[ii]);//4
                     fprintf(fpp, "%lf\t", vp[ii]);//5
                     fprintf(fpp, "%lf\t", wp[ii]);//6
                     fprintf(fpp, "%lf\t", axp[ii]);//7
                     fprintf(fpp, "%lf\t", ayp[ii]);//8
                     fprintf(fpp, "%lf\t", azp[ii]);//9
                     fprintf(fpp, "%lf\t", w1p[ii]);//10
                     fprintf(fpp, "%lf\t", w2p[ii]);//11
                     fprintf(fpp, "%lf\t", w3p[ii]);//12
                     fprintf(fpp, "%lf\t", s11p[ii]);//13
                     fprintf(fpp, "%lf\t", s12p[ii]);//14
                     fprintf(fpp, "%lf\t", s13p[ii]);//15
                     fprintf(fpp, "%lf\t", s22p[ii]);//16
                     fprintf(fpp, "%lf\t", s23p[ii]);//17
                     fprintf(fpp, "%lf\t", s33p[ii]);//18
					 fprintf(fpp, "%lf\t", utp[ii]);//19
                     fprintf(fpp, "%lf\t", vtp[ii]);//20
                     fprintf(fpp, "%lf\t", wtp[ii]);//21
                     fprintf(fpp, "%lf\t", daxdxp[ii]);//22
                     fprintf(fpp, "%lf\t", daxdyp[ii]);//23
                     fprintf(fpp, "%lf\t", daxdzp[ii]);//24
					 fprintf(fpp, "%lf\t", daydxp[ii]);//25
                     fprintf(fpp, "%lf\t", daydyp[ii]);//26
                     fprintf(fpp, "%lf\t", daydzp[ii]);//27
					 fprintf(fpp, "%lf\t", dazdxp[ii]);//28
                     fprintf(fpp, "%lf\t", dazdyp[ii]);//29
                     fprintf(fpp, "%lf\t", dazdzp[ii]);//30
                     fprintf(fpp, "%lf\t", quality);//31 0=good, 1=bad
                     fprintf(fpp, "%lf\n", (double)(ii));//32 age along trajectory
%}


xr=f(:,1);
yr=f(:,2);
zr=f(:,3);

u=f(:,4);
v=f(:,5);
w=f(:,6);

ax=f(:,7);
ay=f(:,8);
az=f(:,9);



w1=f(:,10);
w2=f(:,11);
w3=f(:,12);

s11=f(:,13);
s12=f(:,14);
s13=f(:,15);
s22=f(:,16);
s23=f(:,17);
s33=f(:,18);

alx=f(:,19);
aly=f(:,20);
alz=f(:,21);




axx=f(:,22);
axy=f(:,23);
axz=f(:,24);
ayx=f(:,25);
ayy=f(:,26);
ayz=f(:,27);
azx=f(:,28);
azy=f(:,29);
azz=f(:,30);


reldiv = f(:,31);


age=f(:,32);
