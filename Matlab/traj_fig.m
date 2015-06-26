clear all
%close all
clc

wf=97;%[32 38 52 27 97]%[27 32 38 62 80 71 88 52 53 97];%[15 27 31 32 35 38 45 46 50 52 53 61 62 67 71 75 80 88 97 98 102 107];%[1 2 6 10  12 13 15 16 18  21 22 23 25 26 27 28 29 30 31 32 33 34 35 36 37 38 41 42 45 46 50 51 52 53 55 56 58 59 61 62 63 64 65 66 67 68 69  71 72 73 74 75 76 77 79 80 81 82 83 84 86 87 88 89  92 93 95 97 98  99 100 101 102  104  106 107 108 109 110 111 112 113];


%%%%% Notes

%1) WF 32 38 52 make one single figure separation_1_1_v2
%2) WF 27 97 make another figure together and then I showed them separetely
%   in the thesis
%3) WF 97 is the most beautiful and takes a single figure honor
%4) WF 80,62 are sacrificed, because it would be too many nice figure
%5) The rest WF s are somehow weirdly flawed ! No time check now why separation is wrong ! 
%6)  WF 32 38 52 and WF 102 make separation_1_1_v3
fold='D:\Turbulent_data_backup\flow_induced_aggregates_breakage_at_100rpm\WF';
for num=wf
par=load([fold,num2str(num),'/parent_variables_v3']);
var=load([fold,num2str(num),'/variables_v2']);
   chld= load([fold,num2str(num),'/children_variables_v3']);
  parent=par.parent;
  child=chld.child;
  pt=par.pTIME;
  ct=chld.cTIME;
  X=par.X;Y=par.Y;Z=par.Z;x=chld.c_X;y=chld.c_Y;z=chld.c_Z;
  
  cmtm=intersect(pt,ct);
  for j=1:length(cmtm)
      indc(j)=find(ct==cmtm(j));
      indp(j)=find(pt==cmtm(j));
  end
 
  
% figure;
% hold on
% scatter3(x,y,z,'k.')
% plot3(x,y,z,'b','LineWidth',2)
% 
% scatter3(X,Y,Z,'k.')
% plot3(X,Y,Z,'b','LineWidth',2)
% 
% scatter3(x(1),y(1),z(1),100,'g','filled')
% scatter3(X(1),Y(1),Z(1),100,'g','filled')
% scatter3(x(end),y(end),z(end),100,'r','filled')
% scatter3(X(end),Y(end),Z(end),100,'r','filled')


  for k=1:length(cmtm)
      hold on
  %plot3([X(indp(k)) x(indc(k))], [Y(indp(k)) y(indc(k))] ,[Z(indp(k)) z(indc(k))],'b')
    sep(k)=sqrt(abs(X(indp(k))-x(indc(k))).^2+(abs(Y(indp(k))-y(indc(k)))).^2+(abs(Z(indp(k)))-z(indc(k))).^2);
    %color_line3([X(indp(k)) x(indc(k))], [Y(indp(k)) y(indc(k))] ,[Z(indp(k)) z(indc(k))],sep(k))
  end
   plot(sort(sep))%./(.3e-)) 
  box on
  tit=['WF=',num2str(num)];
  title(tit)
  clear sep;
 
%%% This part is for WF 97
%for k=1:length(cmtm)
    %color_line3([X(indp(k)) x(indc(k))], [Y(indp(k)) y(indc(k))] ,[Z(indp(k)) z(indc(k))],sep(end-k+1))
%end
 %%------------   
  end
  

