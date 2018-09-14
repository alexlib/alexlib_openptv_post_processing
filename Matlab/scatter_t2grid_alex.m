
%This function is used to create a Eulieran data of the 3D position and 
%velocities %of the particles from the scatter data obtained by traj2scatter t.m 
%Its required inputs are the scatter data and the grid paramters  for  the  positions .  
%The  optional  variable ,  q,  can  be  
%used  to  specify  the  number  of  data  files  used .  
%This function  will  output  
%the  3D position  and  velocity  grid  (xg ,  yg ,  zg ,  ug ,  vg ,  
%wg)  and  the  average  
%velocities  (u avg ,  v avg ,  w avg) .  

%Input:[scatter t,x i,x interval,x f,y i,y interval,y f ,z i,z interval ,z f,number of files(optional)] 
%Output:[xg,yg,zg,ug,vg,wg,u avg,v avg,w avg] 
%Jin Tae Kim 2014 University of Illinois at Urbana Champaign 

function [varargout] = scatter_t2grid( scatter_t, x1 , x2 ,x3,y1, y2,y3, z1, z2, z3, q) 
 
if nargin == 10 
    % this is the function with the absence of the optional variable 
    
    [xg,yg,zg] = ndgrid(x1:x2:x3,y1:y2:y3,z1:z2:z3); %creates a grid with the specified variables 
    
     for i=1:length(scatter_t) % uses every data file in the specified directory and tracks how fast the program is running 
            x= scatter_t(i).x; %%to avoid the interpolation around zero 
            y= scatter_t(i).y; 
            z= scatter_t(i).z; 
            u= scatter_t(i).u; 
            v= scatter_t(i).v; 
            w= scatter_t(i).w; 
            
            ind = (find(abs(u) > 10e5& abs(v) > 10e5& abs(w) > 10e5)) ; % for outliers 
            if length(ind)>100 
                ug(:,:,:,i) = griddata(x,y,z,u,xg,yg,zg,' natural ') ; % interpolates scattered data and creates grid for the velocities 
                vg(:,:,:,i) = griddata(x,y,z,v,xg,yg,zg,'natural ') ; 
                wg(:,:,:,i) = griddata(x,y,z,w,xg,yg,zg,'natural ') ; 
            end
     end
end
i
%ind %the above two lines can be used or commented out to track how fast the 
%program is running ( i ) and what the outliers are (ind) 


if (nargin > 10) % this is the function with the presence of the optional variable
    
    [xg,yg,zg] = ndgrid(x1:x2:x3,y1:y2:y3,z1:z2:z3); %creates a grid with the specified variables 
    for i=1:q % uses the specified number of data files in the specified directory and tracks how fast the program is running 
        x= scatter_t(i).x; %to avoid the interpolation around zero 
        y= scatter_t(i).y; 
        z= scatter_t(i).z; 
        u= scatter_t(i).u; 
        v= scatter_t(i).v; 
        w= scatter_t(i).w; 
        ind = (find(abs(u) > 10e5 & abs(v) > 10e5 & abs(w) > 10e5)) ; % for outliers 
        if length(ind)>100 
            ug(:,:,:,i) = griddata(x,y,z,u,xg,yg,zg,'natural ') ; % interpolates scattered data and creates grid for the velocities 
            vg(:,:,:,i) = griddata(x,y,z,v,xg,yg,zg,' natural ') ; 
            wg(:,:,:,i) = griddata(x,y,z,w,xg,yg,zg,' natural ') ; 
        end
    end
end
i
%ind 
%the above two lines can be used or commented out to track how fast the program is running ( i ) 
%and what the outliers are (ind) 

u_avg = nanmean(ug,4); % calculates the mean velocity with respect to time 
v_avg = nanmean(vg ,4); 
w_avg = nanmean(wg,4); 
varargout{1}=xg; % the outputs 
varargout{2}=yg; 
varargout{3}=zg ; 
varargout{4}=ug; 
varargout{5}=vg; 
varargout{6}=wg; 
varargout{7}=u_avg; 
varargout{8}=v_avg; 
varargout{9}=w_avg; 
end 