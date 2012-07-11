%function [ children_x,children_y,children_z,res]= xyz(varargin)

% test case: xyz(19740,583,267)
function  xyz(varargin)
%[ the entire results]= distance(n,xpix1,ypix1)

%n,xpix,ypix
% Now this piece of code can take multiple points from different camera and
% pretty flexible whether the user provide single input or multiple input

% function varargout{1:nargout}= distance(varargin)


n=varargin{1};% This is my frame number
if nargin<2
    disp('You did not choose any point on image')
end
if nargin < 4
    xpix1=varargin{2};
    ypix1=varargin{3};
else
    if nargin < 6
        xpix1=varargin{2};
        ypix1=varargin{3};
        xpix2=varargin{4};
        ypix2=varargin{5};
    else
        if nargin<8
            xpix1=varargin{2};
            ypix1=varargin{3};
            xpix2=varargin{4};
            ypix2=varargin{5};
            xpix3=varargin{6};
            ypix3=varargin{7};
        else
            xpix1=varargin{2};
            ypix1=varargin{3};
            xpix2=varargin{4};
            ypix2=varargin{5};
            xpix3=varargin{6};
            ypix3=varargin{7};
            xpix4=varargin{8};
            ypix4=varargin{9};
        end
    end
end



%%%%%%%%%%%% Target files %%%%%%%%%%
% name_t1=['E:\PTV\Working_folder\WF_1\img_for_full_processing\cam1.',num2Str(n)];
% name_t2=['E:\PTV\Working_folder\WF_1\img_for_full_processing\cam2.',num2Str(n)];
% name_t3=['E:\PTV\Working_folder\WF_1\img_for_full_processing\cam3.',num2Str(n)];
% name_t4=['E:\PTV\Working_folder\WF_1\img_for_full_processing\cam4.',num2Str(n)];
% 


name_t1=['E:\PTV\Working_folder\WF_1\img_for_full_processing\cam1.',num2Str(n)];
name_t2=['E:\PTV\Working_folder\WF_1\img_for_full_processing\cam2.',num2Str(n)];
name_t3=['E:\PTV\Working_folder\WF_1\img_for_full_processing\cam3.',num2Str(n)];
name_t4=['E:\PTV\Working_folder\WF_1\img_for_full_processing\cam4.',num2Str(n)];















for k=1:4

    switch k
        case 1
            name=name_t1;
        case 2
            name=name_t2;
        case 3
            name=name_t3;
        case 4
            name=name_t4;
    end

    name=[name '_targets'];

    fid_target = fopen(name);
    number_target = fscanf(fid_target, '%i', [1 1]);
    g=fscanf(fid_target, '%i %f %f %g %g %g %g %g', [8 inf]);
    g=g';
    fclose(fid_target);
    si_target=size(g);


    switch k
        case 1
            if si_target(1,1)>0
                idt1=round(g(:,1));
                xt1=g(:,2);
                yt1=g(:,3);
                n_px1=g(:,4);
                nx_px1=g(:,5);
                ny_px1=g(:,6);
                sumg1=g(:,7);
                ptr_rt1=g(:,8);
            end
        case 2
            if si_target(1,1)>0
                idt2=round(g(:,1));
                xt2=g(:,2);
                yt2=g(:,3);
                n_px2=g(:,4);
                nx_px2=g(:,5);
                ny_px2=g(:,6);
                sumg2=g(:,7);
                ptr_rt2=g(:,8);
            end
        case 3
            if si_target(1,1)>0
                idt3=round(g(:,1));
                xt3=g(:,2);
                yt3=g(:,3);
                n_px3=g(:,4);
                nx_px3=g(:,5);
                ny_px3=g(:,6);
                sumg3=g(:,7);
                ptr_rt3=g(:,8);
            end
        case 4
            if si_target(1,1)>0
                idt4=round(g(:,1));
                xt4=g(:,2);
                yt4=g(:,3);
                n_px4=g(:,4);
                nx_px4=g(:,5);
                ny_px4=g(:,6);
                sumg4=g(:,7);
                ptr_rt4=g(:,8);
            end
    end
end





name_rt=['E:\PTV\Working_folder\WF_1\res\rt_is.',num2Str(n)];
fid_rt = fopen(name_rt);
number_rt = fscanf(fid_rt, '%i', [1 1]);
f=fscanf(fid_rt, '%i %f %f %f %g %g %g %g', [8 inf]);
f=f';
fclose(fid_rt);
si_rt=size(f);

if si_rt(1,1)>0
    id=round(f(:,1));
    x=f(:,2);
    y=f(:,3);
    z=f(:,4);
    c1=round(f(:,5));
    c2=round(f(:,6));
    c3=round(f(:,7));
    c4=round(f(:,8));
end





% d1=sqrt((xt1-xpix1).^2+(yt1-ypix1).^2);
%
%
%
% parent=ptr_rt1(find(d==min(d)));
% if sign(parent)~=1
%     error('Indentifier in target leading to rt_is is negative !!!!')
% else
%     ind_in_rtis=parent+1;
%     children_x=x(ind_in_rtis)/1000;
%     children_y=y(ind_in_rtis)/1000;
%     children_z=z(ind_in_rtis)/1000;
%
% end

%%%%%--------
if nargin<4
    d1=sqrt((xt1-xpix1).^2+(yt1-ypix1).^2);
    parent1=ptr_rt1(find(d1==min(d1)));
    total_pix1=n_px1(find(d1==min(d1)));
    xpix1=nx_px1(find(d1==min(d1)));
    ypix1=ny_px1(find(d1==min(d1)));
    
else
    if nargin<6

        d1=sqrt((xt1-xpix1).^2+(yt1-ypix1).^2);
        parent1=ptr_rt1(find(d1==min(d1)));
        total_pix1=n_px1(find(d1==min(d1)));
        xpix1=nx_px1(find(d1==min(d1)));
        ypix1=ny_px1(find(d1==min(d1)));

        d2=sqrt((xt2-xpix2).^2+(yt2-ypix2).^2);
        parent2=ptr_rt2(find(d2==min(d2)));
        total_pix2=n_px2(find(d2==min(d2)));
        xpix2=nx_px2(find(d2==min(d2)));
        ypix2=ny_px2(find(d2==min(d2)));

    else
        if nargin<8

            d1=sqrt((xt1-xpix1).^2+(yt1-ypix1).^2);
            parent1=ptr_rt1(find(d1==min(d1)));
            total_pix1=n_px1(find(d1==min(d1)));
            xpix1=nx_px1(find(d1==min(d1)));
            ypix1=ny_px1(find(d1==min(d1)));

            d2=sqrt((xt2-xpix2).^2+(yt2-ypix2).^2);
            parent2=ptr_rt2(find(d2==min(d2)));
            total_pix2=n_px2(find(d2==min(d2)));
            xpix2=nx_px2(find(d2==min(d2)));
            ypix2=ny_px2(find(d2==min(d2)));


            d3=sqrt((xt3-xpix3).^2+(yt3-ypix3).^2);
            parent3=ptr_rt3(find(d3==min(d3)));
            total_pix3=n_px3(find(d3==min(d3)));
            xpix3=nx_px3(find(d3==min(d3)));
            ypix3=ny_px3(find(d3==min(d3)));

            
            
        else
            d1=sqrt((xt1-xpix1).^2+(yt1-ypix1).^2);
            parent1=ptr_rt1(find(d1==min(d1)));
            total_pix1=n_px1(find(d1==min(d1)));
            xpix1=nx_px1(find(d1==min(d1)));
            ypix1=ny_px1(find(d1==min(d1)));


            d2=sqrt((xt2-xpix2).^2+(yt2-ypix2).^2);
            parent2=ptr_rt2(find(d2==min(d2)));
            total_pix2=n_px2(find(d2==min(d2)));
            xpix2=nx_px2(find(d2==min(d2)));
            ypix2=ny_px2(find(d2==min(d2)));

            d3=sqrt((xt3-xpix3).^2+(yt3-ypix3).^2);
            parent3=ptr_rt3(find(d3==min(d3)));
            total_pix3=n_px3(find(d3==min(d3)));
            xpix3=nx_px3(find(d3==min(d3)));
            ypix3=ny_px3(find(d3==min(d3)));

            d4=sqrt((xt4-xpix4).^2+(yt4-ypix4).^2);
            parent4=ptr_rt4(find(d4==min(d4)));
            total_pix4=n_px4(find(d4==min(d4)));
            xpix4=nx_px4(find(d4==min(d4)));
            ypix4=ny_px4(find(d4==min(d4)));
            
            
        end
    end
end



if sign(parent1)==1

    ind_in_rtis1=parent1+1;
    children_x=x(ind_in_rtis1)/1000;
    children_y=y(ind_in_rtis1)/1000;
    children_z=z(ind_in_rtis1)/1000;
    disp('1st cam is successful')
    
    res=inter(children_x,children_y,children_z,total_pix1,xpix1,ypix1);
    return;
    %     error('Indentifier in target leading to rt_is is negative !!!!')
else
    if sign(parent2)==1

        ind_in_rtis2=parent2+1;
        children_x=x(ind_in_rtis2)/1000;
        children_y=y(ind_in_rtis2)/1000;
        children_z=z(ind_in_rtis2)/1000;
        disp('1st cam is blind and 2nd cam is successful')
        res=inter(children_x,children_y,children_z,total_pix2,xpix2,ypix2);
        return;

    else

        if sign(parent3)==1

            ind_in_rtis3=parent3+1;
            children_x=x(ind_in_rtis3)/1000;
            children_y=y(ind_in_rtis3)/1000;
            children_z=z(ind_in_rtis3)/1000;
            disp('1st and 2nd cam are blind,3rd cam is successful')
            res=inter(children_x,children_y,children_z,total_pix3,xpix3,ypix3);
            return;
        else
            if sign(parent4)==1

                ind_in_rtis4=parent4+1;
                children_x=x(ind_in_rtis4)/1000;
                children_y=y(ind_in_rtis4)/1000;
                children_z=z(ind_in_rtis4)/1000;
                disp('1st,2nd and 3rd cam are blind,4th cam is successful')
                res=inter(children_x,children_y,children_z,total_pix4,xpix4,ypix4);
            else
                disp('All the selected points from four camera have no rt_is !!')

            end
        end
    end
end

%inter(children_x,children_y,children_z);













