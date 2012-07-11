%function [] = subtract_divide_rename_modified(varargin)
% Reads all the image between the first and last, provided
% by the user, subtracts from each image the mean image
% Inputs:         none
% Outputs:        none
%
% Usage:
%        >> subtract_mean_image
% point to the desired first and last files
%
% See also: HELP UIGETFILE, IMREAD, IMWRITE

% Author: Alex Liberzon
% Copyright (c) 2004, Alex Liberzon, IHW, ETHZ
% Last modified at: June 17, 2004
%  Version 1.01, on VideoPC
% - takes average of 100 files, than convert and divide all



%%-----Revised by Debashish Saha at: Octerber 9,2010---------%





% [filename1,pathname] = uigetfile('*.tif','First file');
%wd = cd;
%
% cd(pathname);
% [filename2,pathname] = uigetfile('*.tif','Last file');






for i=750:750:750
    for j=1:5


        %     direc=['D:\Deb/Exp_',num2str(k),];
        %     d = ['D:\Deb/Exp_',num2str(k),'/*.tif'];


        direc=['E:\mean_flow_analysis_in_turbulent_box_v2\rpm',num2str(i),'\WF',num2str(j),'\image'];
        d = ['E:\mean_flow_analysis_in_turbulent_box_v2\rpm',num2str(i),'\WF',num2str(j),'\image','\*.tif'];

        dd=dir(d);
        filename1=dd(1).name;
        filename2=dd(end).name;
        [pathstr,name1,ext1,versn] = fileparts(filename1);
        [pathstr,name2,ext2,versn] = fileparts(filename2);
        cd(direc);


        first = eval(name1(end-5:end));
        last = eval(name2(end-5:end));
        num = last - first + 1;
        m=0;
        for k = first:last
            m=m+1;

            tmp = imread([name1(1:end-6),sprintf('%06d',k),'.tif']);
            imwrite(tmp(1:512,1:512),['E:\mean_flow_analysis_in_turbulent_box_v2\rpm',num2str(i),'\WF',num2str(j),'\image','/Cam1.',sprintf('%d',m),'.tif'],'tiff','compression','none');
            imwrite(tmp(1:512,513:1024),['E:\mean_flow_analysis_in_turbulent_box_v2\rpm',num2str(i),'\WF',num2str(j),'\image','/Cam2.',sprintf('%d',m),'.tif'],'tiff','compression','none');
            imwrite(tmp(513:1024,513:1024),['E:\mean_flow_analysis_in_turbulent_box_v2\rpm',num2str(i),'\WF',num2str(j),'\image','/Cam3.',sprintf('%d',m),'.tif'],'tiff','compression','none');
            imwrite(tmp(513:1024,1:512),['E:\mean_flow_analysis_in_turbulent_box_v2\rpm',num2str(i),'\WF',num2str(j),'\image','/Cam4.',sprintf('%d',m),'.tif'],'tiff','compression','none');

        end

        for k=first:last
            tmpname=[direc,'\',name1(1:end-6),sprintf('%06d',k),'.tif'];
            recycle('off');
            delete(tmpname);
        end
[i j]
    end


end







