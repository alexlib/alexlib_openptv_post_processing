first=51;
last=500;
n=50;

count=0;
for i=first-n:first-1
    count=count+1;
    if i<10
        name_load=['F:\mean_flow_analysis_in_turbulent_box\WF14\image/C1S000100000',num2Str(i),'.tif'];
    elseif i<100
        name_load=['F:\mean_flow_analysis_in_turbulent_box\WF14\image/C1S00010000',num2Str(i),'.tif'];
    elseif i<1000
        name_load=['F:\mean_flow_analysis_in_turbulent_box\WF14\image/C1S0001000',num2Str(i),'.tif'];
    else
        name_load=['F:\mean_flow_analysis_in_turbulent_box\WF14\image/C1S000100',num2Str(i),'.tif'];
    end
    A(count,1:1024,1:1024)=imread(name_load);
end

for i=first:last
    if mod(i,10)==0
        i
    end
    A(1,:,:)=[];
    if i<10
        name_load=['F:\mean_flow_analysis_in_turbulent_box\WF14\image/C1S000100000',num2Str(i),'.tif'];
    elseif i<100
        name_load=['F:\mean_flow_analysis_in_turbulent_box\WF14\image/C1S00010000',num2Str(i),'.tif'];
    elseif i<1000
        name_load=['F:\mean_flow_analysis_in_turbulent_box\WF14\image/C1S0001000',num2Str(i),'.tif'];
    else
        name_load=['F:\mean_flow_analysis_in_turbulent_box\WF14\image/C1S000100',num2Str(i),'.tif'];
    end
    A(n,1:1024,1:1024)=imread(name_load);
    B=squeeze(max(A,[],1));
    name_save=['F:\mean_flow_analysis_in_turbulent_box\WF14\image/streak_',num2Str(100000+i),'.tif'];
    
    imwrite(B,name_save,'tif');
end

