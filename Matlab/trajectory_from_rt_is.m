clear all
close all
%close all;
tau_eta=0.14;
frame_rate=250;
factor=tau_eta*frame_rate;
direc=['G:\flow_induced_aggregates_breakage_at_100rpm\WF45\res'];
d = [direc,'/rt_is*'];

% dd=dir(d);
% filename1=dd(1).name;
% 
% first = eval(filename1(7:end));
% 
% if first==1
%     last = size(dd,1);
% else
%     filename2=dd(end).name;
%     last=eval(filename2(7:end));
% end

  first=1;
 last=450;

x_agg=[];
y_agg=[];
z_agg=[];
time=[];
total_px=[];
x_px=[];
y_px=[];
sumgrv=[];
% load_data=0;
% if load_data==0;
    for i=first:last
        name_rt=[direc,'/rt_is_v2.',num2str(i)];
        fid_rt = fopen(name_rt,'r');
        num_points = fscanf(fid_rt, '%i', [1 1]);
        tmp = fscanf(fid_rt, '%i %f %f %f %i %i %i %i %f %f %f %f ', [12 num_points]);
        A=tmp';
        x_a=A(:,2)*0.001;
        y_a=A(:,3)*0.001;
        z_a=A(:,4)*0.001;
        total_px=A(:,9);
        x_px=A(:,10);
        y_px=A(:,11);
        sumgrv=A(:,12);
        x_agg=[x_agg;x_a(abs(x_a)>0)];
        y_agg=[y_agg;y_a(abs(y_a)>0)];
        z_agg=[z_agg;z_a(abs(z_a)>0)];
        
        
        t=repmat(i-first+1,[numel(A(:,2)) 1]);
        time=[time;t];
        fclose(fid_rt);
    end
    
    %save flow_induced_aggregates_breakage_at_100rpm_WF6 x_agg y_agg z_agg time
%else
    %load flow_induced_aggregates_breakage_at_100rpm_WF6
%end
scatter3(x_agg,y_agg,z_agg)
%plot3(x_agg,y_agg,z_agg)

% figure(100)
% for j=1:size(time,1)-1
%     if time(j+1)==time(j)
%         if ((x_agg(time(j))-x_agg(time(j)+1))^2+(y_agg(time(j))-y_agg(time(j)+1))^2+(z_agg(time(j))-z_agg(time(j)+1))^2)^0.5>2e-3
%             
%         scatter3(x_agg(time(j)),y_agg(time(j)),z_agg(time(j)),'b')
%         hold on
%         scatter3(x_agg(time(j+1)),y_agg(time(j+1)),z_agg(time(j+1)),'r')
%         hold on
%         plot3([x_agg(time(j)); x_agg(time(j+1))],[y_agg(time(j)); y_agg(time(j+1))],[z_agg(time(j)); z_agg(time(j+1))],'g');
%         
%         end
%     end
% end
% 
% 
% 
% 
% k=[time x_agg y_agg z_agg]
% 















