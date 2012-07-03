% statistics_of_type2

matdirectory ='C:\Documents and Settings\user.EFDL5-PC\My Documents\Dropbox\resuspension\2011\trajectories(186-194)';
%'/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
large_ones = dir(fullfile(matdirectory,'large*'));
collected_data_large = [];
collected_data_small = [];

% CHANGE EVERY RUN

time_interval_before_event = 10;
close_enough = 30;
% name_of_the_function = 'type2_quantity' % 'type2_abs_u', 'type_2';
% quantity_large = 'xf'
% quantity_small = 'vf'

% if Reynolds stresses are required, one can use:
name_of_the_function ='type2_before_acceleration_correlation'; % 'type2_ReynoldsStress';
quantity_large = 'ayf';
quantity_small = 'ayf'; % or 'uv', 'vw'


hf2 = figure; hold on;

for i = 1:length(large_ones)
    num = large_ones(i).name(6:end-4);
    load(fullfile(matdirectory,large_ones(i).name));
    try
        load(fullfile(matdirectory,strrep(large_ones(i).name,'large','small')));
        disp(['large',num])
        % CHANGE THE NAME OF THE FILE THAT IS RUNNING
        [data_large,data_small] = feval(name_of_the_function,...
             ...quantity_large,quantity_small,...
            eval(['large',num]),eval(['small',num]),...
            time_interval_before_event,close_enough,0,hf2);
        %             [data_large,data_small] = feval('type2_abs_u',eval(['large',num]),eval(['small',num]),20,20,0);
        
        collected_data_large = cat(1,collected_data_large,data_large);
        collected_data_small = cat(1,collected_data_small,data_small);
    catch
        continue
    end
    %keyboard
    
end


if strcmp(name_of_the_function,'type2_ReynoldsStress')
    % if Reynolds stresses are used, then we must collect the data instead of
    % averaging it and show it as histograms
    
    nhist(collected_data_small,11,'b-o')
else
    
    % demo figure - at the moment, type2 returns the u component of the large
    % and small particles
    figure
    valid = ~isnan(collected_data_small);
    plot(collected_data_large(valid),collected_data_small(valid),'o');
    axis equal
    xlabel([quantity_large,'_p']); ylabel(quantity_small)
end

