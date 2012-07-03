if ismac
    matdirectory = 'C:\Users\hadar\Dropbox\resuspension\2011/trajectories(186-194)';
else
    matdirectory ='C:\Documents and Settings\user.EFDL5-PC\My Documents\Dropbox\resuspension\2011\trajectories(186-194)';
end
large_ones = dir(fullfile(matdirectory,'large*'));

for i = 1:length(large_ones)
    plot_type2_large_small_single_figure(large_ones(i).name,'yf','yf',10);
end
