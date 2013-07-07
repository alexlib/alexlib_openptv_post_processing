function [bin_centers,average_in_bins,n] = plot_average_in_bins(a,edges)
% PLOT_AVERAGE_IN_BINS(DATA,EDGES) plots the average in each bin, according to
% the edges provided or the number of bins (N). The DATA is divided in bins
% according to the first (!) column and the average of the second column
% is presented
% 
% Usage:
% a = rand(100,2)
% a(:,2)  = a(:,2) + 100;
% edges = 0:.1:1;
% [bins,ave] = plot_average_in_bins(a,edges);
% % OR
%  [bins,ave] = plot_average_in_bins(a,20); % 20 bins

% Author: Alex Liberzon
% 01-03-2012
% 

if length(edges) == 1
    edges = linspace(min(a(:,1)),max(a(:,1)),edges); 
end

[n,which_bin] = histc(a(:,1),edges);% n is the number of values in each bin, which_bin - which bin the value belongs to
average_in_bins = zeros(size(n)); % size n -  number of bins

for i = 1:length(n)-1, 
    ind = which_bin == i; 
% average_in_bins(i) = mean(a(ind,2)); 
average_in_bins(i) = mean(a(ind,3)); 
%average_in_bins(i) = mean(a(ind,4)); 
end

bin_centers = edges(1:end-1) + diff(edges)/2;

% figure
% plot(bin_centers,average_in_bins(1:end-1),'o')
