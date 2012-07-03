function show_detected_particles(directory,ncam)
% directory
% ncam = scalar for a single camera or vector (up to 4) for all together

if nargin < 1
    directory = '/Users/alex/Documents/PTV/ptv_alex/img_35_try';
    % directory = '/Users/alex/Documents/PTV/ptv_alex/img34_10';
end

for i = ncam
    filenames{i} = getfiles(directory,i);
end

figure,
if ncam > 1
for i = ncam
    ax{i} = subplot(2,2,i);
end
else
    ax{1} = gca;
end

for i = 1:length(filenames{1})
    for j = ncam
        axes(ax{j})
        data = textread([filenames{j}{i},'_targets']);
        imshow(filenames{j}{i},'Parent',ax{j}); hold on
        if length(data) > 1
            scatter(data(2:end,2),data(2:end,3),'r+');
        end
        drawnow
        hold off
    end
end


function filenames = getfiles(d,ncam)
% GETFILES constructs a cell array of file names

% read the files
imlist = dir(sprintf('%s/cam%d.*',d,ncam));
filenames = cell(1,length(imlist));
k = 1;
for i = 1:length(imlist)
    tmp = fullfile(d,imlist(i).name);
    if isempty(findstr(tmp,'_targets')) && isempty(findstr(tmp,'.mat'))
        filenames{k} = tmp;
        k = k + 1;
    end
end
filenames = filenames(1:k-1);


