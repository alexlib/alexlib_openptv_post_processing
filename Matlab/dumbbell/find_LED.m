function find_LED(imname,dumbbels,verbose)
% use template matching and masking

I = imread(imname);

stats = struct('Centroid',[0,0],'Area',0,'sumg',0,...
    'MajorAxisLength',0,'MinorAxisLength',0);

if verbose
figure, hold on
imshow(I); axis tight 
end

for nLED = 1:length(dumbbels)
    T = dumbbels{nLED};
    [sx,sy] = size(T);
    sx = ceil(sx/2);
    sy = ceil(sy/2);
    
    I_SSD = template_matching(T,I);
    
    [y,x] = ind2sub(size(I_SSD),find(I_SSD==max(I_SSD(:))));
    
    if verbose, plot(x,y,'ro'); end
    
%     
    I(y-sy:y+sy,x-sx:x+sx) = uint8(0);
    
    stats(nLED).Centroid = [x,y];
    stats(nLED).Area = sx*sy*4;
    stats(nLED).sumg = sum(T(:));
    stats(nLED).MajorAxisLength = max([sx,sy])*2;
    stats(nLED).MinorAxisLength = min([sx,sy])*2;
end


fname = [imname,'_targets'];
write_targets(fname,stats);

