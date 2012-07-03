function [stats1,stats2] =  colloid_image_segmentation(f1,rect,verbose)
I = imcrop(f1,rect);
I_eq = adapthisteq(I);
% imshow(I_eq)
bw = im2bw(I_eq, .5*graythresh(I_eq));
% imshow(bw)
bw2 = imfill(bw,'holes');
bw3 = imopen(bw2, ones(3,3));
bw4 = bwareaopen(bw3, 10);
bw4_perim = bwperim(bw4);
% overlay1 = imoverlay(I_eq, bw4_perim, [.3 1 .3]);
% imshow(overlay1),
mask_em = imextendedmax(I_eq, 60);
% imshow(mask_em);
mask_em = imclose(mask_em, ones(3,3));
mask_em = imfill(mask_em, 'holes');
mask_em = bwareaopen(mask_em, 10);
overlay2 = imoverlay(I_eq, bw4_perim | mask_em, [.3 1 .3]);
if verbose
    figure;
    imshow(overlay2)
    axis on
end
% I_eq_c = imcomplement(I_eq);
% I_mod = imimposemin(I_eq_c, ~bw4 | mask_em);
L1 = bwlabel(bw4);
L2 = bwlabel(mask_em);

stats1 = regionprops(L1);
stats2  = regionprops(L2);

if verbose
    hold on
for i = 1:length(stats1)
    plot(stats1(i).Centroid(1),stats1(i).Centroid(2),'r.');
end
for i = 1:length(stats2)
    plot(stats2(i).Centroid(1),stats2(i).Centroid(2),'b.');
end
hold off
end
