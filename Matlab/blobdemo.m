%------------------------------------------------------------------------------------------------
% Demo to illustrate simple blob detection, measurement, and filtering.
% Requires the Image Processing Toolbox (IPT) because it demonstates some functions
% supplied by that toolbox, plus it uses the "coins" demo image supplied with that toolbox.
% If you have the IPT (you can check by typing ver on the command line), you should be able to
% run this demo code simply by copying and pasting this code into a new editor window,
% and then clicking the green "run" triangle on the toolbar.
% Running time = 7.5 seconds the first run and 2.5 seconds on subsequent runs.
% A similar Mathworks demo:
% http://www.mathworks.com/products/image/demos.html?file=/products/demos/shipping/images/ipexprops.html
% Code written and posted by ImageAnalyst, July 2009.
%------------------------------------------------------------------------------------------------
% function BlobsDemo()
% echo on;
% Startup code.
tic; % Start timer.
clc; % Clear command window.
clear all; % Get rid of variables from prior run of this m-file.
disp('Running BlobsDemo.m...'); % Message sent to command window.
workspace; % Show panel with all the variables.

% Change the current folder to the folder of this m-file.
% (The line of code below is from Brett Shoelson of The Mathworks.)
if(~isdeployed)
cd(fileparts(which(mfilename)));
end

% Read in standard MATLAB demo image
fullFilename = 'lines.bmp';
rgbOriginalImage = imread(fullFilename);
% It's color, so take just the red channel:
originalImage = rgbOriginalImage(:,:,1);
subplot(3, 3, 1);
imshow(originalImage);
% Maximize the figure window.
set(gcf, 'Position', get(0, 'ScreenSize'));
% Force it to display RIGHT NOW (otherwise it might not display until it's all done, unless you've stopped at a breakpoint.)
drawnow;
caption = sprintf('Original "coins" image showing\n6 nickels (the larger coins) and 4 dimes (the smaller coins).');
title(caption);
axis square; % Make sure image is not artificially stretched because of screen's aspect ratio.

% Just for fun, let's get its histogram.
[pixelCount grayLevels] = imhist(originalImage);
subplot(3, 3, 2);
bar(pixelCount); title('Histogram of original image');
xlim([0 grayLevels(end)]); % Scale x axis manually.

% Threshold the image to get a binary image (only 0's and 1's) of class "logical."
% Method #1: using im2bw()
% normalizedThresholdValue = 0.4; % In range 0 to 1.
% thresholdValue = normalizedThresholdValue *
max(max(originalImage)); % Gray Levels.
% binaryImage = im2bw(originalImage, normalizedThresholdValue); % One way to threshold to binary
% Method #2: using a logical operation.
thresholdValue = 100;
binaryImage = originalImage > thresholdValue; % Bright objects will be the chosen if you use >.
% binaryImage = originalImage < thresholdValue; % Dark objects will be the chosen if you use <.
binaryImage = imclearborder(binaryImage, 4);

% Do a "hole fill" to get rid of any background pixels inside the blobs.
binaryImage = imfill(binaryImage, 'holes');

% Show the threshold as a vertical red bar on the histogram.
hold on;
maxYValue = ylim;
hStemLines = stem(thresholdValue, maxYValue(2), 'r');
children = get(hStemLines, 'children');
set(children(2),'visible', 'off');
% Place a text label on the bar chart showing the threshold.
annotationText = sprintf('Thresholded at %d gray levels', thresholdValue);
% For text(), the x and y need to be of the data class "double" so let's cast both to double.
text(double(thresholdValue + 5), double(0.5 * maxYValue(2)), annotationText, 'FontSize', 10, 'Color', [0 .5 0]);
text(double(thresholdValue - 70), double(0.94 * maxYValue(2)), 'Background', 'FontSize', 10, 'Color', [0 0 .5]);
text(double(thresholdValue + 50), double(0.94 * maxYValue(2)), 'Foreground', 'FontSize', 10, 'Color', [0 0 .5]);


% Display the binary image.
subplot(3, 3, 3); imagesc(binaryImage); colormap(gray(256));
title('Binary Image, obtained by thresholding'); axis square;

labeledImage = bwlabel(binaryImage, 8); % Label each blob so we can make measurements of it
coloredLabels = label2rgb (labeledImage, 'hsv', 'k', 'shuffle'); % pseudo random color labels

subplot(3, 3, 4);
imshow(labeledImage, []);
title('Labeled Image, from bwlabel()');
axis square;
subplot(3, 3, 5);
imagesc(coloredLabels);
axis square;
caption = sprintf('Pseudo colored labels, from label2rgb().\nBlobs are numbered from top to bottom, then from left to right.');
title(caption);

% Get all the blob properties. Can only pass in originalImage in
version R2008a and later.
blobMeasurements = regionprops(labeledImage, originalImage, 'all');
numberOfBlobs = size(blobMeasurements, 1);

% bwboundaries() returns a cell array, where each cell contains the row/column coordinates for an object in the image.
% Plot the borders of all the coins on the original grayscale image using the coordinates returned by bwboundaries.
subplot(3, 3, 6); imagesc(originalImage);
title('Outlines, from bwboundaries()'); axis square;
hold on;
boundaries = bwboundaries(binaryImage);
numberOfBoundaries = size(boundaries);
for k = 1 : numberOfBoundaries
thisBoundary = boundaries{k};
plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);
end
hold off;

fontSize = 14;	% Used to control size of "blob number" labels put atop the image.
labelShiftX = -7;	% Used to align the labels in the centers of the coins.
% Print header line in the command window.
fprintf(1,'Blob # Mean Intensity Area Perimeter Centroid Diameter (ECD) Solidity\n');
% Loop over all blobs printing their measurements to the command window.
for k = 1 : numberOfBlobs % Loop through all blobs.
% Find the mean of each blob. (R2008a has a better way where you can pass the original image
% directly into regionprops. The way below works for all versions including earlier versions.)
thisBlobsPixels = blobMeasurements(k).PixelIdxList; % Get list of pixels in current blob.
meanGL = mean(originalImage(thisBlobsPixels)); % Find mean intensity (in original image!)
meanGL2008a = blobMeasurements(k).MeanIntensity; % Mean again, but only for version >= R2008a

blobArea = blobMeasurements(k).Area;	 % Get area.
blobPerimeter = blobMeasurements(k).Perimeter;	 % Get perimeter.
blobCentroid = blobMeasurements(k).Centroid;	 % Get centroid.
blobSolidity = blobMeasurements(k).Solidity;	 % Get solidity.
blobECD = sqrt(4 * blobArea / pi);	 % Compute ECD - Equivalent Circular Diameter.
fprintf(1,'#%2d %17.1f %11.1f %8.1f %8.1f %8.1f % 8.1f %8.1f\n', k, meanGL, blobArea, blobPerimeter, blobCentroid, blobECD, blobSolidity);
% Put the "blob number" labels on the "boundaries" grayscale image.
text(blobCentroid(1) + labelShiftX, blobCentroid(2), num2str(k),  'FontSize', fontSize, 'FontWeight', 'Bold');
end

% Put the labels on the rgb labeled image also.
subplot(3, 3, 5);
for k = 1 : numberOfBlobs % Loop through all blobs.
blobCentroid = blobMeasurements(k).Centroid;	 % Get centroid.
text(blobCentroid(1) + labelShiftX, blobCentroid(2), num2str(k), 'FontSize', fontSize, 'FontWeight', 'Bold');
end

% Now I'll demonstrate how to select certain blobs based using the ismember function.
% Let's say that we wanted to find only those blobs
% with an intensity between 150 and 220 and an area less than 2000 pixels.
% This would give us the three brightest dimes (the smaller coin type).
allBlobSolidities = [blobMeasurements.Solidity];
% Get a list of the blobs that meet our criteria and we need to keep.
allowableIndexes = allBlobSolidities > 0.8; % Take the round objects.
keeperIndexes = find(allowableIndexes);
% Extract only those blobs that meet our criteria, and
% eliminate those blobs that don't meet our criteria.
% Note how we use ismember() to do this.
keeperBlobsImage = ismember(labeledImage, keeperIndexes);
subplot(3, 3, 7);
imshow(keeperBlobsImage, []);
axis square;
title('"Keeper" blobs (round blobs)');
% Re-label with only the keeper blobs kept.
labeledRoundBlobImage = bwlabel(keeperBlobsImage, 4); % Label each blob so we can make measurements of it
% Now we're done. We have a labeled image of blobs that meet our specified criteria.
subplot(3, 3, 8);
imshow(labeledRoundBlobImage, []);
axis square;
title('Re-labeled blobs (round blobs)');

elapsedTime = toc;
% Alert user that the demo is done and give them the option to save an image.
message = sprintf('Finished running BlobsDemo.m.\n\nElapsed time = %. 2f seconds.', elapsedTime);
message = sprintf('%s\n\nCheck out the figure window for the images.\n Check out the command window for the numerical results.', message);