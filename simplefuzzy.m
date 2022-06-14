function [Area_image, area_list] = simplefuzzy(in, grading, c, i)
% Here the input parameters are-
%in- the directory or name of the image 
%grading which considers the index of the image if it is found to be in a set 
% c which specifies the number of clusters to initialize the fuzzy clustering algorithm.
% i represents the index of the membership map to be used for extracting
% quantitative data from the image.


% Querying whether image exists in a set or on its own.
if grading == 0
    im = imread(in);
else
    im = imread(in,grading);
end

% This function generates membership matrix U5 and an image histogram H5
[~,U,~,H] = FastFCMeans(im,c);

% This function generates membership maps for all the images. These
% membership maps can be visualized as images
Umap = FM2map(im,U,H);

% Umap has the membership map of interest. This is first processed using an imclose operation and then binarized for extracting
% the clusters easier.
s = strel('disk',2);
Ie = Umap(:,:,i);
Ie = imopen(Ie,s);
Ie = imbinarize(Ie);

% Calculates the centroids of each cluster that was identified
Ccss = bwconncomp(Ie); % image with 10 cluster points removed- this needs to be processed
rpss = regionprops(Ccss, 'Centroid');
centroids = cat(1,rpss.Centroid);

%Visualization of centroid regions within the image.
f1 = figure('name', 'Plotting of the centroid regions on the original image')
imshow(im); hold on;
plot(centroids(:,1), centroids(:,2), '*r');

Area_image = 0;
area_list = 0;

%% Area Extraction
% Step 1- Generate an image object and convert centroid coordinates to
% pixel coordinates within the image.

%Step 1.5- Adjust for 0 values around the boundary of the cell within the
%image.

% Step 2- Store the converted centroid pixel coordinates in an array. Assign each centroid an index number. These indices will mark the
% array that contain distances of the centroids from the selected pixel of
% interest. The centroid index number would remain consistent across all
% pixel values. 3 arrays will thus be generated in total- two with the
% centroid x-y coordinates and one with the index locations. it will be a
% n-by-3 cell matrix as the index will be stored as an integer and the
% coordinates will be stored in double format.

% Step 3- Calculate the distance between the pixels and the centroid
% location and store it in an array- the order in which the distances will
% be calculated will be dictated by the centroid index. The array will be
% termed as the optimal distance array.

%Step 4- The optimal distance (OD) array will then be utilized to determine the
%array index of the minimum distance generated and whether this distance is
%repeated twice in the OD array. Use sum and logical operators for this. If
%it is repeated twice, the centroid with the lower index will be
%responsible for that pixel.

%Step 5- Once the minimum distance is determined, assign the pixel the
%index of the OD array. Repeat for all pixels in the array. 


%This step is conducted to ensure that the image object generated can be used to query pixel coordinates from axes coordinates.
Im1 = image(Ie);

[nrows, ncols] = size(get(Im1,'CData'));
XData = get(Im1, 'XData');
YData = get(Im1, 'YData');

icy = size(Ie,1);
icx = size(Ie,2);

Area_image = zeros(size(Ie));

% Used to obtain the number of rows contained within the data container 'Centroids'

Sz = size(centroids,1);

% Generating an empty array with the same size as the array centroids.

pixel_centroids(Sz,2) = 0;

%Converting centroid axes positions into centroid pixel positions- this step outputs the same results and is thus irrelevant.
for k = 1:Sz
    Kx = centroids(k,1);
    Ky = centroids(k,2);
    pixel_centroids(k,1) = axes2pix(ncols, XData, Kx);
    pixel_centroids(k,2) = axes2pix(nrows, YData, Ky);
end


%An additional data container termed as a cell array to organiye the x and y coordinates of the 
%centroids within the image along with the labels stored in the 3rd cell of the array
OD_array= {};
OD_array{1} = pixel_centroids(:,1);
OD_array{2} = pixel_centroids(:,2);

Sz = size(pixel_centroids,1);

lsz = [1:Sz];
tsz = lsz';

OD_array{3} = tsz; % labels of each of the centroid coordinates stored in the 3rd cell of the OD_array

D_array(Sz) = 0;

% Iteration throughout all the pixels of the image. Exemplified by the double for loop statements.
for i = 1:icx
    for j = 1:icy
        if(Ie(j,i) == 0)
            continue
        end
        
        %Distance of particular pixel in question is assayed against each centroid stored in the earlier cell array. All pixel-centroid distances
        % are stored in another container termed D_array.
        
        for ik = 1:Sz
        D = sqrt((i-OD_array{1}(ik))^2+(j-OD_array{2}(ik))^2);
        D_array(OD_array{3}(ik)) = D;
        
        end
        
        %The minimum value of the d_array is assayed and the corresponding label of this particular minimum value is used as an entry instead
        %of the original greyscale value corresponding to the pixel of the image in question.
        mind = min(D_array);
        binD = D_array == mind;
        duplicates = sum(binD, 'all');
        ind = find(binD);
        
        Area_image(i,j) = ind(1);
        
    end
end

%The image Area_image2 is now a label map that can be used to query the areas of each individual centroid region within the image.

area_list = [];
Area_image2 = [];

%This for loop statement is utilized to extract the areas from individual regions within the Area_image2 label map and store it in a list container termed
%as area_list which gives an iteration of all the labels and their associated areas.
for ji = 1:Sz
   area_list(ji,1) = ji;
   Area_image2 = Area_image;
   Area_image2(Area_image2 > ji) = 0;
   Area_image2(Area_image2 < ji) = 0;
   Area_image2(Area_image2 == ji) = 1;
   area_list(ji,2) = sum(Area_image2,'all');
end


