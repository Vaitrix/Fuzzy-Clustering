function [Area_image, area_list] = optifuzzy(in, grading,j)

% j represents the factor by which clusters need to be compressed
% grading represents the index of the picture in the picture set to be
% processed
% in represents the picture set to be processed

 %% Generate points for 5-cluster membership maps
 
[C,Bsp, Fimage, ~, border_mask,ref_image] = fuzzy_clustering_5(in, grading);
 
 if grading == 0
    im = imread(in);
 else
    im = imread(in,grading);
 end

[~,U,~,H] = FastFCMeans(ref_image,10);
[~,U5,~,H5] = FastFCMeans(ref_image,5);

Umap = FM2map(ref_image,U,H);
Umap5 = FM2map(ref_image, U5, H5);

SE = strel('disk', 2);
Ie = Umap(:,:,2);
figure;
imshow(Ie);
IeL = imbinarize(Ie);

Ie1 = imclose(IeL, SE); %changing operation from opening to closing
Ie1 = background_mask(Ie1,0,1, 'all');

%Bimage = Ie1.*Bsp;
Bimage = Ie1;
Bimage = imopen(Bimage, SE); 

ff1 = figure('name', 'protocell centroid positions before/after boundary mask for 10 clusters');
figure(ff1);
subplot(1,2,1);
imshow(Ie1);
subplot(1,2,2);
imshow(Bimage);
%functional to this point
Cc = bwconncomp(Bimage);
Rp11 = regionprops(Cc, 'centroid');
Total_centroids_10 = cat(1,Rp11.Centroid);

%ff05 = figure('name', 'membership maps of all the centroid regions');
ff1p5 = figure('name', 'All operations to generate binarized 5 cluster membership map');
figure(ff1p5)
subplot(2,2,1);
U5 = Umap5(:,:,2);
imshow(U5);
t1 = sprintf('2nd membership map generated from 5 cluster fuzzy clustering');
set(get(gca,'Title'), 'String', t1);
subplot(2,2,2);
U5B = imbinarize(U5);
imshow(U5B);
t1 = sprintf('Image after Binarization');
set(get(gca,'Title'), 'String', t1);
subplot(2,2,3);
U5B = imclose(U5B, SE);
imshow(U5B);
t1 = sprintf('Image after closing operation');
set(get(gca,'Title'), 'String', t1);
subplot(2,2,4);
U5B = imopen(U5B,SE);
imshow(U5B);
t1 = sprintf('Image after opening operation');
set(get(gca,'Title'), 'String', t1);


Ccss = bwconncomp(U5B); % image with 10 cluster points removed- this needs to be processed
rpss = regionprops(Ccss, 'Centroid');
centroid55c = cat(1,rpss.Centroid);

ff2 = figure('name', 'Centroid locations for 5 clusters on a 5 cluster membership map');
figure(ff2);
subplot(1,2,1);
% imshow(Bimage); hold on %Change image with membership map instead
imshow(Umap5(:,:,2)); hold on
plot(centroid55c(:,1), centroid55c(:,2), '*b', 'MarkerSize', 10); 
subplot(1,2,2);
imshow(U5B); hold on
plot(centroid55c(:,1), centroid55c(:,2), '*b', 'MarkerSize', 10); 
%plot(C(:,1), C(:,2), '.r');

%Generated a 2nd set of centroids focused on the 5 cluster membership maps.

ff3 = figure('name', 'Centroid locations of 10 clusters on a 10 cluster membership map');
figure(ff3);
subplot(1,2,1);
%imshow(Fimage); hold on %Change image with membership map instead
imshow(Bimage); hold on 
if(~isempty(Total_centroids_10))
    plot(Total_centroids_10(:,1), Total_centroids_10(:,2), '.r', 'MarkerSize', 12);
end
subplot(1,2,2);
imshow(Umap(:,:,2)); hold on
if(~isempty(Total_centroids_10))
    plot(Total_centroids_10(:,1), Total_centroids_10(:,2), '.r', 'MarkerSize', 12);
end

% Removing areas with associated 10 membership clusters from 5-cluster maps
[g h] = size(Total_centroids_10);
[A b] = bwlabel(U5B); %change from Fimage to U5B


for ii = 1:g
    x5 = ceil(Total_centroids_10(ii,1));
    y5 = ceil(Total_centroids_10(ii,2));
    
    if(A(y5, x5) > 0)
        label = A(y5, x5);
        A(A == label) = 0;
    end
end

A(A>0) = 1;

S_image = A;


ff4 = figure('name', '5 cluster membership map before/after removal of 10-cluster areas');
figure(ff4);
subplot(1,2,1);
imshow(U5B);hold on; plot(Total_centroids_10(:,1), Total_centroids_10(:,2), '*r', 'MarkerSize', 10); hold on; %changed image plotted from Fimage to U5B
subplot(1,2,2);
imshow(S_image); hold on; plot(Total_centroids_10(:,1), Total_centroids_10(:,2), '*r', 'MarkerSize', 10);

Ccs = bwconncomp(S_image); % image with 10 cluster points removed- this needs to be processed
rps = regionprops(Ccs, 'Centroid');
centroid5c = cat(1,rps.Centroid);

centroid_list = cat(1,centroid5c, Total_centroids_10);

[~,U5,~,H5] = FastFCMeans(im,5);

U5map = FM2map(im, U5,H5);

U5b = imbinarize(U5map(:,:,1));

%ff5 = figure('name', 'Final list of centroids captured on a membership map');
%figure(ff5);
%imshow(im); hold on;
%plot(centroid_list(:,1), centroid_list(:,2), '.r');

unclipped_centroids = centroid_list;

%% Calculating whether clusters need to be compressed or not

[lt, bd] = size(centroid_list);

compressed_list_x = [];
compressed_list_y = [];
Dk = j;
for ki = 1:lt
    compressed_list_x = [];
    compressed_list_y = [];
    Dk = j;
    for li = 1:lt
        D = sqrt((centroid_list(ki,1)-centroid_list(li,1))^2+(centroid_list(ki,2)-centroid_list(li,2))^2);
        Dk = D;
        if(D < j)&&(li ~= ki)
            compressed_list_x(1,li) = centroid_list(li,1);
            compressed_list_y(1,li) = centroid_list(li,2);
            centroid_list(li,1) = 0;
            centroid_list(li,2) = 0;
        end
    end
    if(Dk < j)&&(Dk ~= 0)
    compressed_list_x = compressed_list_x(compressed_list_x > 0); %error revised
    compressed_list_y = compressed_list_y(compressed_list_y > 0); %error revised
    Mx = mean(compressed_list_x);
    My = mean(compressed_list_y);
    centroid_list(ki,1) = Mx;
    centroid_list(ki,2) = My;
    end
 
end

centroid_list_x = centroid_list(:,1);
centroid_list_y = centroid_list(:,2);
centroid_list_x = centroid_list_x(centroid_list_x > 0);
centroid_list_y = centroid_list_y(centroid_list_y > 0);

clipped_centroids= cat(2, centroid_list_x, centroid_list_y);

% clipped_centroids = cat(1, clipped_centroids,Bo);

ff6 = figure('name', 'Final list of centroids after compression captured on a membership map');
figure(ff6);
subplot(1,2,1);
imshow(im); hold on;
plot(unclipped_centroids(:,1), unclipped_centroids(:,2), '.r')
subplot(1,2,2);
imshow(im); hold on;
plot(clipped_centroids(:,1), clipped_centroids(:,2), '.g');


%% Plot voronoi tessellation- approach 1

 [v,c] = voronoin([clipped_centroids(:,1),clipped_centroids(:,2)]);
 
  ensind=[];% All the labels of the cells inside the droplet
  figure(ff6); 
  A=[];
  
  for j = 1:length(c)
      v1 = v(c{j},1) ;
      v2 = v(c{j},2) ;
        % Test to check if current Voronoi cell is in the droplet
      if length(find(in==0))==0
          ensind=[ensind,j];
          plot([v1;v1(1)],[v2;v2(1)],'-r','linewidth',2)
          'cell inside';
          A_loc =polyarea(v1,v2);
          A=[A,A_loc] ;
          %pause
       end
  end

   
%% Plot voronoi tessellation- approach 2
im = im2double(im);
border_mask = im2double(border_mask);
im = background_mask(im,0,1, 'all');
border_image = im.*border_mask;

ff7 = figure('name', 'voronoi tessellation- 2nd approach');
figure(ff7);
imshow(border_image); hold on;
voronoi(clipped_centroids(:,1), clipped_centroids(:,2));


%% Calculating the areas of the Voronoi segmentations 
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
Im1 = image(border_image);

[nrows, ncols] = size(get(Im1,'CData'));
XData = get(Im1, 'XData');
YData = get(Im1, 'YData');

icy = size(border_image,1);
icx = size(border_image,2);

Area_image = zeros(size(border_image));

Sz = size(clipped_centroids,1);

for k = 1:Sz
    Kx = clipped_centroids(k,1);
    Ky = clipped_centroids(k,2);
    pixel_centroids(k,1) = axes2pix(ncols, XData, Kx);
    pixel_centroids(k,2) = axes2pix(nrows, YData, Ky);
end

OD_array= {};
OD_array{1} = pixel_centroids(:,1);
OD_array{2} = pixel_centroids(:,2);

Sz = size(pixel_centroids,1);

lsz = [1:Sz];
tsz = lsz';

OD_array{3} = tsz;

D_array(Sz) = 0;

for i = 1:icx
    for j = 1:icy
        if(border_image(j,i) == 0)
            continue
        end
        
        for ik = 1:Sz
        D = sqrt((i-OD_array{1}(ik))^2+(j-OD_array{2}(ik))^2);
        D_array(OD_array{3}(ik)) = D;
        
        end
        
        mind = min(D_array);
        binD = D_array == mind;
        duplicates = sum(binD, 'all');
        ind = find(binD);
        
        Area_image(i,j) = ind(1);
        
    end
end


area_list = [];
Area_image2 = [];

for ji = 1:Sz
   area_list(ji,1) = ji;
   Area_image2 = Area_image;
   Area_image2(Area_image2 > ji) = 0;
   Area_image2(Area_image2 < ji) = 0;
   Area_image2(Area_image2 == ji) = 1;
   area_list(ji,2) = sum(Area_image2,'all');
end

figure(ff7);
imshow(border_image); hold on;
[vx, vy]= voronoi(clipped_centroids(:,1), clipped_centroids(:,2));
plot(clipped_centroids(:,1), clipped_centroids(:,2), 'r*',vx,vy,'b-', 'linewidth', 2); axis equal

end

        

    





























 
 
 
 
 
