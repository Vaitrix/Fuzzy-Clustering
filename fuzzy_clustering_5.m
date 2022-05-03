function [C,Bsp, Fimage, Bo, border_mask, ref_image] = fuzzy_clustering_5(in, grading)
% grading stands for the index of the particular image in the picture set.
%ref_image represents channel clusters with the channel boundaries removed.
%border mask shows binarized region of the channel with white representing
%regions where clusters are analyzed

if grading == 0
    im = imread(in);
else
    im = imread(in,grading);
end

[~,U,~,H] = FastFCMeans(im,5);

Umap = FM2map(im,U,H);
figure('color', 'w');
imshow(Umap(:,:,2)); %% fuzzy clustering based on membership map.
ttl = sprintf('Class %d membership map', 1);
set(get(gca, 'Title'), 'String', ttl);


%% Image Filtering- boundary approach


%Planning to do a closing operation to dilate the image boundaries to
%remove it from the screen

SE = strel('disk', 2);
Ie = Umap(:,:,2);
IeL = imbinarize(Ie);

%Having binarized the image, it can now be subjected to filling conditions-
%logical does not work for clustered images because every pixel has a
%certain value and thus will all point to 1. 
Ie1 = imopen(IeL, SE);
Ie1 = background_mask(Ie1,0,1,'all');



f3 = figure('name', 'Dilation_boundary');
figure(f3);
imshow(Ie1);

%%If Boundary is established to be correct- proceed to omit boundary from
%%the clustered image for seed centers.
%Ie = Umap(:,:,2);
%Iec = 1-Ie;
%for j = 1:size(boundary,1)
    %Ie = imfill(Ie, [boundary(j,2), boundary(j,1)],8);
%end

%f3 = figure('name', 'Xenopus droplet without the cell boundary');
%figure(f3);subplot(1,2,1);
%imshow(Ie);
%figure(f3);subplot(1,2,2);
%imshow(Iec);

% Border extraction was unsuccessful because the contours captured the
% region located slightly outside the pixellated boundaries. Turns out the
% boundaries dont have to be removed for the conncomp function to work so
% trying to directly store centroid datapoints.

%Centroids cannot be calculated without the removal of the boundary region.
%So boundary region needs to be removed again.

%binary size- pixels smaller than a particular value are removed

f5 = figure('name', 'Membership maps of the boundary+outside region  of the channel');
figure(f5);

%Update- boundary region from membership map 4 is not general enough to
%spread across different aspect ratios- its time to generate a 2nd fuzzy
%clustering map using only 2 clusters- this would iterate through the
%protocells well.
[~,U,~,H] = FastFCMeans(im,2);
boundary_map = FM2map(im,U,H);


Uc = boundary_map(:,:,2);
subplot(4,2,1)
imshow(Uc);
t1 = sprintf('2nd Membership map');
set(get(gca,'Title'),'String', t1, 'FontSize', 15);
Uc = imbinarize(Uc);
subplot(4,2,2)
imshow(Uc);
t1 = sprintf('Image after Binarization');
set(get(gca,'Title'),'String', t1);
%Ud = (1-Umap(:,:,5))-Umap(:,:,4);
%Ud = imbinarize(Ud);
SE3 = strel('disk',9);
SE4 = strel('disk',6);
%Uc2 = background_mask(Uc,0,0, 'L-R');
Uc2 = Uc;
%Uc21 = Ud;
%Uc2 = imdilate(Uc2, SE3);
%Uc2 = imdilate(Uc2,SE3);
%Uc2 = imdilate(Uc2, SE3);
Uc2 = 1-Uc2;
subplot(4,2,3)
imshow(Uc2);
t1 = sprintf('Image after inversion');
set(get(gca,'Title'),'String', t1);
Uc2 = imfill(Uc2);
subplot(4,2,4)
imshow(Uc2);
t1 = sprintf('Image after filling operation');
set(get(gca,'Title'),'String', t1);
%ff01 = figure;
%figure(ff01);
%imshow(Uc);
%subplot(1,7,1);
%imshow(Uc);

%subplot(1,7,2);
%imshow(Uc2);

%uc2_colsize = size(Uc2,2);
%Uc2 = Uc2(:,4:(uc2_colsize-3));
Ccp = bwconncomp(Uc2);
rp1 = regionprops(Ccp, 'Area');
indB = find([rp1(:).Area] == max([rp1(:).Area]));
Masked_image = zeros(size(Uc2));
Masked_image(Ccp.PixelIdxList{indB}) = 1;
%Masked_image = 1-Masked_image;

%Border_image = zeros(size(Uc21));
%Ccp2 = bwconncomp(Uc21);
%rp12 = regionprops(Ccp2, 'Area');
%indB2 = find([rp12(:).Area] == max([rp12(:).Area]));
%Border_image(Ccp2.PixelIdxList{indB2}) = 1;
%Border_image = 1-Border_image;


%figure(f5);
%subplot(1,7,3);
%imshow(Masked_image);

SE2 = strel('disk',2);
Masked_image2 = imdilate(Masked_image, SE2);

figure(f5);
subplot(4,2,5);
imshow(Masked_image2);
t1 = sprintf('Image after Boundary Extraction+Image Dilation');
set(get(gca,'Title'),'String', t1, 'FontSize', 15);

Masked_image2 = background_mask(Masked_image2,0,0, 'L-R');
Masked_image2 = imerode(Masked_image2, SE2);

figure(f5);
f5.Position(3:4)
subplot(4,2,6);
imshow(Masked_image2);
t1 = sprintf('Image after background mask(L,R)+erosion');
set(get(gca,'Title'),'String', t1);

Masked_image2 = background_mask(Masked_image2, 0,0, 'T-B');
Masked_image2 = imerode(Masked_image2, SE2);

figure(f5);
subplot(4,2,7);
imshow(Masked_image2);
t1 = sprintf('Image after background mask(T,B)+erosion');
set(get(gca,'Title'),'String', t1);

Masked_image2 = imerode(Masked_image2, SE4);

figure(f5);
subplot(4,2,8);
imshow(Masked_image2);
t1 = sprintf('Image after erosion');
set(get(gca,'Title'),'String', t1);





 %ff5 = figure('name', 'Undilated border of the channel');
 %figure(ff5);
 %imshow(Border_image);

Bimage = Ie1.*Masked_image2;
Bimage = imopen(Bimage, SE);

Bsp = Masked_image2;




f6 = figure('name', 'protocell centroid positions before/after boundary mask');
figure(f6);
subplot(1,2,1);
imshow(Ie1);
subplot(1,2,2);
imshow(Bimage);

Cc = bwconncomp(Bimage);
Sd2 = regionprops(Cc, 'centroid');
Total_centroids = cat(1,Sd2.Centroid);

f4 = figure('name', 'Determining cell center locations before filling'); 
figure(f4);
imshow(Bimage);
hold on;
if(isempty(Total_centroids))
    msg = 'The program is unable to identify distinct protocell clusters!';
    error(msg);
else
plot(Total_centroids(:,1), Total_centroids(:,2), 'r.');
end

f5 = figure('name', 'Cell filling');
Fimage = imfill(Bimage);
%figure(f5);
%imshow(Fimage);

Fcc = bwconncomp(Fimage);
Sfd = regionprops(Fcc, 'centroid');
Total_filled_centroids = cat(1,Sfd.Centroid);

f6 = figure('name', 'Determining cell center locations after filling');
%figure(f6);
%imshow(Fimage); hold on;
%plot(Total_filled_centroids(:,1), Total_filled_centroids(:,2), 'r.');

C = Total_centroids;

%% Get channel boundary- circular channels only

Ccbd = bwconncomp(Masked_image2);
rpbd = regionprops(Ccbd, 'Area');
indd = find([rpbd(:).Area] == max([rpbd(:).Area]));
New_masked_image = zeros(size(Masked_image2));
New_masked_image(Ccbd.PixelIdxList{indd}) = 1;
bboundary = bwboundaries(New_masked_image);
bboundary_index = size(bboundary,1);

lb = [];

for k = 1:bboundary_index
    bb = bboundary{k};    
    lb(k) = size(bb,1);    
end

[G,gi] = sort(lb, 'descend');
ind = gi(1);

Bo = bboundary{ind};
SE4 = strel('diamond',5);
New_masked_image = imerode(New_masked_image,SE4);
border_mask = New_masked_image;

%f7 = figure('name', 'Boundary of the channel');
%figure(f7);
%subplot(1,2,1);
%imshow(im); hold on;
%plot(Bo(:,2),Bo(:,1), '.r');
%subplot(1,2,2);

%imshow(New_masked_image);
New_masked_image = uint8(New_masked_image);


im = background_mask(im,0,1, 'all');
ref_image = im.*New_masked_image;


%f8 = figure('name', 'refined image after boundary removal');
%figure(f8);
%imshow(ref_image);
















% All Bead centers are successfully captured and now it is time to
% calculate the area of the protocell clusters around the bead centers.














    



 
 
 


