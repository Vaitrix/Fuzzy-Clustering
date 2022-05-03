function Bimage = fuzzy_clustering_10(in, grading)

[C,Bsp, Fimage] = fuzzy_clustering_5(in, grading);

if grading == 0
    im = imread(in);
else
    im = imread(in,grading);
end

[~,U,~,H] = FastFCMeans(im,10);

Umap = FM2map(im,U,H);

SE = strel('disk', 2);
Ie = Umap(:,:,2);
IeL = imbinarize(Ie);

Ie1 = imopen(IeL, SE);
Ie1 = background_mask(Ie1,0);

Bimage = Ie1.*Bsp;
Bimage = imopen(Bimage, SE);
