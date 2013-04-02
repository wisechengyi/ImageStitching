T = maketform('affine',[.5 0 0; .5 2 0; 0 0 1]);

I = imread('left.jpg');
I2 = imtransform(I,T);
figure, imshow(I2)

tformfwd([10 20],T)
I3 = imtransform(I,T);
figure
imshow(I3)
