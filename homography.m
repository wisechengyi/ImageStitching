close all;
clear all;
% CS 543 Assignment 1, starter Matlab code
% Adapted from A. Efros
% (http://graphics.cs.cmu.edu/courses/15-463/2010_fall/hw/proj1/)

FEATURENUM=200;
BINSIZE=15;
PATCHSIZE=5;
PairNum=100;
RANSAC_ITERATION=5000;
Threshold = 5;
SAMPLESIZE=6; % # of pairs at a time to compute the parameters

% read in the image
% fullimLeft = imread('left.jpg');
% fullimRight = imread('right.jpg');
fullimLeft = imread('set3A.jpg');
fullimRight = imread('set3B.jpg');
% fullimRight = imresize(fullimRight,[reRow reCol]);

% fullimLeft = imread('stitch_01.png');
% fullimRight = imread('stitch_02.png');

% convert to double matrix (might want to do this later on to same memory)
fullimLeft = im2double(fullimLeft);
fullimRight = im2double(fullimRight);

greyImageLeft = rgb2gray(fullimLeft);
greyImageLeft = im2double(greyImageLeft);

greyImageRight = rgb2gray(fullimRight);
greyImageRight = im2double(greyImageRight);

maxLeft=max(max(greyImageLeft));
maxRight=max(max(greyImageRight));

% % % normalization
greyImageLeft = greyImageLeft./maxLeft;
greyImageRight = greyImageRight./maxRight;

% cornerLeft = corner(greyImageLeft,'QualityLevel',.0005);
% cornerRight = corner(greyImageRight,'QualityLevel',.0005);

% corner feature extraction, 1st column is X, 2nd column is Y
cornerLeft = corner(greyImageLeft,FEATURENUM);
cornerRight = corner(greyImageRight,FEATURENUM);



% for every corner in each image, compute the histogram of the pixel values
% in bins as the descriptor
leftDescriptorArray = getDescriptor(greyImageLeft,cornerLeft,BINSIZE,PATCHSIZE);
rightDescriptorArray = getDescriptor(greyImageRight,cornerRight,BINSIZE,PATCHSIZE);

% calculate the dot product between the corner features into an 1-D array
sizeRight = size(cornerRight,1);
sizeLeft = size(cornerLeft,1);
correlation = zeros(sizeLeft*sizeRight,1);
index=1;
for i=1:sizeLeft
    for j=1:sizeRight
        %         correlation(index) = dot(leftDescriptorArray{i},rightDescriptorArray{j});
        correlation(index) = rms(leftDescriptorArray{i}-rightDescriptorArray{j});
        index=index+1;
    end
end

% the getting top PairNum matches
[sortedValues,sortIndex] = sort(correlation,'ascend');  %# Sort the values in ascending order
minIndex = sortIndex(1:PairNum); %# Get a linear index into A of the 5 smallest values

correspondingCornerIndices = zeros(PairNum,2);
for n=1:PairNum
    pairIndex = minIndex(n);
    leftCornerIndex = floor((pairIndex-1) / sizeRight)+1;
    rightCornerIndex = mod(pairIndex-1,sizeRight)+1;
    correspondingCornerIndices(n,:)= [leftCornerIndex rightCornerIndex];
end

% put two images side by side in order to draw lines between correponding
% corners
[leftImageHeight, leftImageWidth] = size(greyImageLeft);
[rightImageHeight, rightImageWidth] = size(greyImageRight);
totalImage=zeros(max(leftImageHeight,rightImageHeight),leftImageWidth+rightImageWidth);
totalImage(1:leftImageHeight,1:leftImageWidth)=greyImageLeft;
totalImage(1:rightImageHeight,leftImageWidth+1:end)=greyImageRight;

imshow(totalImage);
hold on;
plot(cornerLeft(:,1), cornerLeft(:,2), 'r*');
plot(cornerRight(:,1)+leftImageWidth, cornerRight(:,2), 'r*');
for n=1:PairNum
    CornerPairLeftIndex = correspondingCornerIndices(n,1);
    CornerPairRightIndex = correspondingCornerIndices(n,2);
    leftX=cornerLeft(CornerPairLeftIndex,1);
    leftY=cornerLeft(CornerPairLeftIndex,2);
    rightX=cornerRight(CornerPairRightIndex,1);
    rightY=cornerRight(CornerPairRightIndex,2);
    %     if ( rightX<rightImageWidth/2)
    plot([leftX rightX+leftImageWidth],[leftY rightY],'Color','g');
    %     else
    %          plot([leftX rightX+leftImageWidth],[leftY rightY],'Color','g');
    %     end
end

% run RANSAC for affline model
% first, pick 3 pairs of corners, fit it into the equation on page 13 of
% lecture 12

% obtainedTransformations stores the computed m1 m2 m3 m4 t1 t2 in order
obtainedTransformations= cell(RANSAC_ITERATION,1);
% IORatio stores the inlier and outlier ratio for each particular
% Transformation matrix
IORatio = zeros(RANSAC_ITERATION,1);

x1s = cell(RANSAC_ITERATION,1);
y1s = cell(RANSAC_ITERATION,1);
x2s =cell(RANSAC_ITERATION,1);
y2s=cell(RANSAC_ITERATION,1);


for n= 1: RANSAC_ITERATION
    r=randperm(PairNum);
    pairIndices = r(1:SAMPLESIZE);
    
%   http://cseweb.ucsd.edu/classes/wi07/cse252a/homography_estimation/homography_estimation.pdf
    MatrixA=zeros(2*SAMPLESIZE,9);
    theX1=zeros(4,1);
    theY1=zeros(4,1);
    theX2 = zeros(4,1);
    theY2 = zeros(4,1);
    
    for ii=1:SAMPLESIZE
        pairIdx=pairIndices(ii);
        leftCnrIdx=correspondingCornerIndices(pairIdx,1);
        rightCnrIdx=correspondingCornerIndices(pairIdx,2);
        x1=cornerLeft(leftCnrIdx,1);
        y1=cornerLeft(leftCnrIdx,2);
        x2=cornerRight(rightCnrIdx,1);
        y2=cornerRight(rightCnrIdx,2);
        MatrixA(2*ii-1,:)=[-x1 -y1 -1 0   0   0  x2*x1 x2*y1 x2];
        MatrixA(2*ii,:)  =[0   0   0  -x1 -y1 -1 y2*x1 y2*y1 y2];
        if (ii<=4)
        theX1(ii)=x1;
        theY1(ii)=y1;
        theX2(ii)=x2;
        theY2(ii)=y2;        
        end
    end
    
    x1s{n}=theX1;
    y1s{n}=theY1;
    x2s{n}=theX2;
    y2s{n}=theY2;
    
    if ( cond(MatrixA)>1e9)
        obtainedTransformations{n}=[0 0 0 0 0 0 0 0 0];
        IORatio(n) = -1;
        continue;
    end
    
    [U, S, V]=svd(MatrixA);
    P=V(:,9);
%     P=null(MatrixA);
%     m1=P(1);
%     m2=P(2);
%     m3=P(3);
%     m4=P(4);
%     t1=P(5);
%     t2=P(6);
%     RotScale=[m1 m2;m3 m4];
%     Trans = [t1;t2];
    obtainedTransformations{n}=P;
    
    inlier =1;
    outlier=1;
    for i=1:PairNum
        testX=cornerLeft(correspondingCornerIndices(i,1),1);
        testY=cornerLeft(correspondingCornerIndices(i,1),2);
        testXp=cornerRight(correspondingCornerIndices(i,2),1);
        testYp=cornerRight(correspondingCornerIndices(i,2),2);
        H=[P(1) P(2) P(3);P(4) P(5) P(6);P(7) P(8) P(9)];
        HomoPoint = H*[testX;testY;1];
        NonHomoPoint = HomoPoint./HomoPoint(3,1);
        Residual = NonHomoPoint - [testXp;testYp;1];
        dist = norm(Residual,2);
        if (dist>Threshold)
            outlier=outlier+1;
        else
            inlier=inlier+1;
        end
    end
    IORatio(n) = inlier/outlier;
    
end
[Vmax, I] =max(IORatio)
temp = obtainedTransformations{I};


X1=x1s{I};
Y1=y1s{I};
X2=x2s{I};
Y2=y2s{I};
% T=maketform('projective',[X1 Y1],[X2 Y2]);
H=[temp(1) temp(2) temp(3);temp(4) temp(5) temp(6);temp(7) temp(8) temp(9)];
T=maketform('projective',H');

% t = maketform('projective',H);
% J = imtransform(greyImageLeft,t);
% J = imtransform(greyImageLeft,T);
% figure;
% imshow(J);
im2 = greyImageLeft;
im1 = greyImageRight;

[im2t,xdataim2t,ydataim2t]=imtransform(im2,T);
% now xdataim2t and ydataim2t store the bounds of the transformed im2
xdataout=[min(1,xdataim2t(1)) max(size(im1,2),xdataim2t(2))];
ydataout=[min(1,ydataim2t(1)) max(size(im1,1),ydataim2t(2))];
% let's transform both images with the computed xdata and ydata
im2t=imtransform(im2,T,'XData',xdataout,'YData',ydataout);
im1t=imtransform(im1,maketform('projective',eye(3)),'XData',xdataout,'YData',ydataout);
ims=im1t/2+im2t/2;
figure, imshow(ims)



newImage=im1t+im2t;
[row,col]= size(newImage);
for i=1:row
    for j=1:col
        if (im1t(i,j)*im2t(i,j)==0)
            newImage(i,j)=im2t(i,j)+im1t(i,j);
        else
            newImage(i,j)=(im2t(i,j)+im1t(i,j))/2;
        end
       
    end
   
end
figure;
imshow(newImage);

