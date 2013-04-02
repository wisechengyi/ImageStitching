close all;
clear all;
% CS 543 Assignment 3


FEATURENUM=300;
BINSIZE=15;
PATCHSIZE=5;
PairNum=150;
RANSAC_ITERATION=10000;
Threshold = 3;
SAMPLESIZE=5; % # of pairs at a time to compute the parameters

% read in the image
fullimLeft = imread('left.jpg');
fullimRight = imread('right.jpg');
% fullimLeft = imread('a.png');
% fullimRight = imread('b.png');

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
    
    MatrixA=zeros(2*SAMPLESIZE,6);
    MatrixC=zeros(2*SAMPLESIZE,1);
    
    theX1=zeros(3,1);
    theY1=zeros(3,1);
    theX2 = zeros(3,1);
    theY2 = zeros(3,1);
    
    for ii=1:SAMPLESIZE
        pairIdx=pairIndices(ii);
        leftCnrIdx=correspondingCornerIndices(pairIdx,1);
        rightCnrIdx=correspondingCornerIndices(pairIdx,2);
        x1=cornerLeft(leftCnrIdx,1);
        y1=cornerLeft(leftCnrIdx,2);
        x1p=cornerRight(rightCnrIdx,1);
        y1p=cornerRight(rightCnrIdx,2);
        MatrixA(2*ii-1,1:6)=[x1 y1 0 0 1 0];
        MatrixA(2*ii,1:6)=[0 0 x1 y1 0 1];
        MatrixC(2*ii-1)=x1p;
        MatrixC(2*ii)=y1p;
        if (ii<=3)
            theX1(ii)=x1;
            theY1(ii)=y1;
            theX2(ii)=x1p;
            theY2(ii)=y1p;
        end
    end
    
    x1s{n}=theX1;
    y1s{n}=theY1;
    x2s{n}=theX2;
    y2s{n}=theY2;
    
    if ( cond(MatrixA)>1e7)
        obtainedTransformations{n}=[0 0 0 0 0 0];
        IORatio(n) = -1;
        continue;
    end
    P=MatrixA\MatrixC;
    m1=P(1);
    m2=P(2);
    m3=P(3);
    m4=P(4);
    t1=P(5);
    t2=P(6);
    RotScale=[m1 m2;m3 m4];
    Trans = [t1;t2];
    obtainedTransformations{n}=[m1 m2 m3 m4 t1 t2];
    
    inlier =1;
    outlier=1;
    for i=1:PairNum
        testX=cornerLeft(correspondingCornerIndices(i,1),1);
        testY=cornerLeft(correspondingCornerIndices(i,1),2);
        testXp=cornerRight(correspondingCornerIndices(i,2),1);
        testYp=cornerRight(correspondingCornerIndices(i,2),2);
        Residual = RotScale*[testX;testY]+Trans - [testXp;testYp];
        dist = norm(Residual,2);
        if (dist>Threshold)
            outlier=outlier+1;
        else
            inlier=inlier+1;
        end
    end
    IORatio(n) = inlier/outlier;
    
end
[V I] =max(IORatio)
AffineMatrix = obtainedTransformations{I};
m1=AffineMatrix(1);
m2=AffineMatrix(2);
m3=AffineMatrix(3);
m4=AffineMatrix(4);
t1=AffineMatrix(5);
t2=AffineMatrix(6);
M=[m1 m2;m3 m4]
T=[t1 ; t2]
% transpose the rotation matrix
t = maketform('affine',[m1 m3 0;m2 m4 0;t1 t2 1]);
% J = imtransform(greyImageLeft,t);
% offsetX=floor(t1);
% offsetY=floor(t2);
% offsetX=offsetX
% newImage=zeros(size(J)+size(greyImageRight));
% newImage(1:size(J,1),1:size(J,2))=J;
% row=1;
% for i=-offsetY:-offsetY+rightImageHeight-1
%     col=1;
%     for j=-offsetX:-offsetX+rightImageWidth-1
%         if (newImage(i,j)~=0)
%             newImage(i,j)=(greyImageRight(row,col)+newImage(i,j))/2;
%         else
%             newImage(i,j)=greyImageRight(row,col);
%         end
%         col=col+1;
%     end
%     row=row+1;
% end
% figure;
% imshow(newImage);


% X1=x1s{I};
% Y1=y1s{I};
% X2=x2s{I};
% Y2=y2s{I};
% T=maketform('affine',[X1 Y1],[X2 Y2]);
im2 = greyImageLeft;
im1 = greyImageRight;

[im2t,xdataim2t,ydataim2t]=imtransform(im2,t);
% now xdataim2t and ydataim2t store the bounds of the transformed im2
xdataout=[min(1,xdataim2t(1)) max(size(im1,2),xdataim2t(2))];
ydataout=[min(1,ydataim2t(1)) max(size(im1,1),ydataim2t(2))];
% let's transform both images with the computed xdata and ydata
im2t=imtransform(im2,t,'XData',xdataout,'YData',ydataout);
im1t=imtransform(im1,maketform('affine',eye(3)),'XData',xdataout,'YData',ydataout);
ims=im1t/2+im2t/2;
figure, imshow(ims)


