clear all;
vid = videoinput('macvideo');
set(vid, 'ReturnedColorSpace', 'RGB');
triggerconfig(vid, 'manual');
% start(vid);


background = getsnapshot(vid);
% background = im2double(background);
% background = rgb2gray(background);
% background = background./max(max(background));
imshow(background);