%% IIP Project - (14UEC080, 14UEC109)
%% Cartoonification of an Image
% Main Code

clear all;
close all;

% img3 = double(imread('academy.jpg'))/255;
% img3 = double(imread('einstein.jpg'))/255;
% img3 = double(imread('marissa.jpg'))/255;
% img3 = double(imread('taylor.jpg'))/255;
img3 = double(imread('katiehopkins.jpg'))/255;
cartoon_img3 = cartoon(img3);

% Display color input image and abstracted output.
figure; clf;
set(gcf,'Name','Image Abstraction Input');
subplot(1,2,1);
imagesc(img3); axis image;
title('Input Image');
set(gcf,'Name','Result of Image Abstraction');
subplot(1,2,2);
imagesc(cartoon_img3); axis image;
title('Abstracted Image');