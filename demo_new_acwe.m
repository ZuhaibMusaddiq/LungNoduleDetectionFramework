function u = demo_new_acwe(Img, iterNum)
%Effect: demonstrate the active contour without edge algorithm
%Inputs:
%Img: the input image
%iterNum: number of iteration
%Outputs:
%u: the final level-set function phi
%Author: Su dongcai at 2012/1/12 Email: suntree4152@gmail.com, qq:272973536
Img = double(Img(:, :, 1));
%apply median filter to denoise
Img = medfilt2(Img, [5, 5]);
%setting the initial level set function 'u':
truncated_len = 50;
c0=2;
u = ones(size(Img, 1), size(Img, 2))*c0;
u([truncated_len:end-truncated_len], [truncated_len:end-truncated_len])=-c0; 
%setting the parameters in ACWE algorithm:
mu=1;
lambda1=1; lambda2=1;
timestep = .1; v=1; epsilon=1;
%show the initial 0-level-set contour:
figure;imshow(Img, []);hold on;axis off,axis equal
title('Initial contour');
[c,h] = contour(u,[0 0],'r');
pause(0.1);
I=Img;
count=1;

    [thr I1]=maxentropie(I);
    c11=mean(Img(find(Img>=thr)));
    c22=mean(Img(find(Img<thr)));
% start level set evolution

    u =acwe(u,Img,count,iterNum,);
         
   
    

imshow(Img, []);hold on;axis off,axis equal
[c,h] = contour(u,[-1 -1],'r');


figure;
imagesc(u);axis off,axis equal;
title('Final level set function');
end