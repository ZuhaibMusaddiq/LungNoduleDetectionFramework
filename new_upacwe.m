function[ u] = new_upacwe(u,c11,c22,Img,orginal,count,iterNum,oldth)
    if(count<iterNum)
        
     u=NeumannBoundCond(u);
     mu=1;
     lambda1=1; 
     lambda2=1;
     timestep = .1; v=1; epsilon=1;
     pc=1;
     DrcU=(epsilon/pi)./(epsilon^2+u.^2);                %eq.(9), ref[2] the delta function
     Hu=0.5*(1+(2/pi)*atan(u./epsilon));                 %eq.(8)[2] the character function how large is 'epsilon'?
     %figure, hist(Hu(:), 100);
     th = .5;
     inside_idx = find(Hu(:) < th); outside_idx = find(Hu(:) >= th);
     c1 = mean(orginal(inside_idx)); c2 = mean(orginal(outside_idx));  %local
     K=curvature_central(u);
     G=K-(c11-c22);
     data_force = -DrcU.*(mu*G - v - lambda1*(orginal-c1).^2 + lambda2*(orginal-c2).^2);
     P=pc*(4*del2(u) - G);               %ref[2]
     u = u+timestep*(data_force+P);
     seg = u<=0;
     [B,L] = bwboundaries(seg); 
      measurements= regionprops(L, 'all');
           for a = 1 : length(measurements)
               if(measurements(a).Area>2)
                thisBB = measurements(a).BoundingBox;
                Img= imcrop(Img,[thisBB(1),thisBB(2),thisBB(3),thisBB(4)]);
                [thr I1]=maxentropie(Img); 
                c11=mean(Img(find(Img>=thr)));
                c22=mean(Img(find(Img<thr)));
                oldth=thr;
                u=new_upacwe(u,c11,c22,Img,orginal,count,iterNum,oldth);
                                    
               end
           end
count=count+1
 pause(0.1);
imshow(orginal, []);hold on;axis off,axis equal
        [c,h] = contour(u,[0 0],'r');
    else
        return;
    end

function g = NeumannBoundCond(f)
%Neumann boundary condition
%originally written by Li chunming
%http://www.mathworks.com/matlabcentral/fileexchange/12711-level-set-for-image-segmentation
[nrow, ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  

function k = curvature_central(u)
%compute curvature:
%originally written by Li chunming
%http://www.mathworks.com/matlabcentral/fileexchange/12711-level-set-for-im
%age-segmentation
[ux, uy] = gradient(u);
normDu = sqrt(ux.^2+uy.^2+1e-10);

Nx = ux./normDu; Ny = uy./normDu;
[nxx, junk] = gradient(Nx); [junk, nyy] = gradient(Ny);
k = nxx+nyy;                                          