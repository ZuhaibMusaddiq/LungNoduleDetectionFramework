FM='G:\Thesis\Dataset\LIDC-IDRI\LIDC-IDRI-0001\1.3.6.1.4.1.14519.5.2.1.6279.6001.298806137288633453246975630178\1.3.6.1.4.1.14519.5.2.1.6279.6001.179049373636438705059720603192\';

% FM='G:\Thesis\Dataset\LIDC-IDRI\LIDC-IDRI-0001\1.3.6.1.4.1.14519.5.2.1.6279.6001.298806137288633453246975630178\1.3.6.1.4.1.14519.5.2.1.6279.6001.179049373636438705059720603192\'
srcFiles = dir(strcat(FM,'\*.dcm'));
  load('FINAL.mat')
 %close all;
%% Taking the size of the volumatic images for row, column and slices
filename = strcat(FM,srcFiles(1).name);
orig=dicomread(filename);
header=dicominfo(filename);
[r c]=size(image);
features =  struct('slice',[]);
noduleCandidate=zeros(r,c)
close all;

for i=52
   
  
    %% Step 1 Thresholdin
    filename = strcat(FM,srcFiles(i).name);
     info=dicominfo(filename);
      image=dicomread(filename);
      [r c slices]=size(image);
      imageMed=medfilt2(image);
      %sharp=imsharpen(imageMed);
      
      ind=find((imageMed)>500); % find pixels with HU value >500.
      new=ones(r,c);
      new(ind)=0;  % find pixels with HU value >500.

perimImage=~imfill(~new,'holes'); %% region segmentation can be acheived using imfill here same result. region growing was time consuming
P=(r+c)/2;
P=floor(P*0.30);
n = bwareaopen(~new+perimImage,P); %% remove the areas in images that are very small
yu=imerode(~n,strel('disk',2));
v=imfill(yu,'holes'); %% Rolling ball reconstruction of missing edges
L = bwlabel(v);
s = regionprops(L, 'PixelIdxList', 'PixelList','Area','BoundingBox','Eccentricity' );
Im=zeros(r,c);
%Im=Im.*500;
%image= imsharpen(image) ;

    %Lung Extraction
for k=1:length(s)
     s(k).Eccentricity
     thisBB = s(k).BoundingBox;
    if((s(k).Area>1))
           temp2=image;
           temp=image;
           temp(:,:)=0;
           x=s(k).PixelIdxList;
           temp(x)=temp2(x);
           
           Im(x)=image(x);
      
    end
end

Im=imfill(Im,'holes');
  
 
%        mkdir(strcat('SegmentedLungs\',header.PatientID));
%        
%      
%       imwrite(Im,strcat('SegmentedLungs\',header.PatientID,'\Slice No ' ,num2str(i),'.png'));
%          %   dicomdisp(dicomread(strcat('SegmentedLungs\',header.PatientID,'\Slice No ' ,num2str(i),'.dcm')));     
%         %  
%        %Im=histeq(Im);
      
   %    L=EMSeg(Im,40);
  [threshold I1]=maxentropie(Im);
 
 

u = demo_acwe(Im,5);
       seg = u<=0;
  
  erodedBW=seg;

%      for ukp=0:180
%         se=strel('line',2,ukp);
%    erodedBW = imerode(erodedBW,se);
%    end
  
  
 
 
%  se=strel('line',4,90);
%  
%  erodedBW = imerode(erodedBW,se);
%  
 %view 3D
 close all
  
       [B,L] = bwboundaries(medfilt2(erodedBW)); 
       
       
                                           
                                             
                                            
                                            
                                            
                                            clear Lung3D;
  figure, imshow(image,[]); title('Lungs'); hold on     
   
% subplot(2,2,2); imshow(seg);title('Segmentation');
       measurements= regionprops(L, 'all');
                   for a = 1 : length(measurements)
                       
noduleCandidateMask=zeros(r,c);
            nodule3D=zeros(30,30,slices);
            feat=[];
            glcms=[];
                       if(measurements(a).Area<200)
                       
                                 thisBB = measurements(a).BoundingBox;
                          
                                       temp2=image;
                                       temp=image;
                                       temp(:,:)=0;
                                       x=measurements(a).PixelIdxList;
                                       temp(x)=1;
                                       noduleCandidateMask=imfill(temp,'holes');   
                                      
  
                                      nodule=imcrop(image,[thisBB(1),thisBB(2),thisBB(3),thisBB(4)]);
                                       nodule=imresize(nodule,[30 30]);
                                       noduleMask=imcrop(noduleCandidateMask,[thisBB(1),thisBB(2),thisBB(3),thisBB(4)]);
                                       noduleMask=im2bw(imresize(noduleMask,[30 30]));
                                       %3D
                                       volume=dicom23D(FM);
                                                                               
                                       %featureAndPostion= featureExtractionCandidate(volume,[0.7;0.7;1.25],noduleCandidateMask);
                                       [ROIonly,levels,ROIbox,maskBox] = prepareVolume(volume,noduleCandidateMask,'Other',0.7,0.7,1,1.25,'Matrix','Uniform',32);
                                       [GLCM] = getGLCM(ROIonly,levels);
                                       [textures] = getGLCMtextures(GLCM);
                                       %f=featureAndPostion.feature;
                                       clear featureAndPostion;
                                        clear volume;
                                        clear glcms;
                                        clear nodule3D;
                                        HaralickF=[textures.Energy textures.Contrast textures.Entropy textures.Homogeneity textures.Correlation textures.SumAverage textures.Variance textures.Dissimilarity textures.AutoCorrelation];
                                        clear glcms;
                                        clear BW;
                                        clear nodule3D;
                                        clear haralick3D;
                                       I2=imrotate(nodule,90);
                                       mapping=getmaplbphf(8);
                                       h=lbp(nodule,1,8,mapping,'h');
                                       h=h/sum(h);
                                       histograms(1,:)=h;
                                       h=lbp(I2,1,8,mapping,'h');
                                       h=h/sum(h);
                                       histograms(2,:)=h;
                                       lbp_hf_features=constructhf(histograms,mapping);
                                       hog= hog_feature_vector(nodule);
                                       Circularities = (((4 * pi* measurements(a).Perimeter)/measurements(a).Perimeter)^ 2 );
                                        Ellipse=measurements(a).MinorAxisLength/measurements(a).MajorAxisLength;
                                        rectangularity=measurements(a).Extent;
                                        Differenceofcenter=measurements(a).Extrema;
                                        Centre=measurements(a).Centroid;
                                        for m=1:length(Differenceofcenter)
                                            
                                            xaxis=Differenceofcenter(m,1)-Centre(1,1);
                                            yaxis=Differenceofcenter(m,2)-Centre(1,2);
                                        end
                                        variance=var([xaxis yaxis]); 
                                        
                                     
                                        
                                          X=measurements(a).ConvexHull;  
                                        f= double([Circularities  Ellipse rectangularity measurements(a).EulerNumber variance measurements(a).Orientation measurements(a).Eccentricity mean(mean(image(int16(X))))/(max(max(image)))]);
                                          feat=[f HaralickF lbp_hf_features(1,:) hog];
                                       X=isnan(feat);
                                         feat(X)=0;
                                    predChar1 = B_Geo_Har.predict(feat);   
                                   Group = str2double(predChar1);
                                       
                                      
                                    if(Group==1)
                                     
                                     thisBB = measurements(a).BoundingBox;
                                         rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)],...
                                                     'EdgeColor','g','LineWidth',2 );
                                   
                                    end
        %                        1
                            % end

                       end
                   end
                   clc;
         %  close all;        
end
 

           

   
   
   
  
   
 