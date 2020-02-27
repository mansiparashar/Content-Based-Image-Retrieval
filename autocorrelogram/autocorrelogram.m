clc
D = 'C:\Users\student\Desktop\images';
S = dir(fullfile(D,'*.jpg'));
distance_vector=[1 3];
q=imread('C:\Users\student\Desktop\images\18.jpg');
correlogram_vectorq=color_auto_correlogram(q,distance_vector);
info_table = cell2table(cell(0,2),'VariableNames',{'Imagename' 'd' });
for k = 1:numel(S)
   clear F;
   clear I;
   F= fullfile(D,S(k).name);
   I= imread(F);
   correlogram_vectorq1=color_auto_correlogram(I,distance_vector);
   d=abs(correlogram_vectorq-correlogram_vectorq1);
   
   new_row={S(k).name,sum(d)};
   info_table=[info_table;new_row];
 end
info_table=sortrows(info_table,'d');
writetable(info_table,'mansi.csv');
subplot(3,3,2);
imshow(q);
title('Query Image');
file_names=info_table(:,'Imagename').Imagename;
for i=1:6
    F=fullfile(D,char(file_names(i)));
    I=imread(F);
    subplot(3,3,i+3);
    imshow(I);
end
function [positive_count,total_count]=get_n(n,x,y,color,img_no_dither,X,Y)
% This function is useful to get the validity map of the neighborhood case.
% It can handle any number of neighborhood distances.
% Input
% n=The order of the neighborhood
% x & y= x y co-ordinates of the given pixel
% color= particular quantized color
% img_no_dither= The color quantized image matrix
% X & Y= The original dimensions of the input image
% Output
% positive_count= The number of occurences which have the same color
% total_count= The total number of valid cases for this particular instant
    valid_vector8n=zeros(1,8*n); % This is because of the propoerty of inf-norm. Each distance has 8 times the order
    positive_count=0;   total_count=0;
    nbrs_x=zeros(1,8*n);    nbrs_y=zeros(1,8*n);
    % The counting of the pixels is done in the following manner: From the
    % given pixel, go left-->up-->right-->down-->left-->up
    % Y co-ordinates of nbrs
    nbrs_y(1)=y;
    d=1;
    for k=2:1+n
       nbrs_y(k)=y-d;
       d=d+1;
    end
    
    nbrs_y(1+n:1:3*n+1)=y-n;
    
    d=0;
    for k=3*n+1:5*n+1
       nbrs_y(k)=y-n+d;
       d=d+1;
    end
    
    nbrs_y(5*n+1:1:7*n+1)=y+n;
    
    d=0;
    for k=7*n+1:1:7*n+1+(n-1)
       nbrs_y(k)=y+n-d;
       d=d+1;
    end
    
    % X co-ordinates of nbrs
    nbrs_x(1)=x-n;
    
    nbrs_x(2:1:1+n)=x-n;
    
    d=0;
    for k=1+n:1:3*n+1
        nbrs_x(k)=x-n+d;
        d=d+1;
    end
    
    nbrs_x(3*n+1:5*n+1)=x+n;
    
    d=0;
    for k=5*n+1:7*n+1
        nbrs_x(k)=x+n-d;
        d=d+1;
    end
        
    nbrs_x(7*n+1:7*n+1+(n-1))=x-n;
    
    % Assigning the validity of the neighborhood
    for i=1:8*n
        
        if nbrs_x(i)>0 && nbrs_x(i)<=X && nbrs_y(i)>0 && nbrs_y(i)<=Y
            valid_vector8n(i)=1;
        
        else
            valid_vector8n(i)=0;
    
        end
    
    end
    
    
    % Couting the number of common colors in the valid areas of the
    % neighborhood.
    for j=1:8*n
       if valid_vector8n(j)==1
          data= img_no_dither(nbrs_y(j),nbrs_x(j));
          if (data==color)
              positive_count=positive_count+1;
          end
          total_count=total_count+1;
       end
    end
end
function correlogram_vector=color_auto_correlogram(I,distance_vector)
% This function creates the auto-correlogram vector for an input image of
% any size. The different distances which is assumed apriori can be user-defined in a vector.
% It implements the algorithm as defined in Huang et al. paper 'Image Indexing using color
% autocorelogram'
% Input:
% I=The uint8 matrix representing the color image
% distance_vector= The vector representating the different distances in
% which the color distribution is calculated.
% Output:
% correlogram_vector=This is a straight vector representating the
% probabilities of occurence of 64 quantized colors. Its total dimension is
% 64n X 1; where n is the number of different inf-norm distances
% Usage: (To create the auto-correlogram vector for user-defined distances)
% I=imread('peppers.png'); distance_vector=[1 3];
% correlogram_vector=color_auto_correlogram(I,distance_vector);
% Contact Author:
% Soumyabrata Dev
% E-mail: soumyabr001@e.ntu.edu.sg
% http://www3.ntu.edu.sg/home2012/soumyabr001/
correlogram_vector=[];
[Y,X]=size(rgb2gray(I));
% quantize image into 64 colors = 4x4x4, in RGB space
[img_no_dither, ~] = rgb2ind(I, 64, 'nodither');
% figure, imshow(img_no_dither, map);
%rgb = ind2rgb(img_no_dither, map); % rgb = double(rgb)
[~,d]=size(distance_vector);
count_matrix=zeros(64,d);   total_matrix=zeros(64,d);
prob_dist=cell(1,d);
for serial_no=1:1:d
    for x=1:X
        for y=1:Y
            color=img_no_dither(y,x);
       
            % At the given distance 
            [positive_count,total_count]=get_n(distance_vector(serial_no),x,y,color,img_no_dither,X,Y);
            count_matrix(color+1,serial_no)=count_matrix(color+1,serial_no)+positive_count;
            total_matrix(color+1,serial_no)=total_matrix(color+1,serial_no)+total_count;       
        end
    end
    prob_dist{serial_no}=count_matrix(:,serial_no)./(1+total_matrix(:,serial_no));
end
for serial_no=1:d
    correlogram_vector=cat(1,correlogram_vector,prob_dist{serial_no});
end
end