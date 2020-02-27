clc
clear all
D = 'C:\Users\91984\Desktop\lbp_images';
S = dir(fullfile(D,'*.jpg'));
q=imread('C:\Users\91984\Desktop\lbp_images\4.jpg');
im1=rgb2gray(q);
aq=extractLBPFeatures(im1,'Upright',false);
info_table = cell2table(cell(0,3),'VariableNames',{'Imagename' 'lbp_feature' 'ed' });
for k = 1:numel(S)
   clear F;
   clear I;
   F= fullfile(D,S(k).name);
   I= imread(F);
   testi=rgb2gray(I);
   a=extractLBPFeatures(testi,'Upright',false);
   ed=sqrt((aq-a).^2);
   
   new_row={S(k).name,a,sum(ed)};
   info_table=[info_table;new_row];
 end
info_table=sortrows(info_table,'ed');
writetable(info_table,'lbp.csv');
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