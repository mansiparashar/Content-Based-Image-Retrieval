D = 'C:\Users\student\Desktop\images;
S = dir(fullfile(D,'*.jpg'));
q=imread('C:\Users\student\Desktop\images\18.jpg');
im1=rgb2gray(q);
[count,bin]=imhist(im1);
info_table = cell2table(cell(0,3),'VariableNames',{'Imagename' 'count' 'ed' });
for k = 1:numel(S)
   clear F;
   clear I;
   F= fullfile(D,S(k).name);
   I= imread(F);
   testi=rgb2gray(I);
   [count1,bin1]=imhist(testi);
   ed=sqrt((count-count1).^2);
   
   new_row={S(k).name,count1,sum(ed)};
   info_table=[info_table;new_row];
 end
info_table=sortrows(info_table,'ed');
writetable(info_table,'mansi.csv');
subplot(3,3,2);
imshow(q);
title('Query Image');
file_names=info_table(:,'Imagename').Imagename;
for i=1:6
    F=fullfile(D,char(file_names(k)));
    I=imread(F);
    subplot(3,3,i+3);
    imshow(I);
end
