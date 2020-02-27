clc
clear all
D = 'C:\Users\student\Desktop\images';
S = dir(fullfile(D,'*.jpg'));
I=imread('C:\Users\student\Desktop\images\18.jpg');
info_table = cell2table(cell(0,3),'VariableNames',{'Imagename' 'count' 'ed' });
%partitioning and getting icon
B = imresize(I,[512,512]);
r=B(:,:,1);
g=B(:,:,2);
b=B(:,:,3);
r1 = blockproc(r,[64 64],@(x)round(mean(x.data(:))));
r2 = blockproc(g,[64 64],@(x)round(mean(x.data(:))));
r3 = blockproc(b,[64 64],@(x)round(mean(x.data(:))));
cat=cat(3,r1,r2,r3);
icon=rgb2ycbcr(cat);
y=icon(:,:,1);
cb=icon(:,:,2);
cr=icon(:,:,3);
y1 = dct2(y);
cb1 = dct2(cb);
cr1 = dct2(cr);
a=zigzag(y1);
b=zigzag(cb1);
c=zigzag(cr1);
P=[a,b,c];
for k = 1:numel(S)
   clear F;
   clear I;
   F= fullfile(D,S(k).name);
   I1= imread(F);
   B1 = imresize(I1,[512,512]);
   r1=B1(:,:,1);
   g2=B1(:,:,2);
   b3=B1(:,:,3);
   r11 = blockproc(r1,[64 64],@(x)round(mean(x.data(:))))
   r22 = blockproc(g2,[64 64],@(x)round(mean(x.data(:))))
   r33 = blockproc(b3,[64 64],@(x)round(mean(x.data(:))))
   cat1=cat(3,r11,r22,r33); 
   icon=rgb2ycbcr(cat1);
   y1=icon(:,:,1);
   cb1=icon(:,:,2);
   cr1=icon(:,:,3);
   y11 = dct2(y1);
   cb11 = dct2(cb1);
   cr11 = dct2(cr1);
   a1=zigzag(y11);
   b1=zigzag(cb11);
   c1=zigzag(cr11);
   P1=[a1,b1,c1];
   ed=sqrt((P-P1).^2);
 
   new_row={S(k).name,count1,sum(ed)};
   info_table=[info_table;new_row];
end
info_table=sortrows(info_table,'ed');
writetable(info_table,'colorlayout.csv');
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
 
 
 
function output = zigzag(in)
% initializing the variables
%----------------------------------
h = 1;
v = 1;
vmin = 1;
hmin = 1;
vmax = size(in, 1);
hmax = size(in, 2);
i = 1;
output = zeros(1, vmax * hmax);
%----------------------------------
while ((v <= vmax) & (h <= hmax))
    
    if (mod(h + v, 2) == 0)                 % going up
        if (v == vmin)       
            output(i) = in(v, h);        % if we got to the first line
            if (h == hmax)
          v = v + 1;
        else
              h = h + 1;
            end;
            i = i + 1;
        elseif ((h == hmax) & (v < vmax))   % if we got to the last column
            output(i) = in(v, h);
            v = v + 1;
            i = i + 1;
        elseif ((v > vmin) & (h < hmax))    % all other cases
            output(i) = in(v, h);
            v = v - 1;
            h = h + 1;
            i = i + 1;
     end;
        
    else                                    % going down
       if ((v == vmax) & (h <= hmax))       % if we got to the last line
            output(i) = in(v, h);
            h = h + 1;
            i = i + 1;
        
       elseif (h == hmin)                   % if we got to the first column
            output(i) = in(v, h);
            if (v == vmax)
          h = h + 1;
        else
              v = v + 1;
            end;
            i = i + 1;
       elseif ((v < vmax) & (h > hmin))     % all other cases
            output(i) = in(v, h);
            v = v + 1;
            h = h - 1;
            i = i + 1;
        end;
    end;
    if ((v == vmax) & (h == hmax))          % bottom right element
        output(i) = in(v, h);
        break
    end;
end;
end
