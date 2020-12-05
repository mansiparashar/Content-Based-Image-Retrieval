clc
clear all
imagefiles = dir('C:\Users\91984\Desktop\faces\*.jpg');
D = 'C:\Users\91984\Desktop\faces';
nfiles = length(imagefiles);
for ii=1:nfiles
currentfilename = imagefiles(ii).name;
cimage = double(imread(fullfile(D,currentfilename)));
ccoeff = pca(cimage(:,1:8));
cfea=ccoeff(:);
res = [ii,cfea'];
res=res(:,1:51);
dlmwrite('result2.csv',res,'-append')
end
I = double(imread('subject1.gif'));
X = reshape(I,size(I,1)*size(I,2),3);
coeff = pca(X);
fea=coeff(:);
imagefiles = dir('yalefaces/*.gif');
nfiles = length(imagefiles);
sim = zeros(1,nfiles);
for ii=1:nfiles
 currentfilename = strcat(int2str(ii),'.gif');
 cimage = double(imread(currentfilename));
 cX = reshape(cimage,size(cimage,1)*size(cimage,2),3);
 ccoeff = pca(cX);
 cfea=ccoeff(:);
 sim(ii) = sqrt(sum((fea - cfea) .^ 2));
end
[ASorted, AIdx] = sort(sim);
smallestNElements = ASorted(1:6);
smallestNIdx = AIdx(1:6);
for ii=1:6
 imagename = strcat(int2str( smallestNIdx(ii)),'.gif');
 im = imread(imagename);
 subplot(3,2,ii),imshow(im);
end