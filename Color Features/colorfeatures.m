D=’C:\Users\student\Desktop\images’;
imgFolder = fullfile(D);
queryImage = imread('C:\Users\student\Desktop\images\3.jpg');
R = queryImage(:,:,1);
R= double(R(:));
G = queryImage(:,:,2);
G=double(G(:));
B = queryImage(:,:,3);
B = double(B(:));
QueryVals = [mean(R) std(R) skewness(R) kurtosis(R) mean(G) std(G) skewness(G) kurtosis(G) mean(B) std(B) skewness(B) kurtosis(B)];
clear R G B;
filePattern = fullfile(imgFolder, '*.jpg');
jpegFiles = dir(filePattern);
info_table = cell2table(cell(0,14),'VariableNames',{'fname','RedMean','RedStd','RedSkew','RedKurtosis','GreenMean','GreenStd','GreenSkew','GreenKurtosis','BlueMean','BlueStd','BlueSkew','BlueKurtosis','eucledian'});
for k = 1:length(jpegFiles)
  baseFileName = jpegFiles(k).name;
  fullFileName = fullfile(imgFolder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  img = imread(fullFileName);
  R = img(:,:,1);
  R= double(R(:));
  G = img(:,:,2);
  G=double(G(:));
  B = img(:,:,3);
  B = double(B(:));
  ImgVals = [mean(R) std(R) skewness(R) kurtosis(R) mean(G) std(G) skewness(G) kurtosis(G) mean(B) std(B) skewness(B) kurtosis(B)];
  eucledian = sqrt(sum((QueryVals - ImgVals) .^ 2));
  new_row = {fullFileName,ImgVals(1),ImgVals(2),ImgVals(3),ImgVals(4),ImgVals(5),ImgVals(6),ImgVals(7),ImgVals(8),ImgVals(9),ImgVals(10),ImgVals(11),ImgVals(12),dist};
  info_table=[info_table;new_row];
  clear R G B dist ImgVals;
end

info_table = sortrows(info_table,'eucledian');
writetable(info_table,'18bce1048.csv');

subplot(332)
imshow(queryImage); title('Query Image'); 
for i=1:6
    im1=imread(char(info_table(i,1).name));
    [as,sd,fg] = fileparts(char(info_table(i,1).name));
    subplot(3,3,i+3)
    imshow(im1);
    title(sd);
end  