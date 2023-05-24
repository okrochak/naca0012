%% Self made code to post-process PIV image pair

% Inputs:
img1name = 'Image_0001_a.tif';
img2name = 'Image_0001_b.tif';

Ninter = 32; % interrogation width in pixels
tol = 1.; % min SNR ratio
img1 = imread(img1name);
img2 = imread(img2name);
img1 = img1(1:1024,1:1344);
img2 = img2(1:1024,1:1344);
DoPlot = 0;

[ht, wd] = size(img1); % image dimensions in px
%% 1. Overlay the image over each other
if DoPlot == 1
figure(1)
imshow(0.5*img1 + 0.5*img2)
end
%% 2. Use cross-correlation

Nwd = wd / Ninter; Nht = ht / Ninter;

% Assemble the interrogation windows matrix:
I1 = zeros(Nht,Nwd,Ninter,Ninter); I2 = I1; dispx = zeros(Nht,Nwd); dispy = dispx;
corrMat = zeros(Nht,Nwd,Ninter*2-1,Ninter*2-1);
SNR = dispx;
for i = 1:Nht % row
    for j = 1:Nwd % column
        mat1 = img1((i-1)*Ninter + 1:i*Ninter,(j-1)*Ninter + 1:j*Ninter); 
        mat2 = img2((i-1)*Ninter + 1:i*Ninter,(j-1)*Ninter + 1:j*Ninter);
        I1(i,j,:,:) = mat1;
        I2(i,j,:,:) = mat2;
        mat1 = mat1 - mean(mat1,'all');
        mat2 = mat2 - mean(mat2,'all');        
        corr = xcorr2(mat1,mat2) ./ (std(double(mat1),0,'all')*std(double(mat2),0,'all'));
        corrMat(i,j,:,:) = corr;
        [peak1, id] = max(corr,[],'all');
        [dx, dy] = ind2sub(size(corr),id);
        corr(id) = 0;
        [peak2, id2] = max(corr,[],'all');
        SNR(i,j) = peak1/peak2;
        dispx(i,j) = dx - Ninter;
        dispy(i,j) = dy - Ninter;
        if SNR(i,j) < tol
            dispx(i,j) = nan;
            dispy(i,j) = nan;
        end
    end
end


% Visualize the interrogation window
if DoPlot == 1
figure(2)
subplot(1,3,1)
imshow(squeeze(uint8(I1(10,10,:,:))));
subplot(1,3,2)
imshow(squeeze(uint8(I2(10,10,:,:))));
subplot(1,3,3)
imshow(squeeze(uint8(I1(10,10,:,:))*0.5 + uint8(I2(10,10,:,:))*0.5));

figure(3)
surf(squeeze(corrMat(15,15,:,:)))

figure(4)
quiver(dispy,dispx)
end
% 
% 
% function C = myCorr(mat1,mat2) %updated correlation function
%     [M, N] = size(mat1); 
%     
% end

figure(10)
imshow(0.5*img1 + 0.5*img2); hold on;
quiver(Ninter/2:Ninter:wd,Ninter/2:Ninter:ht,dispy,dispx);
