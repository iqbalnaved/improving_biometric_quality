close all;clear;clc;

% test phase
% test data: unline test images, images in test data will have pose
% variations, varying shadow conditions, blurring as well as alignment
% errors to show quality variations from the ideal state.
% Viola-Jones haar-based face detector is used to detect face, then resized
% to 64x64 to match with training image dimension. Image log transformed to
% normalize intensity. Now, the image is segmented into patches in the same
% way as training phase. Each patch is normalized to zero mean and unit
% variance. DCT coeffecients are extracted. Now to determine the similarity
% of each patch of the face image to the corrresponding patch of the
% generic face model a posterior probbility for each patch is determining
% by vector into a Gaussian probability density function based on the
% generic face model.

% The aim is not to recognize an individual face but to determine ow lose
% the patch of the face from the test set is to an average patch at the
% same location represented by the generic facce model. The formua assigns
% a low probability to patches with DCT coefficients far fromt he mean of
% the probabilistic model and a high probability to patches with DCT
% coefficients close tot he mean of the probabilitic model.

% The quality of the image is determined by combining all the posterior
% probabilities of the patches.

minsz = [75 75];
imSize = [64 64];
patchSize = [8 8];

patchRows = imSize(1) - patchSize(1) + 1; % 64-8+1 = 57
patchCols = imSize(2) - patchSize(2) + 1; % 64-8+1 = 57

load('patch_model');

addpath('blurtools\');

method = 'shock'; % 'shock' 'unsharp' 'histeq'

files = dir('test_data\*.jpg');

total_noisy_quality = zeros(1,6);
total_noisy_psnr = zeros(1,6);

total_denoised_quality = zeros(1,6);
total_denoised_psnr = zeros(1,6);

fc = 0;
for f = 1:length(files)
    fprintf('processing %s\n', files(f).name);
    filename = files(f).name;


    I = imread(['test_data\' filename]);

    if max(size(I))>1200 % avoid too large image
        I = imresize(I, 1200/max(size(I)), 'bilinear');
    end
    I = im2double(I);
    try
        Igr = rgb2gray(I);
    catch exception
    end
    fDect = vision.CascadeObjectDetector(); 
    fDect.ScaleFactor = 1.02; fDect.MergeThreshold = 2; fDect.MinSize = minsz; 
    bbox = step(fDect, Igr);

    %% crop and resize face
    m0 = 150; std0 = 50;
    for b = 1:size(bbox,1)
        % display result
        box = bbox(b,:);
        fsize = max(round((box(3)/60)^(0.5)*10),8);

        crgr = imcrop(Igr, bbox(b,:));
        crgr = crgr*255;
        crgr = (crgr-mean(crgr(:)))*std0/std(crgr(:))+m0;
        crgr(crgr<0)=0; crgr(crgr>255)=255;        

        k=1;
        for t = -1:4
            if t < 0
                itr = 0;
            else
                itr = 2^t;
            end
            crgr_noisy = heat_diffusion(crgr,itr); % itr=0,1,2,4,8,16
            noisy_quality = get_quality_score(crgr_noisy,  patchRows, patchCols, patchSize, featureMean, featureCovariance);
            noisy_quality = noisy_quality/(patchRows*patchCols);
            noisy_quality = min(max(noisy_quality,0), 100);
            total_noisy_quality(k) = total_noisy_quality(k) + round(noisy_quality);
            total_noisy_psnr(k) = total_noisy_psnr(k) + psnr(crgr_noisy/255,crgr/255);

            if strcmp(method,'unsharp') == 1
                lambda = 0.1;
                crgr_denoised = unsharp_masking(crgr_noisy, lambda * (t+1)/10); 
            elseif strcmp(method,'shock') == 1
                crgr_denoised = shock(crgr_noisy, 10); 
            elseif strcmp(method,'histeq') == 1
                crgr_denoised = histeq(crgr_noisy/255);
            end
            denoised_quality = get_quality_score(crgr_denoised,  patchRows, patchCols, patchSize, featureMean, featureCovariance);
            denoised_quality = denoised_quality/(patchRows*patchCols);
            denoised_quality = min(max(denoised_quality,0), 100);                    
            total_denoised_quality(k) = total_denoised_quality(k) + round(denoised_quality);
            total_denoised_psnr(k) = total_denoised_psnr(k) + psnr(crgr_denoised/255, crgr/255);
            k=k+1;
            fc = fc + 1;
        end
    end
end

fprintf('itr=0, tnq=%.2f,  tnp=%.2f, tdq=%.2f, tdp=%.2f\n', total_noisy_quality(1)/fc,total_noisy_psnr(1)/fc,total_denoised_quality(1)/fc,total_denoised_psnr(1)/fc);
fprintf('itr=1, tnq=%.2f,  tnp=%.2f, tdq=%.2f, tdp=%.2f\n', total_noisy_quality(2)/fc,total_noisy_psnr(2)/fc,total_denoised_quality(2)/fc,total_denoised_psnr(2)/fc);
fprintf('itr=2, tnq=%.2f,  tnp=%.2f, tdq=%.2f, tdp=%.2f\n', total_noisy_quality(3)/fc,total_noisy_psnr(3)/fc,total_denoised_quality(3)/fc,total_denoised_psnr(3)/fc);
fprintf('itr=4, tnq=%.2f,  tnp=%.2f, tdq=%.2f, tdp=%.2f\n', total_noisy_quality(4)/fc,total_noisy_psnr(4)/fc,total_denoised_quality(4)/fc,total_denoised_psnr(4)/fc);
fprintf('itr=8, tnq=%.2f,  tnp=%.2f, tdq=%.2f, tdp=%.2f\n', total_noisy_quality(5)/fc,total_noisy_psnr(5)/fc,total_denoised_quality(5)/fc,total_denoised_psnr(5)/fc);
fprintf('itr=16, tnq=%.2f,  tnp=%.2f, tdq=%.2f, tdp=%.2f\n', total_noisy_quality(6)/fc,total_noisy_psnr(6)/fc,total_denoised_quality(6)/fc,total_denoised_psnr(6)/fc);

