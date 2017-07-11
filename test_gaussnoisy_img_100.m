close all;clear;clc;

addpath('patch\BLS-GSM');
addpath('patch\BLS-GSM\denoising_subprograms');
addpath('patch\BLS-GSM\Added_PyrTools');
addpath('patch\BLS-GSM\matlabPyrTools');

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

% method = 'tv'; % 'pm' 'tv'
% addpath('pde\');

% method = 'bls-gsm'; % 'weiner' 'bls-gsm' 
% addpath('patch\');

% method = 'dct'; % 'dct' 'wav'
% addpath('wavelet\');
% 
% method = 'bfilt'; % 'nlm' 'bfilt'
% addpath('nonlocal\');
% 
method = 'bm3d';
addpath('bm3d\');


total_noisy_quality = zeros(1,5);
total_psnr_noisy = zeros(1,5);
total_denoised_quality =  zeros(1,5); 
total_psnr_denoised = zeros(1,5);

files = dir('test_data\*.jpg');
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
        box = bbox(b,:);
        fsize = max(round((box(3)/60)^(0.5)*10),8);

        crgr = imcrop(Igr, bbox(b,:));
        crgr = crgr*255;
        crgr = (crgr-mean(crgr(:)))*std0/std(crgr(:))+m0;
        crgr(crgr<0)=0; crgr(crgr>255)=255;        

        k=1;
        for sigma_w = 0:5:20
            crgr_noisy = crgr+randn(size(crgr))*sigma_w;
            noisy_quality = get_quality_score(crgr_noisy,  patchRows, patchCols, patchSize, featureMean, featureCovariance);
            noisy_quality = noisy_quality/(patchRows*patchCols);
            noisy_quality = min(max(noisy_quality,0), 100);
            total_noisy_quality(k) = total_noisy_quality(k) + round(noisy_quality);
            total_psnr_noisy(k) = total_psnr_noisy(k) + psnr(crgr_noisy/255, crgr/255);
                

            if strcmp(method,'pm') == 1
                crgr_denoised = pm(crgr_noisy, 50); 
            elseif strcmp(method,'tv') == 1
                crgr_denoised = tv(crgr_noisy, 50); 
            elseif strcmp(method,'weiner') == 1
                crgr_denoised=wiener2(crgr_noisy/255,sigma_w)*255;
            elseif strcmp(method,'bls-gsm') == 1
                % Wavelet-based BLS-GSM denoising
                crgr_denoised=BLS_GSM_denoise(crgr_noisy,sigma_w);
            elseif strcmp(method,'sadct') == 1
                % !! SADCT-based  denoising [Doesn't work due to outdated
                % function_AnisLPAICI8.p file]
                addpath('patch\SA-DCT_Demobox_v141');
                addpath('patch\LASIP_Image_Restoration_DemoBox_v113');
                crgr_denoised=SADCT_denoising(crgr_noisy/255,sigma_w)*255;
            elseif strcmp(method,'dct') == 1
                th1=100;
                crgr_denoised=dct_denoise(crgr_noisy,th1);
            elseif strcmp(method,'wav') == 1
                th2=100;
                wname='bior4.4';
                crgr_denoised=wavelet_denoise(crgr_noisy,wname,th2);
            elseif strcmp(method, 'nlm') == 1
                t=11;f=7;h=sigma_w;
                crgr_denoised=NLmeansfilter(crgr_noisy,t,f,h);
            elseif strcmp(method, 'bfilt') == 1
                w=2;
                sigma=[2 .2];
                crgr_denoised=bfilter2(crgr_noisy/255,w,sigma);
            elseif strcmp(method, 'bm3d') == 1
                crgr_denoised = BM3D_denoising(crgr_noisy/255,sigma_w)*255;
            end

            denoised_quality = get_quality_score(crgr_denoised,  patchRows, patchCols, patchSize, featureMean, featureCovariance);
            denoised_quality = denoised_quality/(patchRows*patchCols);
            denoised_quality = min(max(denoised_quality,0), 100);        
            total_denoised_quality(k) =  total_denoised_quality(k) + round(denoised_quality); 
            total_psnr_denoised(k) = total_psnr_denoised(k) + psnr(crgr_denoised/255, crgr/255);
            k=k+1;
            fc = fc + 1; % total face count
        end
    end
end

fprintf('sigma_w=0, tnq=%.2f,  tpn=%.2f, tdq=%.2f, tpd=%.2f\n', total_noisy_quality(1)/fc,total_psnr_noisy(1)/fc,total_denoised_quality(1)/fc,total_psnr_denoised(1)/fc);
fprintf('sigma_w=5, tnq=%.2f,  tpn=%.2f, tdq=%.2f, tpd=%.2f\n', total_noisy_quality(2)/fc,total_psnr_noisy(2)/fc,total_denoised_quality(2)/fc,total_psnr_denoised(2)/fc);
fprintf('sigma_w=10, tnq=%.2f, tpn=%.2f, tdq=%.2f, tpd=%.2f\n', total_noisy_quality(3)/fc,total_psnr_noisy(3)/fc,total_denoised_quality(3)/fc,total_psnr_denoised(3)/fc);
fprintf('sigma_w=15, tnq=%.2f, tpn=%.2f, tdq=%.2f, tpd=%.2f\n', total_noisy_quality(4)/fc,total_psnr_noisy(4)/fc,total_denoised_quality(4)/fc,total_psnr_denoised(4)/fc);
fprintf('sigma_w=20, tnq=%.2f, tpn=%.2f, tdq=%.2f, tpd=%.2f\n', total_noisy_quality(5)/fc,total_psnr_noisy(5)/fc,total_denoised_quality(5)/fc,total_psnr_denoised(5)/fc);
