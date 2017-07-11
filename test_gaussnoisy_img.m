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

filename = 'Jennifer_Aniston_0016.jpg';

% method = 'tv'; % 'pm' 'tv'
% addpath('pde\');

% method = 'weiner'; % 'weiner' 'bls-gsm' 
% addpath('patch\');

% method = 'dct'; % 'dct' 'wav'
% addpath('wavelet\');

% method = 'bfilt'; % 'nlm' 'bfilt'
% addpath('nonlocal\');

% method = 'bm3d';
% addpath('bm3d\');

I = imread(['C:\Users\Iqbal\Documents\face_image_quality\lfw_frontals\0\' filename]);

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
    figure;
    box = bbox(b,:);
    fsize = max(round((box(3)/60)^(0.5)*10),8);

    crgr = imcrop(Igr, bbox(b,:));
    crgr = crgr*255;
    crgr = (crgr-mean(crgr(:)))*std0/std(crgr(:))+m0;
    crgr(crgr<0)=0; crgr(crgr>255)=255;        

%     quality = get_quality_score(crgr, patchRows, patchCols, patchSize, featureMean, featureCovariance);
%     quality = quality/(patchRows*patchCols);
%     quality = min(max(quality,0), 100);
%     imshow(crgr,[]);
%     text(box(1)+box(3)*0.02, box(2)+box(3)*0.99, num2str(round(quality)), 'Color', [0 0 0], 'VerticalAlignment', 'top', ...
%         'BackgroundColor', [1 0 1], 'HorizontalAlignment', 'left', 'FontSize', fsize, 'FontWeight', 'bold');

    k=1;
    for sigma_w = 0:5:20
        crgr_noisy = crgr+randn(size(crgr))*sigma_w;
        noisy_quality = get_quality_score(crgr_noisy,  patchRows, patchCols, patchSize, featureMean, featureCovariance);
        noisy_quality = noisy_quality/(patchRows*patchCols);
        noisy_quality = min(max(noisy_quality,0), 100);
        subplot(2, 5, k); imshow(crgr_noisy,[]);
        text(box(1)+box(3)*0.02, box(2)+box(3)*0.99,sprintf('%d,%.2f', round(noisy_quality), psnr(crgr_noisy/255, crgr/255)), 'Color', [0 0 0], 'VerticalAlignment', 'top', ...
            'BackgroundColor', [1 0 1], 'HorizontalAlignment', 'left', 'FontSize', fsize, 'FontWeight', 'bold');
    
        if strcmp(method,'pm') == 1
            crgr_denoised = pm(crgr_noisy, 50); 
        elseif strcmp(method,'tv') == 1
            crgr_denoised = tv(crgr_noisy, 50); 
        elseif strcmp(method,'weiner') == 1
            crgr_denoised=wiener2(crgr_noisy/255,sigma_w)*255;
        elseif strcmp(method,'bls-gsm') == 1
            % Wavelet-based BLS-GSM denoising
            addpath('patch\BLS-GSM');
            addpath('patch\BLS-GSM\denoising_subprograms');
            addpath('patch\BLS-GSM\Added_PyrTools');
            addpath('patch\BLS-GSM\matlabPyrTools');
            crgr_denoised=BLS_GSM_denoise(crgr_noisy,sigma_w);
        elseif strcmp(method,'sadct') == 1
            % SADCT-based  denoising [Doesn't work due to outdated
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
        subplot(2, 5, k+5); imshow(crgr_denoised,[]);
        text(box(1)+box(3)*0.02, box(2)+box(3)*0.99, sprintf('%d,%.2f', round(denoised_quality), psnr(crgr_denoised/255, crgr/255)), 'Color', [0 0 0], 'VerticalAlignment', 'top', ...
            'BackgroundColor', [1 0 1], 'HorizontalAlignment', 'left', 'FontSize', fsize, 'FontWeight', 'bold');
        k=k+1;
    end
%     uicontrol('Style', 'text', 'String', filename, 'FontSize', fsize, 'FontWeight', 'bold', 'Position',[200 375 220 20]);
end
drawnow;    
