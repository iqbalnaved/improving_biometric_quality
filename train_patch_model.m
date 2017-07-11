% Implementation of the paper titled "Wong, Yongkang, et al. "Patch-based 
% probabilistic image quality assessment for face selection and improved
%  video-based face recognition." Computer Vision and Pattern Recognition 
% Workshops (CVPRW), 2011 IEEE Computer Society Conference on. IEEE, 2011.

% training phase
% preprocessed train data: training images are well aligned, centered faces, frontal
% face images under appropriate lighiting with neutral facial expression.
% We have used the FERET dataset.

% training method: the method is initialized by constructing a generic face
% model. One single model is constructed based on the multiple training
% faces. firstly an image patch is moved in single steps accross the image
% along a predefined path. The size of the image is 64x64 pixels, size of
% the patch is 8x8 , with overlap by 7 pixels (in both row and column). As
% a result there are 64-8+1=57 patches in each line of the image and
% 57*57=3249 patches in the entire image.

% For each of the patch of each training images a 2D-DCT is applied and 8x8
% coffecients are determined. Only the frequences of a 2x2 squer at the top
% lefet of the matrix is considered because lower frequency components lie
% towards the upper left of the matrix and they characterise the coarse
% appearance and descriptive properties of a face. Ignoring the DC
% component , the top left 2x2 element of matrix results in a feature
% vector of dimension 3. From all the training images of this patch
% location mean vector and 3x3 covariance matrix of the set of vectors is
% determined. This is done for all the patch locations for all the training
% images. 

% These mean vector and covariance matrix is used at test phase.

trainPath = 'train_data\';
files = dir([trainPath '\*.png']);
trainSize = length(files);

imSize = [64 64];
patchSize = [8 8];

patchRows = imSize(1) - patchSize(1) + 1; % 64-8+1 = 57
patchCols = imSize(2) - patchSize(2) + 1; % 64-8+1 = 57

featureStore = cell(patchRows,patchCols,trainSize);

for k = 1:trainSize
    I = double(imread([trainPath files(k).name]));
    I = imresize(I, [imSize(1) imSize(2)]); 
    
    for i = 1:patchRows
        for j = 1:patchCols
            patch = I(i:i+patchSize(1)-1, j:j+patchSize(2)-1);
            J = dct2(patch);
            feature = [J(1,2); J(2,1); J(3,1)];
            featureStore{i,j,k} = feature;
        end
    end
end

featureMean = cell(patchRows,patchCols);
featureCovariance = cell(patchRows,patchCols); 
for i = 1:patchRows
    for j = 1:patchCols
        patchFeatures = featureStore(i,j,:);
        patchFeatures = cell2mat(patchFeatures);
        featureMean{i,j} = mean(patchFeatures,3);
        patchFeatures = reshape(patchFeatures,[trainSize 3]);
        featureCovariance{i,j} = cov(patchFeatures);
    end
end

save('patch_model.mat', 'featureMean', 'featureCovariance');
