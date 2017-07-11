# Improving Biometric Quality of Noisy Images

This project impements a face image quality assessment method proposed by Wong et al. (Wong, 2011). Then applies noises on test images and tries to improve the quality by applying various state-of-the-art denoising methods.

> (Wong,2011) 'Wong, Yongkang, et al. "Patch-based probabilistic image quality assessment for face selection and improved video-based face 
recognition." Computer Vision and Pattern Recognition Workshops (CVPRW), 2011 IEEE Computer Society Conference on. IEEE, 2011.

## Data:

Train data : FERET images
Test data: LFW images

## File description:

train_patch_model.m - trains using the images in train_data folder (already preprocesssed)
test_patch_model.m - show demo of quality assessment using a single image from test_data
patch_model.mat - trained models
get_quality_score.m - 
test_blur_img - Apply noise to a single image
test_blur_img_100.m - Apply to a 100 image
test_gaussnoisy_img.m - Apply gaussian noise to a single image
test_gaussnoisy_img_100.m - Apply to a 100 image

## 3rd party codes (to apply different types of noise/denoising methods):

blurtools, bm3d, nonlocal, patch, pde, wavelet


## For more details:

slides.pptx - Explanation of the project
result.xlsx - Experiment results


