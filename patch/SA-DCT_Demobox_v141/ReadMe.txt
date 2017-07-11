Pointwise Shape-Adaptive DCT Demobox  -  Public release v1.41 (3 July 2012) 
------------------------------------------------------------------------------------------------------------
(c) 2005-2012
Alessandro Foi
Kostadin Dabov
Tampere University of Technology, Tampere, FINLAND

http://www.cs.tut.fi/~foi/SA-DCT
------------------------------------------------------------------------------------------------------------

The package comprises five demos

  demo_SADCT_denoising.m          (denoising of grayscale images),
  demo_SADCT_color_denoising.m    (denoising of color images),
  demo_SADCT_deblocking.m         (deblocking and deringing of grayscale and color B-DCT compressed images),
  demo_SADCT_deblurring.m         (deblurring of grayscale images),
  demo_SADCT_inverse_halftoning.m (reconstruction of grayscale images from binary halftones)
  
a function

  SADCT_denoising.m               (denoising of grayscale images; similar to demo_SADCT_denoising.m), 

and a function which shows the adaptive-shape neighborhoods which are used in the demos as supports of the
Shape-Adaptive DCT transform

  function_shape_explorer.m      (visualization of the adaptive-shape neighborhoods).


The main demo files are open-source, and may be modified and tuned to be exploited with other data.

Each demo calls a pair of functions in order to perform hard-thresholding and Wiener filtering in SA-DCT domain.

In particular, the denoising and deblocking demos use

  function_SADCT_thresholding_fast.dll     (hard-thresholding in SA-DCT domain),
  function_SADCT_wiener_fast.dll           (Wiener-filtering in SA-DCT domain),

while the deblurring demo uses the pair
 
  function_SADCT_RI_thresholding.dll      (hard-thresholding in SA-DCT domain for regularized inverse),
  function_SADCT_RW_wiener.dll            (Wiener-filtering in SA-DCT domain for regularized-Wiener inverse).

[Corresponding mexw64, mexa64, mexmaci, mexmaci64, and mexglx files are provided for various platforms].
  
Two more functions are called by the color denoising demo in order to perform color-space transformations:

function_rgb2LumChrom.m    (inverse transformation),
function_LumChrom2rgb.m    (forward transformation).


Starting from v1.10, a new function function_AnisLPAICI8.p is introduced to perform the Anisotropic LPA-ICI
adaptive-scale selection. It provides four to five times faster processing for the LPA-ICI part and requires
much less memory than the previous general code based on conv2 and function_ICI code.



Requirements:

This demobox is developed for Matlab version 7 and the dll's are compiled for x86 (Pentium3, Athlon, or better)+Windows.

This demobox requires a few functions from the Anisotropic Nonparametric Image Restoration Demobox
of LASIP Project (available for download at http://www.cs.tut.fi/~lasip/). It is suggested that the
Pointwise Shape-Adaptive DCT Demobox and the Anisotropic Nonparametric Image Restoration Demobox are
installed in the same folder.

The deblocking demo demo_SADCT_deblocking.m uses the function jpeg_read.dll from the Matlab JPEG Toolbox by Phil Sallee.
A complete release of this toolbox is included in our distribution in the file JPEG_Toolbox_MATLAB_1.4.zip.


What's new in this release:
changes from version 1.41 to 1.40:
 + added:    SA_DCT_denoising.m           implements the denoising as a function
 + added:    mex-files for Windows 64 bit, Linux 32 bit, and Linux 64 bit, MacOSX 32 bit, MacOSX 64 bit.
changes from version 1.40 to 1.31:
 + added:    demo_SADCT_inverse_halftoning.m    (demo for reconstruction of grayscale images from binary halftones)
 * updated:  function_rgb2LumChrom.m      (cleaner code)
 * updated:  function_LumChrom2rgb.m      (cleaner code)
 * updated:  demo_SADCT_deblurring.m      (cleaner code)
 * updated:  demo_SADCT_deblocking.m      (added option to disable coefficient alignment)
 * updated:  demo_SADCT_denoising.m       (added option to disable coefficient alignment)
 * updated:  demo_SADCT_color_denoising.m (added option to disable coefficient alignment)
 * updated:  function_SADCT_RW_wiener.dll (minor fix to regularization term)
changes from version 1.31 to 1.30:
 ! bugfix:   demo_SADCT_deblurring.m      (fixed: calculation of RW variance)
changes from version 1.30 to 1.20:
 * updated:  demo_SADCT_deblocking.m      (added option to disable coefficient alignment)
 * updated:  demo_SADCT_denoising.m       (added option to disable coefficient alignment)
 * updated:  demo_SADCT_color_denoising.m (added option to disable coefficient alignment)
 * updated:  function_SADCT_thresholding_fast.dll    (now works with Pentium3/Athlon CPU, coeff. alignment can be disabled)
 * updated:  function_SADCT_wiener_fast.dll     (now works with Pentium3/Athlon CPU, coeff. alignment can be disabled)
 * updated:  function_SADCT_RI_thresholding.dll (now works with Pentium3/Athlon CPU)
 * updated:  function_SADCT_RW_wiener.dll       (now works with Pentium3/Athlon CPU)
changes from version 1.10 to 1.20:
 + added:    demo_SADCT_deblocking.m      (demo for deblocking and deringing of grayscale and color images)
 * updated:  demo_SADCT_denoising.m       (added blocksize parameter to constrain the maximum size of the filtering block)
 * updated:  demo_SADCT_color_denoising.m (added blocksize parameter to constrain the maximum size of the filtering block)
 * updated:  function_shape_explorer.m    (added blocksize parameter to constrain the maximum size of the filtering block)
changes from version 1.04 to 1.10:
 + added:    function_AnisLPAICI8.p       (enables much faster calculation of h_opt_Q, lower memory requirements)
 * updated:  demo_SADCT_denoising.m       (uses function_AnisLPAICI8, added speed-up factor parameter)
 * updated:  demo_SADCT_color_denoising.m (uses function_AnisLPAICI8, added speed-up factor parameter)
changes from version 1.03 to 1.04:
 ! bugfix:   demo_SADCT_deblurring.m      (fixed: boundary handling for calculation of variance of LPA estimates)
changes from version 1.02 to 1.03:
 * updated:  demo_SADCT_denoising.m       (improved comments, cleaner code, new image cropping, and more visualization options)
 * updated:  demo_SADCT_color_denoising.m (improved comments, cleaner code, new image cropping, and more visualization options)
 * updated:  function_rgb2LumChrom.m      (improved comments and cleaner code)
 * updated:  function_LumChrom2rgb.m      (improved comments and cleaner code)
 * updated:  function_shape_explorer.m    (improved comments)
 * updated:  function_SADCT_thresholding_fast.dll   (faster processing of constant blocks)
changes from version 1.01 to 1.02:  
 ! bugfix:   function_SADCT_thresholding_fast.dll   (fixed: number of harmonics for square blocks)
changes from version 1.00 to 1.01:
 ! bugfix:   function_SADCT_thresholding_fast.dll   (fixed: preservation of mean for square blocks)
 + added:    function_shape_explorer.m    (shows adaptive-shape supports of SA-DCT)


Future additions:

Updates, examples, publications, presentations, etc. can be found at  http://www.cs.tut.fi/~foi/SA-DCT


Disclaimer:

Any unauthorized use of these routines for industrial or profit-oriented activities is expressively prohibited.
By downloading and/or using any of these files, you implicitly agree to all the terms of the TUT limited license
(included in the file Legal_Notice.txt).


Questions/feedback:

If you have any comment, suggestion, or question, please do 
contact   Alessandro Foi  at  firstname.lastname@tut.fi


