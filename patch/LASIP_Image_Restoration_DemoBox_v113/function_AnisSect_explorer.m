function function_AnisSect_explorer(z,h_opt_Q,kernels_higher_order,h1)
% Shows the anisotropic neighborhoods (i.e. union of the supports of the adaptive-scale kernels) which are used for estimation in the Anisotropic LPA-ICI algorithms
%
% Alessandro Foi - Tampere University of Technology -   March 2006
% -----------------------------------------------------------------------------------------------------------------------------------------------
%
%  SYNTAX:
%    function_AnisSect_explorer(z,h_opt_Q,kernels_higher_order,h1);
%
%  INPUTS:
%    z         :  noisy/degraded observation  (used as background image only)
%    h_opt_Q   :  array with the adaptive-scale indexes for the kernels used for the Anisotropic LPA-ICI estimation
%    kernels_higher_order  :  cell array containing all the LPA kernels which are used for estimation in the Anisotropic LPA-ICI algorithms
%    h1        :  set of scales (OPTIONAL)
%
%
%  USAGE:
%    Move the mouse over the image. The corresponding shapes are drawn.
%    Clicking with the mouse cycles through different zoom factors.
%    Closing the figure causes the function to end.
%
%  NOTES:
%    This function is to be used AFTER the
%    denoising/deblurring/etc. routines.
%
%  EXAMPLES:
%    The above syntax example corresponds to the exact syntax to use
%    in order to visualize the neighborhoods used in the demo
%    demo_DenoisingGaussian.m .   A similar syntax applies to all
%    other LASIP filtering routines.
%
%
%      NO FILTERING OR SCALE SELECTION IS PERFORMED WITHIN THIS FUNCTION!
%         This function only visualizes the neighborhoods/sectors
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% there are no parameters to be modified below  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4
    no_h1=1;
else
    no_h1=0;
end

%% image
[size_z_1,size_z_2,size_z_3]=size(z);
if size_z_3>1
    disp(' input image must be grayscale !!! ')
    return
end
z=z-max(min(z(:)),0);  %% normalization
z=z/min(1,max(z(:)));  %% normalization
FrameZ=0.6+0.4*z;  %% whitening for better visibility of shapes


% checks for order availability among sizes and directions
directional_resolution=size(kernels_higher_order,1);


%%% constructs kernels2 cell array where all kernels have the same "support
%%% size" (these are used for the spatial domain plots)
lenh=size(kernels_higher_order,2);
for s1=1:directional_resolution
    kernsizes(s1,:)=size(kernels_higher_order{s1,lenh});
end
h_max=(max(kernsizes(:))+1)/2;
for s1=1:directional_resolution
    for s2=1:lenh
        largest=zeros(max(kernsizes(:)));
        diffe1=(size(largest,1)-size(kernels_higher_order{s1,s2},1))/2;
        diffe2=(size(largest,2)-size(kernels_higher_order{s1,s2},2))/2;
        largest(1+diffe1:size(largest,1)-diffe1,1+diffe2:size(largest,2)-diffe2)=kernels_higher_order{s1,s2};
        kernels_higher_order{s1,s2}=largest(end:-1:1,end:-1:1);
    end
end

for indexh=1:lenh
    for indexq=1:directional_resolution
        Sector{indexq,indexh}=kernels_higher_order{indexq,indexh,1}~=0;
        SectorPerim{indexq,indexh}=bwperim(kernels_higher_order{indexq,indexh,1});
    end
end


h_opt_Q = min(uint8(h_opt_Q),uint8(lenh));
clear kernels_higher_order;

minFigWidth = 300; % don't try to display a figure smaller than this (necessary to have room for text)
minFigHeight = 200;

% What are the screen dimensions
screenSize = get(0, 'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4);
if ((screenWidth <= 1) || (screenHeight <= 1))
    screenWidth = Inf;
    screenHeight = Inf;
end

% choose default zoom_ratio (larger than 1x for high-resolution displays)
default_zoom_ratio=round(1+double((screenWidth>1279)&(screenHeight>959)));
zoom_ratio=default_zoom_ratio;
zoom_ratio_old=-3;
set(0,'userdata',zoom_ratio)   % in this way the zoom_ratio variable is accessible from the subfunction (mouse click)

% size of the frame to be drawn
small_fig_size_1=min(round(0.65*screenHeight/zoom_ratio),size_z_1);
small_fig_size_2=min(round(0.6*screenWidth/zoom_ratio),size_z_2);
small_shift_1=round((size_z_1-small_fig_size_1)/2);
small_shift_2=round((size_z_2-small_fig_size_2)/2);

%% some variables
i1old=-screenWidth;
i2old=-2;
small_shift_2_old=-screenWidth;
small_shift_1_old=0;
figHandle=figure;
colormap(gray)
set(gcf,'windowbuttondownfcn',@wbdfcn);   %% sub-function executed on mouse-click

axHandle=gca;
                 h_max_l=h_max;
while ishandle(axHandle)   %% runs loop as long as the image exists (closing the figure causes the function to end)
    if figHandle==get(0,'CurrentFigure')  %% do not execute if focus is on another window / figure
        zoom_ratio=get(0,'userdata');   %% get possibly updated zoom_ratio (updated from subfunction)

        small_fig_size_1=min(round(0.65*screenHeight/zoom_ratio),size_z_1);  %% set the size of frame
        small_fig_size_2=min(round(0.6*screenWidth/zoom_ratio),size_z_2);

        set(0,'units','pixels');    %%% use pixels as units
        set(gcf,'units','pixels');
        set(axHandle,'units','pixels');

        axPos=get(axHandle,'position');
        figPos=get(gcf,'position');

        %% gets position of the mouse pointer
        mouse_displace=ceil(min([(h_max+2/zoom_ratio),small_fig_size_1/10,small_fig_size_2/10]));  %% diagonal displacement of the mouse with respect to the neighborhood center
        MouseXYcoord=get(0,'pointerlocation');   %% grabs mouse pointer coordinates (on the screen)
        i2_small=round((MouseXYcoord(1)-axPos(1)-figPos(1))/zoom_ratio-mouse_displace);
        i1_small=round((figPos(2)+axPos(4)-MouseXYcoord(2)+axPos(2))/zoom_ratio-mouse_displace);

        %% define shift in image depending on mouse relative position within the frame (useful for panning the frame)
        bound_speed=max(h_max+1,15);   %% this defines maximum panning speed when mouse is on frame border
        small_shift_1=small_shift_1-(max(0,bound_speed-i1_small));
        small_shift_1=small_shift_1+(max(0,bound_speed+i1_small-small_fig_size_1));
        small_shift_1=min(max(0,small_shift_1),size_z_1-small_fig_size_1);
        small_shift_2=small_shift_2-(max(0,bound_speed-i2_small));
        small_shift_2=small_shift_2+(max(0,bound_speed+i2_small-small_fig_size_2));
        small_shift_2=min(max(0,small_shift_2),size_z_2-small_fig_size_2);
        i2=  i2_small+small_shift_2;   %% compensate point coordinates with shifts
        i1=  i1_small+small_shift_1;


        if i1~=i1old||i2~=i2old||zoom_ratio~=zoom_ratio_old  %% new position or new zoom_ratio
            i1old=i1;
            i2old=i2;
            update_bitmap=1;
            if i1<1||i2<1||i1>size_z_1||i2>size_z_2   %% pointer is outside frame borders
                if small_shift_1==small_shift_1_old&&small_shift_2== small_shift_2_old
                    if zoom_ratio==zoom_ratio_old
                        update_bitmap=0;
                    end
                else
                    small_shift_1_old=small_shift_1;
                    small_shift_2_old=small_shift_2;
                end
                outside_border=1;
                pause(0.02)
                correct_back=0;
            else   %% pointer is on the frame
                outside_border=0;
                %% constructs shape and shape subframe
                INPUT_MASK=zeros(h_max+h_max-1);
                INPUT_MASK_OUTLINE=INPUT_MASK;
                for ii=1:directional_resolution
                    h_opt=h_opt_Q(i1,i2,ii);
                    INPUT_MASK=INPUT_MASK|Sector{ii,h_opt};
                    INPUT_MASK_OUTLINE=INPUT_MASK_OUTLINE|SectorPerim{ii,h_opt};
                end
%                 h_max_l=h_max;
                INPUT_MASK=INPUT_MASK(h_max-h_max_l+1:h_max+h_max_l-1,h_max-h_max_l+1:h_max+h_max_l-1);
                INPUT_MASK_OUTLINE=INPUT_MASK_OUTLINE(h_max-h_max_l+1:h_max+h_max_l-1,h_max-h_max_l+1:h_max+h_max_l-1);
                ym=max(1,i1-h_max_l+1);  yM=min(size_z_1,i1+h_max_l-1);  xm=max(1,i2-h_max_l+1);     xM=min(size_z_2,i2+h_max_l-1);   % BOUNDS FOR SLIDING WINDOW
                if yM-ym+4+xM-xm<h_max_l+h_max_l+h_max_l+h_max_l  %% near boundaries
                    INPUT_DATA=zeros(h_max_l+h_max_l-1);
                    INPUT_DATA(h_max_l-i1+ym:h_max_l-i1+yM,h_max_l-i2+xm:h_max_l-i2+xM)=z(ym:yM,xm:xM);    % EXPANDS INPUT DATA TO A SQUARE MASK THAT MAY GO BEYOND IMAGE BOUNDARIES
                else
                    INPUT_DATA=z(ym:yM,xm:xM);
                end
                ShapePatch=~bwperim(INPUT_MASK).*(0.4+~INPUT_MASK_OUTLINE.*(0.2+0.4*INPUT_DATA));  %%% this is the shape (number control the graylever intensities of the different segments in the neighborhood structure)
                ShapePatch(h_max_l,h_max_l)=0;   %% marks central pixel black
                FrameZ(ym:yM,xm:xM)=ShapePatch(h_max_l-i1+ym:h_max_l-i1+yM,h_max_l-i2+xm:h_max_l-i2+xM);   %% put patch on frame
                correct_back=1;  %% turn on flag so that frame will be restored back to original after drawing on the screen
            end

            %% updates bitmap on the screen (draws frame)
            if update_bitmap&&figHandle==get(0,'CurrentFigure')
                if zoom_ratio~=1  %% use nearest-neighbor interp. to produce zoomed image
                    image(FrameZ(floor(small_shift_1+1:1/zoom_ratio:small_shift_1+small_fig_size_1+(zoom_ratio-1)/zoom_ratio),floor(small_shift_2+1:1/zoom_ratio:small_shift_2+small_fig_size_2+(zoom_ratio-1)/zoom_ratio))*64);
                else
                    image(FrameZ(small_shift_1+1:small_shift_1+small_fig_size_1,small_shift_2+1:small_shift_2+small_fig_size_2)*64);
                end
                axis off
            end

            %% title for figure, with coordinates and scales / zoom_ratio
            if figHandle==get(0,'CurrentFigure')
                titleHandle = get(axHandle,'title');
                if zoom_ratio_old~=zoom_ratio
                    set(titleHandle, 'string',['Zoom ',num2str(zoom_ratio),'x']);
                else
                    if  outside_border
                        set(titleHandle, 'string',['(',num2str(i1),',',num2str(i2),')    move pointer on image!']);
                    else
                        if no_h1
                            set(titleHandle, 'string',['\itx\rm=(',num2str(i1),',',num2str(i2),')    h\_opt\_Q(\itx\rm,:)=\{', num2str(double(h_opt_Q(i1,i2,:))),'\}']);
                        else
                            set(titleHandle, 'string',['\itx\rm=(',num2str(i1),',',num2str(i2),')    \{\ith\rm^+(\itx\rm,\it\theta_k\rm)\}_{\itk\rm=1,...,8}=\{', num2str(h1(h_opt_Q(i1,i2,:))),'\}']);
                        end
                    end
                end
            end

            %% updates sizes/borders/etc. of figure
            if zoom_ratio_old~=zoom_ratio&&figHandle==get(0,'CurrentFigure')
                set(axHandle, 'Units', 'pixels');
                axPos = get(axHandle, 'Position');
                set(figHandle, 'Units', 'pixels');
                set(0, 'Units', 'pixels');
                figPos = get(figHandle, 'Position');
                defAxesPos = get(0,'DefaultAxesPosition');
                gutterWidth = round((1 - defAxesPos(3)) * small_fig_size_2*zoom_ratio / defAxesPos(3));
                gutterHeight = round((1 - defAxesPos(4)) * small_fig_size_1*zoom_ratio / defAxesPos(4));
                newFigWidth = small_fig_size_2*zoom_ratio + gutterWidth;
                newFigHeight = small_fig_size_1*zoom_ratio + gutterHeight;
                newFigWidth = max(newFigWidth, minFigWidth);
                newFigHeight = max(newFigHeight, minFigHeight);

                figPos(1) = max(1, figPos(1) - floor((newFigWidth - figPos(3))/2));
                figPos(2) = min(screenHeight-newFigHeight-60,max(1, figPos(2) - floor((newFigHeight - figPos(4))/2)));
                figPos(3) = newFigWidth;
                figPos(4) = newFigHeight;

                gutterWidth = figPos(3) - small_fig_size_2*zoom_ratio;
                gutterHeight = figPos(4) - small_fig_size_1*zoom_ratio;
                gutterLeft = floor(gutterWidth/2);
                gutterBottom = floor(gutterHeight/2);

                axPos(1) = gutterLeft + 1;
                axPos(2) = gutterBottom + 1;
                axPos(3) = max(small_fig_size_2*zoom_ratio,1);
                axPos(4) = max(small_fig_size_1*zoom_ratio,1);

                set(figHandle, 'Position', figPos);  %% updates to screen
                set(axHandle, 'Position', axPos);
            end
            if zoom_ratio_old~=zoom_ratio  % pause when zoom ratio changes
                pause(0.3)
                if zoom_ratio==1  % extra pause when zoom cycles back to 1
                    pause(0.2)
                end
            end
        end

        %% corrects back the image (replace the patch with background image)
        if correct_back
            FrameZ(ym:yM,xm:xM)=0.6+0.4*z(ym:yM,xm:xM);
        end
        zoom_ratio_old=zoom_ratio;
        if update_bitmap==0
            pause(0.03)
        end
    end
    pause(0.03)    %% PAUSE TO ALLOW REFRESH OF FIGURES AND TO SLOW DOWN LOOP CYCLING WHEN NO UPDATE IS NECESSARY
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wbdfcn(varargin)   % change zoom_ratio when mouse button is pressed
zoom_ratio=get(0,'userdata');
zoom_ratio=ceil((zoom_ratio+0.1)^1.1);   %% increase zoom ratios in a sensible manner (1 2 3 4 5 7 9 12 16 22 31 44 65 99 etc.)
screenSize = get(0, 'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4);
if ((screenWidth <= 1) || (screenHeight <= 1))
    screenWidth = Inf;  screenHeight = Inf;
end
if zoom_ratio>min(screenWidth/100,30);  %% cycles back to 1x when zoom exceeds a certain ratio
    zoom_ratio=1;
end
set(0,'userdata',zoom_ratio);
return