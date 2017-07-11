function [cpc,cmc,ccp,ccm,cpcp,cmcp,cpcm,cmcm]=eightshifts(u)
    [cpc,cmc,ccp,ccm] = fourshifts(u);
    cpcp = shiftImage( u, 1, 1 );
    cmcp = shiftImage( u, -1 , 1);
    cpcm = shiftImage( u, 1, -1 );
    cmcm = shiftImage( u, -1, -1 );
