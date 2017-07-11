function [cpc,cmc,ccp,ccm]=fourshifts(u)
    cpc = shiftImage( u, 1, 0 );
    cmc = shiftImage( u, -1 , 0);
    ccp = shiftImage( u, 0, 1 );
    ccm = shiftImage( u, 0, -1 );
