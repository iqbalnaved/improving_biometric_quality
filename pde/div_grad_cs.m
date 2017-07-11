function out = div_grad_cs(u, h, r)

[cpc,cmc,ccp,ccm,cpcp,cmcp,cpcm,cmcm]=eightshifts(u);
%[cpc,cmc,ccp,ccm]=fourshifts(u);

v1S = (cpc-u)./sqrt((cpc-u).^2+((ccp+cpcp-ccm-cpcm).^2)/16 + r);%/h/h;
v1N = (u-cmc)./sqrt((u-cmc).^2+((ccp+cmcp-ccm-cmcm).^2)/16 + r);%/h/h;

v1E = (ccp-u)./sqrt((ccp-u).^2+((cpc+cpcp-cmc-cmcp).^2)/16 + r);%/h/h;
v1W = (u-ccm)./sqrt((u-ccm).^2+((cpc+cpcm-cmc-cmcm).^2)/16 + r);%/h/h;

out = (v1S-v1N)/h + (v1E-v1W)/h;