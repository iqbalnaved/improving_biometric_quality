function quality = get_quality_score(crgr,  patchRows, patchCols, patchSize, featureMean, featureCovariance)
        crgr = imresize(crgr/255.0, [64 64], 'bilinear');                    
        crgr = log2(crgr+1);

        quality = 0;
        for i = 1:patchRows
            for j = 1:patchCols
                patch = crgr(i:i+patchSize(1)-1, j:j+patchSize(2)-1);
                patch = (patch - mean(patch(:))) / var(patch(:));
                J = dct2(patch);
                feature = [J(1,2); J(2,1); J(3,1)];
                patchMean = featureMean{i,j}; 
                patchCov = featureCovariance{i,j};
                prob = exp(-0.5*(feature-patchMean)' * inv(patchCov) * (feature - patchMean) ) / (2*pi)^(length(feature)/2) * det(patchCov)^0.5;
                quality = quality + log2(prob);                
%                 fprintf('prob=%f, log2(prob)=%f\n', prob, log2(prob));
            end
        end
end