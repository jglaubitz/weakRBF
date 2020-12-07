% Author: Jan Glaubitz 
% Date: 06.11.2018

function DM = DistanceMatrix( X, evalpoints )

    [s, N] = size(X);
    [s, K] = size(evalpoints);
    DM = zeros(K,N);

    for n = 1:K 
        DM(n,:) = abs( evalpoints(n) - X(:) );
    end

end