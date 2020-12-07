% Author: Jan Glaubitz 
% Date: 06.11.2018

function d = DifferenceMatrix(x)

    [s, N] = size(x);
    d = zeros(N,N);

    for n = 1:N 
        d(n,:) = x(n) - x(:);
    end

end