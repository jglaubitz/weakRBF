%% generate_initBasis_1d
% Jan Glaubitz, Oct 30, 2020 
%
% Generates an intial basis of monomials and corresponding
% moments on a smooth subinterval [a,b]
% INPUT 
%  a, b : boundaries of the smooth subinterval 
%  d : maximum degree 
% OUTPUT 
%  y_min : single vector containing the result pf the minmod function 

function [basis, m] = generate_initBasis_1d(a, b, d)

    alpha = (0:d)'; % exponent vector 
    basis = @(x) x.^alpha; % initial basis 
    
    % compute the moments 
    m = zeros(d+1,1);
    for k=0:d 
        m(k+1) = ( b^(k+1) - a^(k+1) )/(k+1); 
    end

end 