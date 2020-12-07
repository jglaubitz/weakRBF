function p = PolynomialBasis(d)
    % ---------------------------------------------------------------
    % Function returning the handle of the functions evaluating the 
    % P1 2D polynomial basis over a given set of points.
    %
    % Input:  d,  integer,         the desired P polynomial order 
    %
    % Output: p, function hangle,  functions evaluating the P1 2D polynomial 
    %                              basis over a given set of points. When applied to a (NbPoints x 2)
    %                              coordinates vector, returns NbasisFunc x NbPoints
    % ---------------------------------------------------------------
    d
    % Switch on the wished polynomial order and construct hardly the basis
    if d == 0
        p = @(x) x(:,1)'.^0 + x(:,2)'.^0; 
    
    elseif d == 1
        p = @(x) [x(:,1)'.^0 ; x(:,1)'; x(:,2)']; 
    
    elseif d == 2
        p = @(x) [x(:,1)'.^0 ; x(:,1)'; x(:,2)'; x(:,1)'.*x(:,2)'; x(:,1)'.^2; x(:,2)'.^2]; 
    
    elseif d == 3
        p = @(x) [x(:,1)'.^0 ; x(:,1)'; x(:,2)'; x(:,1)'.*x(:,2)'; x(:,1)'.^2; x(:,2)'.^2; x(:,1)'.^3; x(:,2)'.^3]; 
    
    else 
        error('Wished order for the polynomial basis function is not yet implemented. Abortion.');
    end
        
end
