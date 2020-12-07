function [basis, dx_basis, dy_basis] = Solve_EvaluateBasis(rbf, ep, X, d, V)
    % ----------------------------------------------------------------------------------
    % Returning Radial Basis Functions evaluated at a list of points
    %
    % Input:
    %        rbf,      function handle,               the handle of the selected rbf function
    %        ep,       float,                         the epsilon parameter to pass to the rbf function
    %        X,        (NbPoints x 2) float,          the list of xy coordinates of points on which to evaluate the rbf
    %        d,        integer,                       the degree of the considered polynomials when building the basis elements (19)
    %        V,        (NbPoints x NbPoints) float,   the Vandermonde matrix associated to the chosen rbf function and reference nodal points X
    %
    % Output:
    %       basis,     function handle,               function handle that when applied to a list of points evaluates the list of rbf at each point (NbBasis x NbPoints)
    %       dx_basis,  function handle,               function handle that when applied to a list of points evaluates the list of rbf differences in x at each point (NbBasis x NbPoints)
    %       dy_basis,  function handle,               function handle that when applied to a list of points evaluates the list of rbf differences in y at each point (NbBasis x NbPoints)
    %
    % ----------------------------------------------------------------------------------

    %% If no polynomials mollification is included
    if d < 0 
        
        %% Compute directly the RBF handle according to the given set of points X 
        % initial basis and derivatives 
        basis_init = @(x) rbf(ep,sqrt((x(:,1)-X(:,1)').^2+(x(:,2)-X(:,2)').^2)');
        syms x y
        dx_basis_init = matlabFunction(diff(basis_init([x,y]), x)); 
        dy_basis_init = matlabFunction(diff(basis_init([x,y]), y)); 
        
        % Compute the coefficients alpha associated to the rbf and beta associated to the polynomial
        N    = size(X,1);
        y = eye(N); 
        alpha = V\y; 
        
        % Compute nodal basis 
        basis = @(x) (alpha')*basis_init(x);
        dx_basis = @(x,y) (alpha')*dx_basis_init(x,y);
        dy_basis = @(x,y) (alpha')*dy_basis_init(x,y);
        
        
    %% If polynomials mollification is included, compute the polynomials and the coefficients of the rbf
    else
        
        %% Retrieve the handle of the P1 polynomial basis functions according to the desired degree
        
        % initial basis and derivatives 
        basis_init = @(x) rbf(ep,sqrt((x(:,1)-X(:,1)').^2+(x(:,2)-X(:,2)').^2)');
        syms x y
        dx_basis_init = matlabFunction(diff(basis_init([x,y]), x)); 
        dy_basis_init = matlabFunction(diff(basis_init([x,y]), y));          
        
        % Construct the polynomial matrix P over the given list of points 
        Poly = PolynomialBasis(d); % Polynomials
        P = Poly(X);
        syms x y 
        dx_Poly = matlabFunction(diff(Poly([x,y]), x)); 
        dy_Poly = matlabFunction(diff(Poly([x,y]), y));
        
        % Construct the matrix A = ( V P ; P^T 0 )
        % Retrievinf the number of basis functions (P(2d) case) and the number of nodes
        NbFs = (d+2)*(d+1)/2;
        N    = size(X,1);
        A = [V,  P';
             P,  zeros(NbFs, NbFs)];
        
        % Compute the coefficients alpha associated to the rbf and beta associated to the polynomial
        y = zeros(N+NbFs,N); y(1:N,1:N) = eye(N); 
        gamma = A\y; 
        alpha = gamma(1:N,:);
        beta  = gamma(N+1:N+NbFs,:);
       
        % Compute nodal basis 
        basis = @(x) (alpha')*basis_init(x) + (beta')*Poly(x);
        if d==0 
            dx_basis = @(x,y) (alpha')*dx_basis_init(x,y);
            dy_basis = @(x,y) (alpha')*dy_basis_init(x,y);
        else
            dx_basis = @(x,y) (alpha')*dx_basis_init(x,y) + (beta')*dx_Poly(x,y);
            dy_basis = @(x,y) (alpha')*dy_basis_init(x,y) + (beta')*dy_Poly(x,y);
        end
        
    end
end
