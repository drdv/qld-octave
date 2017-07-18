% [x, obj, info, lambda] = qld (x0, H)
% [x, obj, info, lambda] = qld (x0, H, g)
% [x, obj, info, lambda] = qld (x0, H, g, A, b)
% [x, obj, info, lambda] = qld (x0, H, g, A, b, lb, ub)
% [x, obj, info, lambda] = qld (x0, H, g, A, b, lb, ub, Ain, bin)
% [x, obj, info, lambda] = qld (..., options)
%
% Solve the quadratic program
%
%      min 0.5 x'*H*x + x'*g
%       x
%
% subject to
%
%      A*x = b
%      lb <= x <= ub
%      A_in*x <= bin
%
% Input
%
%   x0 -- initial guess (not sure if it makes any difference TODO)
%
%   Variables x0, g, A, b, lb, ub, Ain, bin can be empty matrices.
%
%   options -- a structure containing solver options:
%       1) 'hessian_mode': 1 -- normal Hessian, 0 -- H is upper
%       triangular factor of Cholesky decomposition of the Hessian.
%
%       2) 'tolerance' -- tolerance.
%
%
% Output:
%
%   x -- the solution
%
%   obj -- value of the objective function at the optimal point
%
%   info -- a structure containing the following fields:
%       info -- exit status
%           0 -- ok
%           1 -- reached maximum number of iteration
%           2 -- insufficient accuracy
%           3 -- internal error
%           4 -- insufficient memory
%           5 -- inconsistent constraints
%           6 -- unknown error
%       inconsistent_constraint_index -- index of the inconsistent
%           constraint if the exit status is equal to 5.
%
%   lambda -- Lagrange multipliers in the following order
%       1) equlity constraints
%       2) inequality constraints
%       3) lower bounds
%       4) upper bounds
%
%

function [x, obj, info, lambda] = qld(varargin)
    x0 = [];
    H = [];
    g = [];
    A = [];
    b = [];
    lb = [];
    ub = [];
    Ain = [];
    bin = [];

    options = [];

    x = [];
    lambda = [];
    obj = [];
    info.info = 6;


    QLD_INF = 1e30;


% Check number of input parameters
    if (2 <= nargin)
        x0 = varargin{1};
        H = varargin{2};

        if (3 == nargin)
            if (isstruct(varargin{3}))
                options = varargin{3};
            else
                g = varargin{3};
            end
        else
            if (5 <= nargin)
                g = varargin{3};
                A = varargin{4};
                b = varargin{5};

                if (6 == nargin)
                    options = varargin{6};
                else
                    lb = varargin{6};
                    ub = varargin{7};

                    if (8 == nargin)
                        options = varargin{8};
                    else
                        if (9 <= nargin)
                            Ain = varargin{8};
                            bin = varargin{9};

                            if (10 == nargin)
                                options = varargin{10};
                            end

                            if (10 < nargin)
                                error('Incorrect number of input variables.');
                            end
                        else
                            error('Incorrect number of input variables.');
                        end
                    end
                end
            else
                error('Incorrect number of input variables.');
            end
        end
    else
        error('Incorrect number of input variables.');
    end

% Check number of output parameters
    if (( 4 < nargout ) || ( 1 > nargout))
        error('Incorrect number of input variables.');
    end


% Check input correctness
    % objective
    if (isempty(H))
        error('Hessian is not initialized.');
    end

    if (size(H,1) ~= size(H,2))
        error('Hessian is not square.');
    end

    num_var = size(H, 1);


    if (isempty(g))
        g = zeros(num_var, 1);
    else
        if ( (size(g, 1) ~= num_var) || (size(g, 2) ~= 1) )
            error('Incorrect size of the constant vector.');
        end
    end


    % bounds
    if (isempty(lb))
        lb = -QLD_INF * ones(num_var, 1);
    else
        if ((size(lb, 1) ~= num_var) || (size(lb, 2) ~= 1))
            error('Incorrect size of the vector of lower bounds.');
        end
    end

    if (isempty(ub))
        ub =  QLD_INF * ones(num_var, 1);
    else
        if ((size(ub, 1) ~= num_var) || (size(ub, 2) ~= 1))
            error('Incorrect size of the vector of upper bounds.');
        end
    end


    % initial guess
    if (isempty(x0))
        x0 = zeros(num_var, 1);
    else
        if ((size(x0, 1) ~= num_var) && (size(x0, 2) ~= 1))
            error('Incorrect size of the initial guess.');
        end
    end


    % equality constraints
    if ( (isempty(A) && ~isempty(b)) || (~isempty(A) && isempty(b)))
        error('Partially initialized equality constraints.');
    end
    if (~isempty(b) && ~isempty(A))
        if ( (size(A, 1) ~= size(b, 1)) || (size(b, 2) ~= 1) || (size(A, 2) ~= num_var))
            error('Incorrect size of matrix or vector of equality constraints.');
        end
    end

    % inequality constraints
    if ( (isempty(Ain) && ~isempty(bin)) || (~isempty(Ain) && isempty(bin)) )
        error('Partially initialized inequality constraints.');
    end
    if (~isempty(bin) && ~isempty(Ain))
        if ( (size(Ain, 1) ~= size(bin, 1)) || (size(bin, 2) ~= 1) || (size(Ain, 2) ~= num_var) )
            error('Incorrect size of matrix or vector of inequality constraints.');
        end
    end


% Check options
    if (~isempty(options))
        if (~isstruct(options))
            error('Parameter "options" must be a struct.');
        end
    end


% Solve
    %minimize      1/2 x^ C x + c^x
    %subject to    a_j^x + b_j  =  0  ,  j=1,...,m_e
    %              a_j^x + b_j  >= 0  ,  j=m+1,...,m
    %              x_l <= x <= x_u

    b = -b;
    Ain = -Ain;

    num_eq = size(b, 1);

    A = [A; Ain];
    b = [b; bin];

    [x, info, lambda] = qld_interface(x0, H, g, lb, ub, A, b, num_eq, options);

    if (0 == info.info)
        obj = x'*H*x + g'*x;
    else
        obj = 0;
    end
end
