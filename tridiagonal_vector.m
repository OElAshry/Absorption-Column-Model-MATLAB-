%   --> The vectors a, b, and c define the diagonals of the tridiagonal
%       matrix, A, as shown below:
%
%               ⌈ ⋱ ⋱       ⌉
%               | ⋱ ⋱ c     |
%           A = |   ⋱ b ⋱   |
%               |     a ⋱ ⋱ |
%               ⌊       ⋱ ⋱ ⌋

% Reference: https://github.com/tamaskis/tridiagonal-MATLAB/blob/main/tridiagonal_vector.m

function x = tridiagonal_vector(a,b,c,d)
    
    % determines n
    n = length(d);
    
    % preallocates solution vector
    x = zeros(n,1);
    
    % forward elimination
    for i = 2:n
        w = a(i-1)/b(i-1);
        b(i) = b(i)-w*c(i-1);
        d(i) = d(i)-w*d(i-1);
        
    end
    
    % backward substitution
    x(n) = d(n)/b(n);
    for i = (n-1):(-1):1
        x(i) = (d(i)-c(i)*x(i+1))/b(i);
    end
    
end
