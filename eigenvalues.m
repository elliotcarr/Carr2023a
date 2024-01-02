function lambda = eigenvalues(aL,bL,R,N)
% Solves equation (3.14) from paper for first N positive eigenvalues.

f = @(lambda) bL*R*lambda.*cos(lambda*R) - bL*sin(lambda*R) + R*aL*sin(lambda*R);
lambda = zeros(N,1);
maxiters = 1000;

% Bisection method
for n = 1:N
    a = (2*n-1)*pi/(2*R); b = (2*n+1)*pi/(2*R);
    if f(a)*f(b) > 0
        error('Inverval does not bracket root');
    end
    i = 0;
    while (b-a)/a > 1e-13 && i < maxiters% changed from absolute (b-a) to relative
        c = (a+b)/2; fc = f(c); i = i + 1;
        if sign(fc) == sign(f(a))
            a = c;
        else
            b = c;
        end
    end
    lambda(n) = c;
end

if i == maxiters && (b-a) > 1e-13
    error('Eigenvalue didn''t converge.\n')
end