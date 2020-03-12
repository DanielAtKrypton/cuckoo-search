function z = mantegna(alpha, c, n, matrix)
%MANTEGNA Stable random number generator.
% Z = MANTEGNA(ALPHA, C, n, N)
%
% Based on the method of R. N. Mantegna "Fast, accurate algorithm for
% numerical simulation of Lévy stable stochastic processes"
% Physical Review E 49 4677?83 (1994).
%
% alpha: defines the index of the distribution and controls the scale pro-
% perties of the stochastic process {z}.
%
% c: selects the scale unit of the process. It's the gamma letter inside
% Mantegna's paper.
%
% n: number of independent stochastic variables w calculated by equation
% (15).
%
% Z: Stochastic process.
matrixSize = size(matrix);
% Errortraps:
if (alpha < 0.3 || alpha > 1.99)
    disp('Valid range for alpha is [0.3;1.99].')
    z = NaN * zeros(matrixSize);
    return
end
if (c <= 0)
    disp('c must be positive.')
    z = NaN * zeros(matrixSize);
    return
end
if (n < 1)
    disp('n must be positive.')
    z = NaN * zeros(matrixSize);
    return
    
end
if nargin<4
    matrixSize = size(zeros(1,1));
end
if (matrixSize(2) <= 0)
    disp('N must be positive.')
    z = NaN;
    return
end
invalpha = 1/alpha;
%% sigmaX calculation in function of alpha. Equation (12).
sigx = ((gamma(1+alpha)*sin(pi*alpha/2))/(gamma((1+alpha)/2)...
    *alpha*2^((alpha-1)/2)))^invalpha;
%% v calculation in function of alpha. Equation (6).
sigy = 1;
x = sigx*randn(matrixSize);
y = sigy*randn(matrixSize);
v = x./abs(y).^invalpha;
%% kappa calculation in function of alpha. Equation (20).
kappa = (alpha*gamma((alpha+1)/(2*alpha)))/gamma(invalpha)...
    *((alpha*gamma((alpha+1)/2))/(gamma(1+alpha)*sin(pi*alpha/2)))^invalpha;
%% C estimation
p = [-17.7767 113.3855 -281.5879 337.5439 -193.5494 44.8754];
C_a = polyval(p, alpha);
w = ((kappa-1)*exp(-v/C_a)+1).*v; % equation (15)
if(n>1)
    z = (1/n^invalpha)*sum(w);
else
    z = w;
end
z = c^invalpha*z;