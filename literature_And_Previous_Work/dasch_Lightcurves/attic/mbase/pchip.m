function v = pchip(x,y,u)
% PCHIP  Piecewise Cubic Hermite Interpolating Polynomial.
% X is a row or column vector.  Y is a row or column vector of the same
% length as X, or a matrix with length(X) columns.
% YI = PCHIP(X,Y,XI) evaluates an interpolant at the elements of XI.
% PP = PCHIP(X,Y) returns a piecewise polynomial structure for use by PPVAL.
%
% The PCHIP interpolating function, P(x), satisfies:
%   On each subinterval,  x(k) <= x <= x(k+1),  P(x) is the cubic Hermite
%     interpolant to the given values and certain slopes at the two endpoints.
%   Therefore, P(x) interpolates y, i.e., P(x(j)) = y(j), and the first
%     derivative, P'(x), is continuous, but P''(x) is probably not continuous;
%     there may be jumps at the x(j).
%   The slopes at the x(j) are chosen in such a way that P(x) is
%     "shape preserving" and "respects monotonicity". This means that,
%     on intervals where the data are monotonic, so is P(x);
%     at points where the data have a local extremum, so does P(x).
%
% Comparing PCHIP with SPLINE:
%   The function S(x) supplied by SPLINE is constructed in exactly the same
%     way, except that the slopes at the x(j) are chosen differently, namely
%     to make even S''(x) continuous. This has the following effects.
%   Spline is smoother, i.e., S''(x) is continuous.
%   Spline is more accurate if the data are values of a smooth function.
%   Pchip has no overshoots and less oscillation if the data are not smooth.
%   Pchip is less expensive to set up.
%   The two are equally expensive to evaluate.
%
% Example:
%
%   x = -3:3;
%   y = [-1 -1 -1 0 1 1 1];
%   t = -3:.01:3;
%   plot(x,y,'o',t,[pchip(x,y,t); spline(x,y,t)])
%   legend('data','pchip','spline',4)
%
% See also: INTERP1, SPLINE, PPVAL.

% References:
%   F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
%   Interpolation", SIAM J. Numerical Analysis 17, 1980, 238-246.
%   David Kahaner, Cleve Moler and Stephen Nash, Numerical Methods
%   and Software, Prentice Hall, 1988.
%
%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.7 $  $Date: 2002/04/09 00:14:24 $

% Check dimensions.  Make x a monotonically increasing row vector.

if size(y,2) == 1, y = y.'; end
[m,n] = size(y);
if length(x) ~= n
   error('Y should have length(X) columns.')
end
if n <= 1
   error('There should be at least two data points.')
end
if ~isreal(x)
   error('The data abscissae should be real.')
end
x = x(:)';
h = diff(x);
if any(h < 0)
   [x,p] = sort(x);
   y = y(:,p);
   h = diff(x);
end
if any(h == 0)
   error('The data abscissae should be distinct.');
end

if nargin == 3

   if ~isreal(u)
      error('The interpolation points should be real.')
   end
   v = zeros(m*size(u,1),size(u,2));
   u = u(:)';
   q = length(u);
   if any(diff(u) < 0)
      [u,p] = sort(u);
   else
      p = 1:q;
   end

   % Find indices of subintervals, x(k) <= u < x(k+1).

   if isempty(u)
      k = u;
   else
      [ignore,k] = histc(u,x);
      k(u<x(1) | ~isfinite(u)) = 1;
      k(u>=x(n)) = n-1;
   end

   s = u - x(k);
   for r = 1:m

      % Compute slopes and other coefficients.
   
      del = diff(y(r,:))./h;
      d = pchipslopes(x,y(r,:),del);
      c = (3*del - 2*d(1:n-1) - d(2:n))./h;
      b = (d(1:n-1) - 2*del + d(2:n))./h.^2;
   
      % Evaluate interpolant.
   
      v(m*(p-1)+r) = y(r,k) + s.*(d(k) + s.*(c(k) + s.*b(k)));
   end
   
else
   
   % Generate piecewise polynomial structure.
   
   coefs = zeros(4,m*(n-1));
   for r = 1:m
      del = diff(y(r,:))./h;
      d = pchipslopes(x,y(r,:),del);
      j = r:m:m*(n-1);
      coefs(1,j) = (d(1:n-1) - 2*del + d(2:n))./h.^2;
      coefs(2,j) = (3*del - 2*d(1:n-1) - d(2:n))./h;
      coefs(3,j) = d(1:n-1);
      coefs(4,j) = y(r,1:n-1);
   end
   v.form = 'pp';
   v.breaks = x;
   v.coefs = coefs.';
   v.pieces = n-1;
   v.order = 4;
   v.dim = m;
end

% -------------------------------------------------------

function d = pchipslopes(x,y,del);
% PCHIPSLOPES  Derivative values for Piecewise Hermite Cubic Interpolation.
% d = pchipslopes(x,y,del) computes the first derivatives, d(k) = P'(x(k)).

%  Special case n=2, use linear interpolation.

   n = length(x);
   if n == 2  
      d = zeros(size(y));
      d(1) = del(1);
      d(2) = del(1);
      return
   end

%  Slopes at interior points.
%  d(k) = weighted average of del(k-1) and del(k) when they have the same sign.
%  d(k) = 0 when del(k-1) and del(k) have opposites signs or either is zero.

   d = zeros(size(y));
   if isreal(del)
      k = find(sign(del(1:n-2)).*sign(del(2:n-1)) > 0);
   else
      k = 1:n-2;
   end
   h = diff(x);
   hs = h(k)+h(k+1);
   w1 = (h(k)+hs)./(3*hs);
   w2 = (hs+h(k+1))./(3*hs);
   dmax = max(abs(del(k)), abs(del(k+1)));
   dmin = min(abs(del(k)), abs(del(k+1)));
   d(k+1) = dmin./conj(w1.*(del(k)./dmax) + w2.*(del(k+1)./dmax));

%  Slopes at end points.
%  Set d(1) and d(n) via non-centered, shape-preserving three-point formulae.

   d(1) = ((2*h(1)+h(2))*del(1) - h(1)*del(2))/(h(1)+h(2));
   if isreal(d) & (sign(d(1)) ~= sign(del(1)))
      d(1) = 0;
   elseif (sign(del(1)) ~= sign(del(2))) & (abs(d(1)) > abs(3*del(1)))
      d(1) = 3*del(1);
   end
   d(n) = ((2*h(n-1)+h(n-2))*del(n-1) - h(n-1)*del(n-2))/(h(n-1)+h(n-2));
   if isreal(d) & (sign(d(n)) ~= sign(del(n-1)))
      d(n) = 0;
   elseif (sign(del(n-1)) ~= sign(del(n-2))) & (abs(d(n)) > abs(3*del(n-1)))
      d(n) = 3*del(n-1);
   end
