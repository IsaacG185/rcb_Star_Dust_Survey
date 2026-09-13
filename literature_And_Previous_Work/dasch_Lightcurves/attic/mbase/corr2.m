function r = corr2(varargin)
%CORR2 Compute 2-D correlation coefficient.
%   R = CORR2(A,B) computes the correlation coefficient between A
%   and B, where A and B are matrices or vectors of the same size.
%
%   Class Support
%   -------------
%   A and B can be numeric or logical. 
%   R is a scalar double.
%
%   See also CORRCOEF, STD2.

%   Copyright 1993-2002 The MathWorks, Inc.  
%   $Revision: 5.18 $  $Date: 2002/03/15 15:26:54 $

[a,b] = ParseInputs(varargin{:});

a = a - mean2(a);
b = b - mean2(b);
r = sum(sum(a.*b))/sqrt(sum(sum(a.*a))*sum(sum(b.*b)));

%--------------------------------------------------------
function [A,B] = ParseInputs(varargin)

checknargin(2,2,nargin, mfilename);

A = varargin{1};
B = varargin{2};

checkinput(A, {'logical' 'numeric'}, 'real', mfilename, 'A', 1);
checkinput(B, {'logical' 'numeric'}, 'real', mfilename, 'B', 2);

if any(size(A)~=size(B))
    messageId = 'Images:corr2:notSameSize';
    message1 = 'A and B must be the same size.';
    error(messageId, '%s', message1);
end

if (~isa(A,'double'))
    A = double(A);
end

if (~isa(B,'double'))
    B = double(B);
end










