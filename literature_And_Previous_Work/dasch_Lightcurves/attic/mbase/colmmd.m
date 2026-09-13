function p = colmmd(S)
%COLMMD Column minimum degree permutation.
%   P = COLMMD(S) returns the column minimum degree permutation vector
%   for the sparse matrix S. For a non-symmetric matrix S, S(:,P)
%   tends to have sparser LU factors than S.
%
%   See also SYMMMD, SYMRCM, COLPERM.

%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 5.11 $  $Date: 2002/04/09 00:26:11 $

p = sparsfun('colmmd',S);
[ignore,q] = sparsfun('coletree',S(:,p));
p = p(q);
