function TF = isobject(a)
%ISOBJECT True for MATLAB objects.
%   ISOBJECT(A) returns 1 if A is a MATLAB object and 0 otherwise. 
%
%   See also ISCELL, ISSTRUCT, ISNUMERIC, ISJAVA.

%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.19 $  $Date: 2002/04/08 19:22:48 $

TF = ~isnumeric(a) && ~islogical(a) && ~ischar(a) && ...
     ~isa(a,'cell') && ~isa(a,'struct') && ...
     ~isa(a,'opaque') && ~isa(a,'function_handle');
