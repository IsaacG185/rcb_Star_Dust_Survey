function varargout = strread(varargin);
%STRREAD Read formatted data from string.
%    A = STRREAD('STRING')
%    A = STRREAD('STRING','',N)
%    A = STRREAD('STRING','',param,value, ...)
%    A = STRREAD('STRING','',N,param,value, ...) reads numeric data from
%    the STRING into a single variable.  If the string contains any text data,
%    an error is produced.
%
%    [A,B,C, ...] = STRREAD('STRING','FORMAT')
%    [A,B,C, ...] = STRREAD('STRING','FORMAT',N)
%    [A,B,C, ...] = STRREAD('STRING','FORMAT',param,value, ...)
%    [A,B,C, ...] = STRREAD('STRING','FORMAT',N,param,value, ...) reads
%    data from the STRING into the variables A,B,C,etc.  The type of each
%    return argument is given by the FORMAT string.  The number of return
%    arguments must match the number of conversion specifiers in the FORMAT
%    string.  If there are fewer fields in the string than matching conversion
%    specifiers in the format string, an error is produced.
%
%    If N is specified, the format string is reused N times.  If N is -1 (or
%    not specified) STRREAD reads the entire string.
%
%    Example
%
%      s = sprintf('a,1,2\nb,3,4\n');
%      [a,b,c] = strread(s,'%s%d%d','delimiter',',')
%
%   See TEXTREAD for more examples and definition of terms.
%
%   See also TEXTREAD, SSCANF, FILEFORMATS.

%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.6 $ $Date: 2002/06/05 20:10:15 $

%   Implemented as a mex file.

% do some preliminary error checking
error(nargchk(1,inf,nargin));

% allow empty string to pass through untouched
if isempty(varargin{1})
    return;
end

if nargout == 0
    nlhs = 1;
else
    nlhs = nargout;
end

[varargout{1:nlhs}]=dataread('string',varargin{:});
