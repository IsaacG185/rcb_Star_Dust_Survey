function [f,msg] = fcnchk(fun,varargin)
%FCNCHK Check FUNFUN function argument.
%   FCNCHK(FUN,...) returns an inline object based on FUN if FUN
%   is a string containing parentheses, variables, and math
%   operators.  FCNCHK simply returns FUN if FUN is a function handle, 
%   a name string (e.g. 'sin'), or a matlab object with an feval method 
%   (such as an inline object). 
%
%   FCNCHK is a helper function for FMINBND, FMINSEARCH, FZERO, etc. so they
%   can compute with string expressions in addition to m-file functions.
%
%   FCNCHK(FUN,...,'vectorized') processes the string (e.g., replacing
%   '*' with '.*') to produce a vectorized function.
%
%   When FUN contains an expression then FCNCHK(FUN,...) is the same as
%   INLINE(FUN,...) except that the optional trailing argument 'vectorized'
%   can be used to produce a vectorized function.
%
%   [F,MSG] = FCNCHK(...) returns an empty string in MSG if successful
%   or an error message string if not.
%
%   See also INLINE, @.

%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.29 $  $Date: 2002/04/08 20:26:45 $
%   Jan 23, 2012 Edward J. Los - fix short-circuit operators

msg = '';

nin = nargin;
if (nin>1) && strcmp(varargin{end},'vectorized')
    vectorizing = 1;
    nin = nin-1;
else
    vectorizing = 0;
end

if ischar(fun)
    fun = strtrim(fun);
    % Check for non-alphanumeric characters that must be part of an
    % expression.
    if isempty(fun),
        f = inline('[]');
    elseif ~vectorizing && isidentifier(fun)
        f = fun; % Must be a function name only
        if isequal('x',deblank(fun))
            warning('MATLAB:fcnchk:AmbiguousX', ...
                ['Ambiguous expression or function input.\n The string ''x'' will be ',...
                    'interpreted as the name of a ',...
                'function called ''x'' \n (e.g., x.m) and not as the mathematical expression ''x'' (i.e., f(x)=x). \n ',...
                'Use inline(''x'') if you meant the mathematical expression ''x''.']);
        end
    else
        if vectorizing
            f = inline(vectorize(fun),varargin{1:nin-1});
            var = argnames(f);
            f = inline([formula(f) '.*ones(size(' var{1} '))'],var{1:end});
        else
            f = inline(fun,varargin{1:nin-1});
        end 
    end
elseif isa(fun,'function_handle') 
    f = fun; 
    % is it a matlab object with a feval method?
elseif isobject(fun)
    % delay the methods call unless we know it is an object to avoid runtime error for compiler
    meths = methods(class(fun));
    if any(strmatch('feval',meths,'exact'))
       if vectorizing && any(strmatch('vectorize',meths,'exact'))
          f = vectorize(fun);
       else
          f = fun;
       end
    else % no feval method
        f = '';
        msg = ['If FUN is a MATLAB object, it must have an feval method.'];
    end
else
    f = '';
    msg = ['FUN must be a function, a valid string expression, ', ...
            sprintf('\n'),'or an inline function object.'];
end
if nargout < 2, error(msg); end


%------------------------------------------
function s1 = strtrim(s)
%STRTRIM Trim spaces from string.

if ~isempty(s) && ~isstr(s)
    warning('MATLAB:fcnchk:InputNotString','Input must be a string.')
end

if isempty(s)
    s1 = s;
else
    % remove leading and trailing blanks (including nulls)
    c = find(s ~= ' ' & s ~= 0);
    s1 = s(min(c):max(c));
end

%-------------------------------------------
function tf = isidentifier(str)

tf = 0;

if ~isempty(str)
    first = str(1);
    if (isletter(first))
        letters = isletter(str);
        numerals = (48 <= str) & (str <= 57);
        underscore = (95 == str);
        if (all(letters | numerals | underscore))
            tf = 1;
        end
    end
end

tf = logical(tf);
