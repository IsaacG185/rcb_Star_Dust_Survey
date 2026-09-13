% Copyright the President and Fellows of Harvard College.
% Licensed under the MIT License

function [allfcns,msg] = fprefcnchk(funstr,caller,lenVarIn,gradflag)

%%funstr
%%caller
%%lenVarIn
%%gradflag

%PREFCNCHK Pre- and post-process function expression for FUNCHK.
%   [ALLFCNS,MSG] = PREFUNCHK(FUNSTR,CALLER,lenVarIn,GRADFLAG) takes
%   the (nonempty) expression FUNSTR from CALLER with LenVarIn extra arguments,
%   parses it according to what CALLER is, then returns a string or inline
%   object in ALLFCNS.  If an error occurs, this message is put in MSG.
%
%   ALLFCNS is a cell array: 
%    ALLFCNS{1} contains a flag 
%    that says if the objective and gradients are together in one function 
%    (calltype=='fungrad') or in two functions (calltype='fun_then_grad')
%    or there is no gradient (calltype=='fun'), etc.
%    ALLFCNS{2} contains the string CALLER.
%    ALLFCNS{3}  contains the objective function
%    ALLFCNS{4}  contains the gradient function (transpose of Jacobian).
%  
%    NOTE: we assume FUNSTR is nonempty.
% Initialize
msg='';
allfcns = {};
funfcn = [];
gradfcn = [];

if gradflag
    calltype = 'fungrad';
else
    calltype = 'fun';
end

% {fun}
isa(funstr,'cell')
if isa(funstr,'cell') 
  if length(funstr)==1
    % take the cellarray apart: we know it is nonempty
    if gradflag
        calltype = 'fungrad';
    end
    [funfcn, msg] = fcnchk(funstr{1},lenVarIn);
    if ~isempty(msg)
        error(msg);
    end
    
    % {fun,[]}      
  elseif length(funstr)==2 & isempty(funstr{2})
    if gradflag
        calltype = 'fungrad';
    end
    [funfcn, msg] = fcnchk(funstr{1},lenVarIn);
    if ~isempty(msg)
        error(msg);
    end  
    
    % {fun, grad}   
  elseif length(funstr)==2 % and ~isempty(funstr{2})
    
    [funfcn, msg] = fcnchk(funstr{1},lenVarIn);
    if ~isempty(msg)
        error(msg);
    end  
    [gradfcn, msg] = fcnchk(funstr{2},lenVarIn);
    if ~isempty(msg)
        error(msg);
    end
    calltype = 'fun_then_grad';
    if ~gradflag
        warnstr = ...
            sprintf(['Jacobian function provided but FITOPTIONS.Jacobian=''off'';\n', ...
            '  ignoring Jacobian function and using finite-differencing.\n', ...
            '  Rerun with FITOPTIONS.Jacobian=''on'' to use Jacobian function.\n']);
        warning(warnstr);
        calltype = 'fun';
    end   
  else
    errmsg = sprintf(['FUN must be a function object or an inline object;\n', ...
        ' or, FUN may be a cell array that contains these type of objects.']);
    error(errmsg)
  end
  
else  %Not a cell; is a string expression, function name string, function handle, or inline object
funfcn = funstr;
%%% [funfcn, msg] = fcnchk(funstr,lenVarIn);
%%% if ~isempty(msg)
%%%     error(msg);
%%% end   
    if gradflag % gradient and function in one function/M-file
        gradfcn = funfcn; % Do this so graderr will print the correct name
    end  
end

allfcns{1} = calltype;
allfcns{2} = caller;
allfcns{3} = funfcn;
allfcns{4} = gradfcn;
allfcns{5}=[];
allfcns
isa(allfcns,'cell')
endfunction;


        FUN=inline('c(1)*(x-c(2)).^2+c(3)','c','x');
        x = [0.400000,  -0.040750,   0.153544];
YDATA = [...
0.15402,0.15402,0.16712,0.16179,0.15714,0.15366,0.16743,0.15493,...
0.16911,0.21429,0.16509,0.15667,0.17583,0.15350,0.15459,0.15849,...
0.16010,0.19695,0.16050,0.17868,0.20193,0.19467,0.18327,0.16269,...
0.15526,0.15350,0.15871,0.15647,0.15436,0.15710,0.15546,0.15972,...
0.17885,0.15821,0.16522,0.15513,0.15376,0.15349,0.15351,0.15609,...
0.17658,0.16171,0.19565,0.18052,0.17118,0.21807,0.15991,0.18906,...
0.17870,0.18346];


   LB = [];
   UB = [];
   options = [];
defaultopt = struct('Display','final','LargeScale','on', ...
   'TolX',1e-6,'TolFun',1e-6,'DerivativeCheck','off',...
   'Diagnostics','off',...
   'Jacobian','off','JacobMult',[],...% JacobMult set to [] by default
   'JacobPattern','sparse(ones(Jrows,Jcols))',...
   'MaxFunEvals','100*numberOfVariables',...
   'DiffMaxChange',1e-1,'DiffMinChange',1e-8,...
   'PrecondBandWidth',0,'TypicalX','ones(numberOfVariables,1)',...
   'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
   'TolPCG',0.1,'MaxIter',400,...
   'LineSearchType','quadcubic','LevenbergMarquardt','on'); 
caller = 'lsqcurvefit';
computeLambda = 1;
zzz  = ...
[0.00000, 0.00000, 0.16069, 0.11231, 0.06025,-0.05938, 0.16307, 0.02314,...
 0.17307, 0.41340, 0.14456,-0.12820,-0.29028,-0.04455,-0.08834,-0.15282,...
-0.17125,-0.40694,-0.17550,-0.30817,-0.43090,-0.39571,-0.33547,-0.19687,...
 0.02988,-0.03283, 0.08015, 0.05051, 0.00987, 0.05972, 0.03373, 0.09135,...
 0.23726, 0.07411, 0.14558,-0.10101,-0.01054,-0.03622,-0.03018,-0.11863,...
 0.22354,-0.18756, 0.32720, 0.24720, 0.18783, 0.42957,-0.16920,-0.36702,...
 0.23639,-0.33657];


varargin{1:1} = [zzz];

%%%FUN
%%%x
%%%YDATA
%%%LB
%%%UB
%%%options
%%%defaultopt
%%%caller
%%%computeLambda
%%%varargin



%LSQNCOMMON Solves non-linear least squares problems.
%   [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN]=...
%      LSQNCOMMON(FUN,X0,YDATA,LB,UB,OPTIONS,DEFAULTOPT,CALLER,COMPUTELAMBDA,XDATA,VARARGIN...) 
%   contains all the setup code common to both LSQNONLIN and LSQCURVEFIT to call either the 
%   large-scale SNLS or the medium-scale NLSQ.

%   Copyright 2001-2002 The MathWorks, Inc. 
%   $Revision: 1.11 $  $Date: 2002/03/12 20:26:38 $

% Note: none of these warnings should go to the command line. They've
% all been taken care of in FIT. If they appear from here (not FIT), 
% it is an error.

msg = '';
xstart=x(:);
numberOfVariables=length(xstart);
% Note: XDATA is bundled with varargin already for lsqcurvefit 
lenVarIn = length(varargin);

large = 'large-scale';
medium = 'medium-scale';

% Options setup
switch optimget(options,'Display',defaultopt,'fast')
case {'notify'}
   verbosity = 1;
case {'off','none'}
   verbosity = 0;
case 'iter'
   verbosity = 3;
case 'final'
   verbosity = 2;
case 'testing'
    verbosity = Inf;
otherwise
   verbosity = 2;
end

l = LB; u = UB;
lFinite = ~isinf(l);
uFinite = ~isinf(u);
l
u
xstart


%%%if min(min(u-xstart),min(xstart-l)) < 0
%%%    xstart = startx(u,l); 
%%%end

diagnostics = isequal(optimget(options,'Diagnostics',defaultopt,'fast'),'on');
gradflag =  strcmp(optimget(options,'Jacobian',defaultopt,'fast'),'on');
line_search = strcmp(optimget(options,'LargeScale',defaultopt,'fast'),'off'); % 0 means large-scale, 1 means medium-scale
mtxmpy = optimget(options,'JacobMult',[]); % use old



if isequal(mtxmpy,'atamult')
    warnstr = sprintf(['Potential function name clash with a Toolbox helper function:\n',...
                      'Use a name besides ''atamult'' for your JacobMult function to\n',...
                      'avoid errors or unexpected results.\n']);
    warning(warnstr)
end

% Convert to inline function as needed
if ~isempty(FUN)  % will detect empty string, empty matrix, empty cell array
    [funfcn, msg] = fprefcnchk(FUN,caller,lenVarIn,gradflag);
    
else
    errmsg = sprintf(['FUN must be a function or an inline object;\n', ...
        ' or, FUN may be a cell array that contains these type of objects.']);
    error(errmsg)
end

fuser = [];  
JAC = [];
x(:) = xstart;
try
    switch funfcn{1}
    case 'fun'
        fuser = feval(funfcn{3},x,varargin{:});
    case 'fungrad'
        [fuser,JAC] = feval(funfcn{3},x,varargin{:});
    case 'fun_then_grad'
        fuser = feval(funfcn{3},x,varargin{:}); 
        JAC = feval(funfcn{4},x,varargin{:});
    otherwise
        errmsg = sprintf('Undefined calltype in %s\n',upper(caller));
        error(errmsg);
    end
catch
    if ~(isequal(funfcn{1},'fun_then_grad')) | isempty(fuser)
        badfunfcn = funfcn{3};
    else % 'fun_then_grad & ~isempty(fuser) (so it error'ed out on the JAC call)
        badfunfcn = funfcn{4};
    end   
    if isa(badfunfcn,'inline')
        errmsg = sprintf(['User supplied expression or inline function ==> %s\n' ...
                'failed with the following error:\n\n%s'],...
                 formula(badfunfcn),lasterr);
    else % function (i.e. string name of), function handle, inline function, or other
        % save last error in case call to "char" errors out
        lasterrmsg = lasterr;
        % Convert to char if possible
        try
            charOfFunction = char(badfunfcn);
        catch
            charOfFunction = '';
        end
        if ~isempty(charOfFunction)
            errmsg = sprintf(['User supplied function ==> %s\n' ...
                    'failed with the following error:\n\n%s'],...
                     charOfFunction,lasterrmsg);
        else
            errmsg = sprintf(['User supplied function\n' ...
                    'failed with the following error:\n\n%s'],...
                    lasterrmsg);
        end
    end
    error(errmsg)
end

if isequal(caller,'lsqcurvefit')
    if ~isequal(size(fuser), size(YDATA))
        error('Function value and YDATA sizes are incommensurate.')
    end
    fuser = fuser - YDATA;  % preserve fuser shape until after subtracting YDATA 
end

f = fuser(:);
nfun=length(f);

if gradflag
    % check size of JAC
    [Jrows, Jcols]=size(JAC);
    if isempty(mtxmpy) 
        % Not using 'JacobMult' so Jacobian must be correct size
        if Jrows~=nfun | Jcols ~=numberOfVariables
            rrstr = sprintf(['User-defined Jacobian is not the correct size:\n',...
                '    the Jacobian matrix should be %d-by-%d\n'],nfun,numberOfVariables);
            error(errstr);
        end
    end
else
    Jrows = nfun; 
    Jcols = numberOfVariables;   
end

% trustregion and enough equations (as many as variables) 
if ~line_search & nfun >= numberOfVariables 
    OUTPUT.algorithm = large;
    
    % trust region and not enough equations -- switch to line_search
elseif ~line_search & nfun < numberOfVariables 
   warnstr = sprintf(['Large-scale method requires at least as many equations as variables; \n',...
        '   switching to line-search method instead.  Upper and lower bounds will be ignored.\n']);
    warning(warnstr);
    OUTPUT.algorithm = medium;
    
    % line search and no bounds  
elseif line_search & isempty(l(lFinite)) & isempty(u(uFinite))
    OUTPUT.algorithm = medium;
    
    % line search and  bounds  and enough equations, switch to trust region 
elseif line_search & (~isempty(l(lFinite)) | ~isempty(u(uFinite))) & nfun >= numberOfVariables
    warnstr = sprintf(['Levenberg-Marquardt and Gauss-Newton algorithms do not handle\n',...
        'bound constraints; switching to trust-region algorithm instead.\n']);
    warning(warnstr);
    OUTPUT.algorithm = large;
    
    % can't handle this one:   
elseif line_search & (~isempty(l(lFinite)) | ~isempty(u(uFinite)))  & nfun < numberOfVariables
    errstr = sprintf(['Line-search method does not handle bound constraints \n',...
        '   and large-scale method requires at least as many equations as variables; \n',...
        '   aborting.\n']);

    error(errstr);
end

% Execute algorithm
if isequal(OUTPUT.algorithm,large)
    if ~gradflag % provide sparsity of Jacobian if not provided.
        Jstr = optimget(options,'JacobPattern',[]);
        if isempty(Jstr)  
            % Put this code separate as it might generate OUT OF MEMORY error
            Jstr = sparse(ones(Jrows,Jcols));
        end
        if ischar(Jstr) 
            if isequal(lower(Jstr),'sparse(ones(jrows,jcols))')
                Jstr = sparse(ones(Jrows,Jcols));
            else
                error('Option ''JacobPattern'' must be a matrix if not the default.')
            end
        end
    else
        Jstr = [];
    end
    warnstate = warning('off');
    [x,FVAL,LAMBDA,JACOB,EXITFLAG,OUTPUT,msg]=...
        snls(funfcn,x,l,u,verbosity,options,defaultopt,f,JAC,YDATA,caller,Jstr,computeLambda,varargin{:});
    warning(warnstate);
    % So nlsq and snls return the same thing when no bounds exist:
    if all(~isfinite(l)) & all(~isfinite(u))
        LAMBDA.upper=[]; LAMBDA.lower=[];   
    end        
else
    warnstate = warning('off');
    [x,FVAL,JACOB,EXITFLAG,OUTPUT,msg] = ...
        nlsq(funfcn,x,verbosity,options,defaultopt,f,JAC,YDATA,caller,varargin{:});
    LAMBDA.upper=[]; LAMBDA.lower=[];   
    warning(warnstate);
end
Resnorm = FVAL'*FVAL;
if verbosity > 1 | ( verbosity > 0 & EXITFLAG <=0 )
    disp(msg);
end


% Reset FVAL to shape of the user-function, fuser
FVAL = reshape(FVAL,size(fuser));


printf("--end of lsqncommon--\n");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

