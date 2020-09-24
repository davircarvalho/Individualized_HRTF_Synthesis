function [X,ndx,dbg] = natsort(X,rgx,varargin)
% Alphanumeric / Natural-Order sort the strings in a cell array of strings (1xN char).
%
% (c) 2012-2019 Stephen Cobeldick
%
% Alphanumeric sort a cell array of strings: sorts by character order and
% also by the values of any number substrings. Default: match all integer
% number substrings and perform a case-insensitive ascending sort.
%
%%% Example:
% >> X = {'x2', 'x10', 'x1'};
% >> sort(X)
% ans =   'x1'  'x10'  'x2'
% >> natsort(X)
% ans =   'x1'  'x2'  'x10'
%
%%% Syntax:
%  Y = natsort(X)
%  Y = natsort(X,rgx)
%  Y = natsort(X,rgx,<options>)
% [Y,ndx,dbg] = natsort(X,...)
%
% To sort filenames or filepaths use NATSORTFILES (FEX 47434).
% To sort the rows of a cell array of strings use NATSORTROWS (FEX 47433).
%
%% Number Substrings %%
%
% By default consecutive digit characters are interpreted as an integer.
% Specifying the optional regular expression pattern allows the numbers to
% include a +/- sign, decimal digits, exponent E-notation, quantifiers,
% or look-around matching. For information on defining regular expressions:
% http://www.mathworks.com/help/matlab/matlab_prog/regular-expressions.html
%
% The number substrings are parsed by SSCANF into numeric values, using
% either the *default format '%f' or the user-supplied format specifier.
%
% This table shows examples of regular expression patterns for some common
% notations and ways of writing numbers, with suitable SSCANF formats:
%
% Regular       | Number Substring | Number Substring              | SSCANF
% Expression:   | Match Examples:  | Match Description:            | Format Specifier:
% ==============|==================|===============================|==================
% *         \d+ | 0, 123, 4, 56789 | unsigned integer              | %f  %i  %u  %lu
% --------------|------------------|-------------------------------|------------------
%      [-+]?\d+ | +1, 23, -45, 678 | integer with optional +/- sign| %f  %i  %d  %ld
% --------------|------------------|-------------------------------|------------------
%     \d+\.?\d* | 012, 3.45, 678.9 | integer or decimal            | %f
% (\d+|Inf|NaN) | 123, 4, NaN, Inf | integer, Inf, or NaN          | %f
%  \d+\.\d+e\d+ | 0.123e4, 5.67e08 | exponential notation          | %f
% --------------|------------------|-------------------------------|------------------
%  0[0-7]+      | 012, 03456, 0700 | octal notation & prefix       | %o  %i
%   [0-7]+      |  12,  3456,  700 | octal notation                | %o
% --------------|------------------|-------------------------------|------------------
%  0X[0-9A-F]+  | 0X0, 0X3E7, 0XFF | hexadecimal notation & prefix | %x  %i
%    [0-9A-F]+  |   0,   3E7,   FF | hexadecimal notation          | %x
% --------------|------------------|-------------------------------|------------------
%  0B[01]+      | 0B1, 0B101, 0B10 | binary notation & prefix      | %b   (not SSCANF)
%    [01]+      |   1,   101,   10 | binary notation               | %b   (not SSCANF)
% --------------|------------------|-------------------------------|------------------
%
%% Debugging Output Array %%
%
% The third output is a cell array <dbg>, to check if the numbers have
% been matched by the regular expression <rgx> and converted to numeric
% by the SSCANF format. The rows of <dbg> are linearly indexed from <X>,
% even columns contain numbers, odd columns contain split substrings:
%
% >> [~,~,dbg] = natsort(X)
% dbg =
%    'x'    [ 2]
%    'x'    [10]
%    'x'    [ 1]
%
%% Examples %%
%
%%% Multiple integers (e.g. release version numbers):
% >> A = {'v10.6', 'v9.10', 'v9.5', 'v10.10', 'v9.10.20', 'v9.10.8'};
% >> sort(A)
% ans =   'v10.10'  'v10.6'  'v9.10'  'v9.10.20'  'v9.10.8'  'v9.5'
% >> natsort(A)
% ans =   'v9.5'  'v9.10'  'v9.10.8'  'v9.10.20'  'v10.6'  'v10.10'
%
%%% Integer, decimal, NaN, or Inf numbers, possibly with +/- signs:
% >> B = {'test+NaN', 'test11.5', 'test-1.4', 'test', 'test-Inf', 'test+0.3'};
% >> sort(B)
% ans =   'test' 'test+0.3' 'test+NaN' 'test-1.4' 'test-Inf' 'test11.5'
% >> natsort(B, '[-+]?(NaN|Inf|\d+\.?\d*)')
% ans =   'test' 'test-Inf' 'test-1.4' 'test+0.3' 'test11.5' 'test+NaN'
%
%%% Integer or decimal numbers, possibly with an exponent:
% >> C = {'0.56e007', '', '43E-2', '10000', '9.8'};
% >> sort(C)
% ans =   ''  '0.56e007'  '10000'  '43E-2'  '9.8'
% >> natsort(C, '\d+\.?\d*([eE][-+]?\d+)?')
% ans =   ''  '43E-2'  '9.8'  '10000'  '0.56e007'
%
%%% Hexadecimal numbers (with '0X' prefix):
% >> D = {'a0X7C4z', 'a0X5z', 'a0X18z', 'a0XFz'};
% >> sort(D)
% ans =   'a0X18z'  'a0X5z'  'a0X7C4z'  'a0XFz'
% >> natsort(D, '0X[0-9A-F]+', '%i')
% ans =   'a0X5z'  'a0XFz'  'a0X18z'  'a0X7C4z'
%
%%% Binary numbers:
% >> E = {'a11111000100z', 'a101z', 'a000000000011000z', 'a1111z'};
% >> sort(E)
% ans =   'a000000000011000z'  'a101z'  'a11111000100z'  'a1111z'
% >> natsort(E, '[01]+', '%b')
% ans =   'a101z'  'a1111z'  'a000000000011000z'  'a11111000100z'
%
%%% Case sensitivity:
% >> F = {'a2', 'A20', 'A1', 'a10', 'A2', 'a1'};
% >> natsort(F, [], 'ignorecase') % default
% ans =   'A1'  'a1'  'a2'  'A2'  'a10'  'A20'
% >> natsort(F, [], 'matchcase')
% ans =   'A1'  'A2'  'A20'  'a1'  'a2'  'a10'
%
%%% Sort order:
% >> G = {'2', 'a', '', '3', 'B', '1'};
% >> natsort(G, [], 'ascend') % default
% ans =   ''   '1'  '2'  '3'  'a'  'B'
% >> natsort(G, [], 'descend')
% ans =   'B'  'a'  '3'  '2'  '1'  ''
% >> natsort(G, [], 'num<char') % default
% ans =   ''   '1'  '2'  '3'  'a'  'B'
% >> natsort(G, [], 'char<num')
% ans =   ''   'a'  'B'  '1'  '2'  '3'
%
%%% UINT64 numbers (with full precision):
% >> natsort({'a18446744073709551615z', 'a18446744073709551614z'}, [], '%lu')
% ans =       'a18446744073709551614z'  'a18446744073709551615z'
%
%% Input and Output Arguments %%
%
%%% Inputs (*==default):
% X   = CellArrayOfCharRowVectors, to be sorted into natural-order.
% rgx = Regular expression to match number substrings, '\d+'*
%     = [] uses the default regular expression, which matches integers.
% <options> can be entered in any order, as many as required:
%     = Sort direction: 'descend'/'ascend'*
%     = NaN/number order: 'NaN<num'/'num<NaN'*
%     = Character/number order: 'char<num'/'num<char'*
%     = Character case handling: 'matchcase'/'ignorecase'*
%     = SSCANF number conversion format, e.g.: '%f'*, '%x', '%li', '%b', etc.
%
%%% Outputs:
% Y   = CellArrayOfCharRowVectors, <X> sorted into natural-order.
% ndx = NumericArray, such that Y = X(ndx). The same size as <X>.
% dbg = CellArray of the parsed characters and number values.
%       Each row is one input char vector, linear-indexed from <X>.
%
% See also SORT NATSORTFILES NATSORTROWS CELLSTR REGEXP IREGEXP SSCANF
%% Input Wrangling %%
%
assert(iscell(X),'First input <X> must be a cell array.')
tmp = cellfun('isclass',X,'char') & cellfun('size',X,1)<2 & cellfun('ndims',X)<3;
assert(all(tmp(:)),'First input <X> must be a cell array of char row vectors (1xN char).')
%
if nargin<2 || isnumeric(rgx)&&isempty(rgx)
    rgx = '\d+';
else
    assert(ischar(rgx)&&ndims(rgx)<3&&size(rgx,1)==1,...
        'Second input <rgx> must be a regular expression (char row vector).') %#ok<ISMAT>
end
%
% Optional arguments:
tmp = cellfun('isclass',varargin,'char') & cellfun('size',varargin,1)<2 & cellfun('ndims',varargin)<3;
assert(all(tmp(:)),'All optional arguments must be char row vectors (1xN char).')
% Character case:
ccm = strcmpi(varargin,'matchcase');
ccx = strcmpi(varargin,'ignorecase')|ccm;
% Sort direction:
sdd = strcmpi(varargin,'descend');
sdx = strcmpi(varargin,'ascend')|sdd;
% Char/num order:
chb = strcmpi(varargin,'char<num');
chx = strcmpi(varargin,'num<char')|chb;
% NaN/num order:
nab = strcmpi(varargin,'NaN<num');
nax = strcmpi(varargin,'num<NaN')|nab;
% SSCANF format:
sfx = ~cellfun('isempty',regexp(varargin,'^%([bdiuoxfeg]|l[diuox])$'));
%
nsAssert(1,varargin,sdx,'Sort direction')
nsAssert(1,varargin,chx,'Char<->num')
nsAssert(1,varargin,nax,'NaN<->num')
nsAssert(1,varargin,sfx,'SSCANF format')
nsAssert(0,varargin,~(ccx|sdx|chx|nax|sfx))
%
% SSCANF format:
if nnz(sfx)
    fmt = varargin{sfx};
    if strcmpi(fmt,'%b')
        cls = 'double';
    else
        cls = class(sscanf('0',fmt));
    end
else
    fmt = '%f';
    cls = 'double';
end
%
%% Identify Numbers %%
%
[mat,spl] = regexpi(X(:),rgx,'match','split',varargin{ccx});
%
% Determine lengths:
nmx = numel(X);
nmn = cellfun('length',mat);
nms = cellfun('length',spl);
mxs = max(nms);
%
% Preallocate arrays:
bon = bsxfun(@le,1:mxs,nmn).';
bos = bsxfun(@le,1:mxs,nms).';
arn = zeros(mxs,nmx,cls);
ars =  cell(mxs,nmx);
ars(:) = {''};
ars(bos) = [spl{:}];
%
%% Convert Numbers to Numeric %%
%
if nmx
    tmp = [mat{:}];
    if strcmp(fmt,'%b')
        tmp = regexprep(tmp,'^0[Bb]','');
        vec = cellfun(@(s)sum(pow2(s-'0',numel(s)-1:-1:0)),tmp);
    else
        vec = sscanf(sprintf(' %s',tmp{:}),fmt);
    end
    assert(numel(vec)==numel(tmp),'The %s format must return one value for each input number.',fmt)
else
    vec = [];
end
%
%% Debugging Array %%
%
if nmx && nargout>2
    dbg = cell(mxs,nmx);
    dbg(:) = {''};
    dbg(bon) = num2cell(vec);
    dbg = reshape(permute(cat(3,ars,dbg),[3,1,2]),[],nmx).';
    idf = [find(~all(cellfun('isempty',dbg),1),1,'last'),1];
    dbg = dbg(:,1:idf(1));
else
    dbg = {};
end
%
%% Sort Columns %%
%
if ~any(ccm) % ignorecase
    ars = lower(ars);
end
%
if nmx && any(chb) % char<num
    boe = ~cellfun('isempty',ars(bon));
    for k = reshape(find(bon),1,[])
        ars{k}(end+1) = char(65535);
    end
    [idr,idc] = find(bon);
    idn = sub2ind(size(bon),boe(:)+idr(:),idc(:));
    bon(:) = false;
    bon(idn) = true;
    arn(idn) = vec;
    bon(isnan(arn)) = ~any(nab);
    ndx = 1:nmx;
    if any(sdd) % descending
        for k = mxs:-1:1
            [~,idx] = sort(nsGroup(ars(k,ndx)),'descend');
            ndx = ndx(idx);
            [~,idx] = sort(arn(k,ndx),'descend');
            ndx = ndx(idx);
            [~,idx] = sort(bon(k,ndx),'descend');
            ndx = ndx(idx);
        end
    else % ascending
        for k = mxs:-1:1
            [~,idx] = sort(ars(k,ndx));
            ndx = ndx(idx);
            [~,idx] = sort(arn(k,ndx),'ascend');
            ndx = ndx(idx);
            [~,idx] = sort(bon(k,ndx),'ascend');
            ndx = ndx(idx);
        end
    end
else % num<char
    arn(bon) = vec;
    bon(isnan(arn)) = ~any(nab);
    if any(sdd) % descending
        [~,ndx] = sort(nsGroup(ars(mxs,:)),'descend');
        for k = mxs-1:-1:1
            [~,idx] = sort(arn(k,ndx),'descend');
            ndx = ndx(idx);
            [~,idx] = sort(bon(k,ndx),'descend');
            ndx = ndx(idx);
            [~,idx] = sort(nsGroup(ars(k,ndx)),'descend');
            ndx = ndx(idx);
        end
    else % ascending
        [~,ndx] = sort(ars(mxs,:));
        for k = mxs-1:-1:1
            [~,idx] = sort(arn(k,ndx),'ascend');
            ndx = ndx(idx);
            [~,idx] = sort(bon(k,ndx),'ascend');
            ndx = ndx(idx);
            [~,idx] = sort(ars(k,ndx));
            ndx = ndx(idx);
        end
    end
end
%
ndx  = reshape(ndx,size(X));
X = X(ndx);
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%natsort
function nsAssert(val,inp,idx,varargin)
% Throw an error if an option is overspecified.
if nnz(idx)>val
    tmp = {'Unknown input arguments',' option may only be specified once. Provided inputs'};
    error('%s:%s',[varargin{:},tmp{1+val}],sprintf('\n''%s''',inp{idx}))
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nsAssert
function grp = nsGroup(vec)
% Groups of a cell array of strings, equivalent to [~,~,grp]=unique(vec);
[vec,idx] = sort(vec);
grp = cumsum([true,~strcmp(vec(1:end-1),vec(2:end))]);
grp(idx) = grp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nsGroup