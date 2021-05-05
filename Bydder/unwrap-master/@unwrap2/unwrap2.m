function output = unwrap2(varargin)
%
%Phase unwrapping in 2D using mex wrapper to
%code at https://github.com/geggo/phase-unwrap
%
%output = unwrap2(input)
%output = unwrap2(input,mask)
%
%Inputs
%-input: wrapped phase between -pi and +pi
%-mask: binary quality mask (1=keep 0=reject)
%
%Outputs
%-output: unwrapped phase

mexbin = [mfilename('fullpath') '.' mexext];

if ~exist(mexbin,'file')
    mexcpp = [mfilename('fullpath') '.cpp'];
    mex('-v',mexcpp,'-output',mexbin);
end

output = unwrap2(varargin{:});
