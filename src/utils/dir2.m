function listing = dir2(varargin)
% Lists directory contents without the annoying '.' and '..' folders.

if nargin == 0
    name = '.';
elseif nargin == 1
    name = varargin{1};
else
    error('Too many input arguments.')
end

listing = dir(name);

listing = listing(~ismember({listing.name}, {'.', '..', '.DS_Store'}))
