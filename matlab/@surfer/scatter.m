function scatter(obj,s,c,varargin)
% SCATTER  Plots values of a surfer object, with colors specified by data.
%
% scatter(obj,s,c,...) plots surfer object obj nodes with points of size s and
%  optional color values given by vector c.
%  Other args are passed to scatter3.
%
% See also: PLOT_NODES

if nargin < 2, s = 20; end

% Handle arrays of surfers: scatter each one, splitting any color vector.
if numel(obj) > 1
    ih = ishold();
    hold on
    ipt = 1;
    ntotal = sum([obj.npts]);
    for k = 1:numel(obj)
        nk = obj(k).npts;
        if nargin < 3
            scatter(obj(k), s, varargin{:});
        else
            ck = c;
            if isnumeric(c) && length(c(:)) == ntotal
                ck = c(ipt:ipt+nk-1);
            end
            scatter(obj(k), s, ck, varargin{:});
        end
        ipt = ipt + nk;
    end
    if ~ih, hold off; end
    return
end

ifhold = ishold();

x = obj.r(1,:);
y = obj.r(2,:);
z = obj.r(3,:);

if (nargin <2)
    s = 20;
    c = z;
elseif (nargin <3)
    c = z;
end

if (isempty((varargin)))
    scatter3(x,y,z,s,c,'filled');
else
    scatter3(x,y,z,s,c,varargin{:});
end

hold off

if ifhold
    hold on
end
