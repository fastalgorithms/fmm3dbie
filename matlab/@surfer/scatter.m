function scatter(obj,s,c,varargin)
% plots values of a surfer object, with colors specified by data

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
