function f = mtimes(f, g)
% * Matrix multiplication for kernel3d class
%
% Currently only supports scalars: returns c*F or F*c for kernel F
% and scalar c.

f = times(f, g);
end
