function obj = tangent_kern(kernin, ids, jds)
%KERNEL3D.TANGENT_KERN   Wrap a kernel so vector blocks are expressed in
%   the local tangent basis (du, dv).
%
%   OBJ = KERNEL3D.TANGENT_KERN(KERNIN, IDS, JDS)
%
%   IDS - 3 x ki matrix; each column gives the indices of a source vector
%         block.  Unlisted indices are scalar and pass through unchanged.
%   JDS - 3 x ko matrix; same for target.
%
%   Equivalent to:
%       f_right(s) = tangent expansion pages (q_cart x p_tan x ns)
%       f_left(t)  = tangent contraction pages (p_tan x m_cart x nt)
%       OBJ = f_left .* KERNIN .* f_right
%
%   OBJ.opdims = [tgt_opdim - ko, src_opdim - ki].
%
%   See also TANGENT_BLOCK, KERNEL3D/TIMES

ki = size(ids, 2);
ko = size(jds, 2);

src_opdim_in  = kernin.opdims(2);
tgt_opdim_in  = kernin.opdims(1);
src_opdim_out = src_opdim_in - ki;
tgt_opdim_out = tgt_opdim_in - ko;

src_scalar_rows = setdiff((1:src_opdim_in)', ids(:));
tgt_scalar_rows = setdiff((1:tgt_opdim_in)', jds(:));

% f_right(s): expansion pages (src_opdim_in x src_opdim_out x ns)
    function F = f_right(s)
        F = tangent_block(s, ids, src_scalar_rows, src_opdim_out, false);
    end

% f_left(t): contraction pages (tgt_opdim_out x tgt_opdim_in x nt)
    function F = f_left(t)
        F = tangent_block(t, jds, tgt_scalar_rows, tgt_opdim_out, true);
    end

% Right multiply by f_right  ->  kernin * f_right
% then left  multiply by f_left   ->  f_left * (kernin * f_right)
obj = times(@f_left, times(kernin, @f_right));

% Restore metadata from the inner kernel
obj.name         = kernin.name;
obj.type         = [kernin.type, '_tangent'];
obj.zk           = kernin.zk;
obj.ifcomplex    = kernin.ifcomplex;
obj.kernel_order = kernin.kernel_order;
obj.src_fields   = union(kernin.src_fields,  {'du', 'dv'});
obj.targ_fields  = union(kernin.targ_fields, {'du', 'dv'});
if isa(kernin.get_overs_orders, 'function_handle')
    obj.get_overs_orders = kernin.get_overs_orders;
end

end
