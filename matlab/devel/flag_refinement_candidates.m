function [refinelist,errqlist] = flag_refinement_candidates( ...
    S,fun,funlp,tol_amr,lptype,eta_dens,norderhead,functional)

 tail_errs = testfuns.surf_fun_tail_val_error(S, fun, lptype, norderhead,...
     functional);

[refinelist,errqlist] = find_patch_refine_list( ...
    S,tail_errs,funlp,tol_amr,[],lptype,eta_dens);

refinelist = refinelist(:);

end
