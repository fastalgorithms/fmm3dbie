function A = gsmatgen(S, zpars, eps)
% Superseded by surfwave.capillary.matgen(S,'gs',zpars,eps).
% Retained for backwards compatibility only.
A = surfwave.capillary.matgen(S, 'gs', zpars, eps).';
end
