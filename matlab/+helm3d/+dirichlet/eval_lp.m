function p = eval(S,targinfo,zpars,sigma,varargin)
%
%  helm.dir.eval
%    Evaluates the helmholtz dirichlet layer potential at a collection 
%    of targets
%
%  Syntax
%   pot = helm.dir.eval(S,targinfo,zpars,sigma)
%   pot = helm.dir.eval(S,targinfo,zpars,sigma,opts)
%
%  Integral representation
%     pot = \alpha S_{k} [\sigma] + \beta D_{k} [\sigma]
%
%  S_{k}, D_{k}: helmholtz single and double layer potential
%  
%  k, \alpha, beta = zpars(1:3)
%    
%
%

