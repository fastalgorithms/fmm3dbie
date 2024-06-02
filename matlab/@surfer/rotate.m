function [objout,varargout] = rotate(obj, euls)

%%%%%
%
%       .   .   .   this routine rotates a surfer object
%          according to the three euler angles specified by euls.
%          in fact, it does so by calling the more general affine 
%          transform routine.
%
%%%%%

    p1 = euls(1);
    p2 = euls(2);
    p3 = euls(3);
    t1 = [cos(p3),sin(p3),0;-sin(p3),cos(p3),0;0,0,1];
    t2 = [1,0,0;0,cos(p2),sin(p2);0,-sin(p2),cos(p2)];
    t3 = [cos(p1),sin(p1),0;-sin(p1),cos(p1),0;0,0,1];
    objout = affine_transf(obj,t1*t2*t3,[0,0,0].');
    
end
