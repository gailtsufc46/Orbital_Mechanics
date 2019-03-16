function [alpha,beta] = LambertAlphaBeta( a, c, r1r2 )

% Compute the angles alpha and beta for Lambert's theorem
%
%   Inputs: 
%     a     Semi major axis
%     c     Chord (distance between r1 and r2)
%     r1r2  Sum of r1 and r2
% 
% 	Outputs:
%     alpha   
%     beta    


s = (r1r2+c)/2;

alpha = 2*asin( [1 -1]*sqrt( s/(2*a) ) );
beta  = 2*asin( [1 -1]*sqrt( (s-c)/(2*a) ) );
