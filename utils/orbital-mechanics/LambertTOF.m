function t = LambertTOF( a, c, r1r2, mu )

% Compute the time of flight using Lambert's theorem
%
%   Inputs: 
%     a     Semi major axis
%     c     Chord (distance between r1 and r2)
%     r1r2  Sum of r1 and r2
%     mu    Gravitational constant
% 
% 	Outputs:
%     t     Time of flight between r1 and r2    


% use MU for earth by default if not given
if( nargin<4 )
  mu = 396800.44;
end

[alpha,beta] = LambertAlphaBeta( a, c, r1r2 );

n = sqrt(mu/a^3);

t = zeros(1,4);
t(1) = (1/n) * ( alpha(1) - beta(1) - (sin(alpha(1))-sin(beta(1))) );
t(2) = (1/n) * ( alpha(1) - beta(2) - (sin(alpha(1))-sin(beta(2))) );
t(3) = (1/n) * ( alpha(2) - beta(1) - (sin(alpha(2))-sin(beta(1))) );
t(4) = (1/n) * ( alpha(2) - beta(2) - (sin(alpha(2))-sin(beta(2))) );

t( t<0 ) = t( t<0 ) + 2*pi/n;
