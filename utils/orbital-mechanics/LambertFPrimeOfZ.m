function dF = LambertFPrimeOfZ( r1m, r2m, dTheta, z )

% Compute the derivative of F(z) 
%
%   Given two position vectors, r1 and r2, and the time of flight between
%   them, TOF, find the corresponding Keplerian orbit. 
%   
%   Inputs: 
%     r1m     Magnitude of position vector r1  
%     r2m     Magnitude of position vector r2
%     dTheta  Angle between vector r1 and r2
%     z       Universal variable "z"
% 
% 	Outputs:
%     dF      Derivative of F(z) at z

A = sin(dTheta)*sqrt(r1m*r2m/(1-cos(dTheta)));

if( abs(z)<eps )
  y0 = r1m+r2m-sqrt(2)*A;
  dF = sqrt(2)/40*y0^1.5 + A/8*( sqrt(y0) + A*sqrt(1/2/y0));
else
  S = stumpS(z);
  C = stumpC(z);
  y = r1m+r2m+A*( z*S-1 )/sqrt(C);  
  dF = (y/C)^1.5 * ( 1/(2*z) * (C-3*S/(2*C))+3*S^2 / (4*C) ) + ...
    (A/8)*( 3*S/C*sqrt(y) + A*sqrt(C/y) );
end


