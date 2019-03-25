function dF = DFofZ(r1mag, r2mag, A, z)
% Written by Garrett Ailts
%
% Usage: dF = @(z) DFofZ(r1mag, r2mag, A, z)
%
% Inputs: r1mag - magnitude of r1
%         r2mag - magnitude of r2
%             A - constant of universal variable equation
%             z - universal variable
% 
% Outputs: dF - the derivative of the Universal Kepler's Equation
%               evaluated at z

%% Calculate df(z)
% if z is close to zero
if abs(z)<eps % eps - spacing between precision numbers (2^-52)
    y = r1m+r2m-sqrt(2)*A;
    dF = sqrt(2)/y^(3/2)/40+(A*(sqrt(y)+Asqrt(1/2/y))/8);
else
    S = stumpffS(z);
    C = stumpffC(z);
    y = r1mag+r2mag+A*(z*S-1)/sqrt(C);
    dF = (y/C)^(3/2)*((C-(3*S/2/C)/2/z)+3*S^2/4/C)+ ...
         A*(3*S*sqrt(y)/C+A*sqrt(C/y))/8;
end
    
    
