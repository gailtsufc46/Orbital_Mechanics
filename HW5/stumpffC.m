function C = stumpffC(z)
% Written by Garrett Ailts
%
% Usage: C = @(z) stumpffS(z); (Then use the function as a coeff dependent
%        on the universal variable z
% 
% Description: Stumpff function for the C coeffciient in the universal
% Kepler equation
%
% Inputs: z - universal variable
%
% Outputs: C - C coefficient

%% Calculate C(z) If z is a Vector
if length(z)>1
   C = zeros(1,length(z));
   gidx = find(z>0);
   lidx = find(z<0);
   eqidx = z==0; % faster than find
   C(gidx) = (1-cos(sqrt(z(gidx))))./z(gidx);
   C(lidx) = (cosh(sqrt(-z(lidx)))-1)./-z(lidx);
   C(eqidx) = 1/2;
   
%% Calculate C(z) If z is a Scalar    
else
    if z>0
        C = (1-cos(sqrt(z)))/z;
    elseif z<0
        C = (cosh(sqrt(-z))-1)/-z;
    else
        C = 1/2;
    end
end
end