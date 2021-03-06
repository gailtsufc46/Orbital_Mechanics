function S = stumpffS(z)
% Written by Garrett Ailts
%
% Usage: S = @(z) stumpffS(z); (Then use the function as a coeff dependent
%        on the universal variable z
% 
% Description: Stumpff function for the S coeffciient in the universal
% Kepler equation
%
% Inputs: z - universal variable
%
% Outputs: S - Stumpff S coefficient

%% Calculate S(z) If z is a Vector
if length(z)>1
   S = zeros(1,length(z));
   gidx = find(z>0);
   lidx = find(z<0);
   eqidx = z==0; % faster than find
   S(gidx) = (sqrt(z(gidx))-sin(sqrt(z(gidx))))./sqrt(z(gidx)).^3;
   S(lidx) = (sinh(sqrt(-z(lidx)))-sqrt(-z(lidx)))./sqrt(-z(lidx)).^3;
   S(eqidx) = 1/6;
   
%% Calculate S(z) If z is a Scalar  
else
    if z>0
        S = (sqrt(z)-sin(sqrt(z)))/sqrt(z).^3;
    elseif z<0
        S = (sinh(sqrt(-z))-sqrt(-z))/sqrt(-z).^3;
    else
        S = 1/6;
    end
end
    
    