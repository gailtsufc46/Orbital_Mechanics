function [thMin, thMax] = TrueAnomalyRange( e, tol )

if( e<1 )
  thMin = -pi;
  thMax = pi;
elseif e==1
  thMin = -pi+tol;
  thMax = pi-tol;
elseif e>1 
  thMin = -acos(-1/e)+tol;
  thMax = acos(-1/e)-tol;
end
