function S = stumpS( z )

if( length(z)>1 )
  
  % vectorized
  kp = find(z>0);
  kn = find(z<0);
  ke = z==0;
  sz = zeros(size(z));
  S = zeros(size(z));
  sz(kp) = sqrt(z(kp));
  sz(kn) = sqrt(-z(kn));
  sz(ke) = 0;
  S(kp) = (sz(kp)-sin(sz(kp)))./(sz(kp).^3);
  S(kn) = (sinh(sz(kn))-sz(kn))./(sz(kn).^3);
  S(ke) = 1/6;

else
  
  if( z>0 )
    sz = sqrt(z);
    S = (sz-sin(sz))/(sz^3);
  elseif( z<0 )
    sz = sqrt(-z);
    S = (sinh(sz)-sz)/(sz^3);
  else
    S = 1/6;
  end

end