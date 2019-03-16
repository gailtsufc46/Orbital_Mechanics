function C = stumpC( z )

if( length(z) > 1 )
  
  % vectorized
  kp = find(z>0);
  kn = find(z<0);
  ke = z==0;
  C = zeros(size(z));
  C(kp) = (1-cos(sqrt(z(kp))))./z(kp);
  C(kn) = (cosh(sqrt(-z(kn)))-1)./(-z(kn));
  C(ke) = 0.5;

else
  
  if( z>0 )
    C = (1-cos(sqrt(z)))./z;
  elseif( z<0 )
    C = (cosh(sqrt(-z))-1)/(-z);
  else
    C = 0.5;
  end
  
end