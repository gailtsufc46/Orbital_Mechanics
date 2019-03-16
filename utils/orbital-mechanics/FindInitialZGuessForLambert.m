function z0 = FindInitialZGuessForLambert( yf, Ff )

zmin = -1e3;
zmax = 1e3;

% first find any y(z)=0 crossings, if they exists
zz = linspace(zmin,zmax,1e4);
yz = yf(zz);
y1 = yz(1:end-1);
y2 = yz(2:end);
kyz = find(sign(y1).*sign(y2)<0);

if( isempty(kyz) )
  %disp('No y(z)=0 crossing found...')
  zz = linspace(zmin,zmax,1e4);
  z0 = FindApproxFZero( Ff, zz );
  
else
  %disp(sprintf('%d y(z)=0 crossings found...',length(kyz)))
  z0 = [];
  for i=1:length(kyz)
    
    zy0 = fzero(yf,kyz(i));
    zz = linspace(zy0,zy0+zmax,1e4);
    z0i = FindApproxFZero( Ff, zz );
    if( ~isempty(z0i) )
      z0 = [z0, z0i];
    end
    
  end
end

%disp(z0)

function z0 = FindApproxFZero( Ff, zz )

  f = Ff(zz);                  % F(z) where y(z) > 0
  
  % find approx. F(z)=0 crossing
  f1 = f(1:end-1);
  f2 = f(2:end);
  ks = find(sign(f1).*sign(f2)<0);
  z0 = zz(ks);


