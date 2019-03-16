function arc = GenArc( b, bp, n )

% Generate an arc between two vectors. It will have unit length 1.
% Call "GenArc" with no inputs for an example.
%
% Inputs:
%   b     First vector, must be size:  3 x 1
%   bp    Second vector, must be size: 3 x 1
%   n     Number of points for the arc. Optional. Default is 100.
%
% Outputs:
%   arc   Array of vectors connecting b and bp in an arc. Size: 3 x n
%

if( nargin<1 )
  b = unit(randn(3,1))*2;
  bp = unit(randn(3,1))*2;
  arc = GenArc( b, bp );
  figure, quiver3(0,0,0,b(1),b(2),b(3)), hold on
  quiver3(0,0,0,bp(1),bp(2),bp(3)), grid on
  plot3(arc(1,:),arc(2,:),arc(3,:),'k.')
  axis equal, rotate3d on
  arc=[];
  return
end


if( nargin<3 )
  n = 100;
end

if( norm(b)<eps || norm(bp)<eps )
  arc = zeros(3,1);
  return
end

% force b and bp to be unit column vectors.
b = b(:)/norm(b);
bp = bp(:)/norm(bp);

a = cross(b,bp);
if( norm(a)<1e-10 && b'*bp<0 )
  b = b + randn(3,1)*1e-8;
%   arc = zeros(3,1);
%   return;
  a = cross(b,bp);
end

a = unit(a);
c = unit(cross(a,b));

rot = [a';b';c'];

% the angles from b to bp range from zero to the max angle "ang"
dotProd = b'*bp;
% prevent acos from being used if the abs of dot product exceeds 1 (due to
% rounding errors)
if( dotProd>1 )   
  ang = 0;
elseif( dotProd < -1 )
  ang = pi;
else
  ang = acos(b'*bp);
end
thx = linspace(0,ang,n);

% define the arc of points in the abc frame
arc_abc = [zeros(1,n); cos(thx); sin(thx)];

% rotate back to the original frame
arc = rot'*arc_abc;

if( any(imag(arc(:))) )
  1;
end

function u=unit(v)
u=v/sqrt(v'*v);

