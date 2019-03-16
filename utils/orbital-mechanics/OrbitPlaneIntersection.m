function u = OrbitPlaneIntersection( el1, el2 )

% Compute the line of intersection (unit vector) between 2 orbital planes.

if( nargin<1 )
  el1 = [8000, 0, 0, 0, 0, 0];
  el2 = [9000, pi/4, pi/2, pi/3, 0.1, 0];
end

mu = 1;
[r1,v1] = RVFromCOE(el1,mu);
h1 = cross(r1,v1);

[r2,v2] = RVFromCOE(el2,mu);
h2 = cross(r2,v2);

v = cross(h1,h2);
if( norm(v)>eps )
  u = unit( cross(h1, h2) );
else
  error('The planes are co-planar.');
end


function uu = unit(vv)
uu = vv/sqrt(vv'*vv);