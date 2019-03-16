function m = EulerRot( axis, angle )

% Compute a rotation matrix for Euler angle rotation about an axis
%
% Inputs:
%   axis    1 for x-axis, 2 for y-axis, 3 for z-axis
%   angle   The angle of rotation about the axis (rad)
%

c = cos(angle);
s = sin(angle);

if( axis == 1 )
  m = [1 0 0;0 c s; 0 -s c];
elseif( axis == 2 )
  m = [c 0 -s; 0 1 0; s 0 c];
elseif( axis == 3 )
  m = [c s 0; -s c 0; 0 0 1];
else
  error('Axis must be 1, 2, or 3.');
end
