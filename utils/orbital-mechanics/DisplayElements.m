function DisplayElements( el )

%  Display an orbital elements vector

a = el(1);
inc = el(2)*180/pi;
W = el(3)*180/pi;
w = el(4)*180/pi;
e = el(5);
th = el(6)*180/pi;
fprintf(1,' a  = %2.1f km\n i  = %2.1f deg\n W  = %2.1f deg\n w  = %2.1f deg\n e  = %2.5f\n TA = %2.1f deg\n\n',...
  a,inc,W,w,e,th);
