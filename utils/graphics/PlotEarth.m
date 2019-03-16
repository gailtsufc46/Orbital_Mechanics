function [f,globe] = PlotEarth( fig, plane )

%PlotEarth    
%
%   f = PlotEarth( fig, plane )  will create a plot of the Earth in the
%   figure "fig". In this case, f=fig. 
%   If plane=1, it will include the orbital plane.
%   
%   f = PlotEarth will plot the Earth in a new figure window. By default
%   the orbit plane is not drawn.
%   
%   f = PlotEarth( [], plane) will plot the Earth in a new figure window. 
%   If plane=1, it will include the orbital plane.
%
%   Inputs:
%     fig  
%     plane  
%
%   Outputs: 
%     f

%% Draw the Earth

if( nargin < 1 || ~ishandle(fig) )
  f = figure('color', 'k');
else
  figure(fig);
  f=fig;
end

image_file = '1024px-Land_ocean_ice_2048.jpg';
set(gca, 'NextPlot','add', 'Visible','off');

axis equal;
axis auto;
axis vis3d;
hold on;

Re = 6378.14; % equatorial radius (km)
[x, y, z] = ellipsoid(0, 0, 0, Re, Re, Re, 180);
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);

cdata = imread(image_file);
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
set(f,'color','k');

% equatorial plane
if( nargin>1 && plane )
  xc = 2*a*cos(0:.01:2*pi); yc = 2*a*sin(0:.01:2*pi);
  ep = fill3(xc,yc,xc*0,...
    'c','facealpha',.2);
end

