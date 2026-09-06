function [] = drawMars(units)

% Draws latitude and longitude lines in 3D, in specified units.

if strcmp(units, 'sm') == 1
    r = 2110.177; % Mars radius in statute miles.
elseif strcmp(units, 'km') == 1
    r = 3396; % Mars radius in km.
elseif strcmp(units, 'm') == 1
    r = 3396*1000; % Mars radius in m.
elseif strcmp(units, 'nm') == 1
    r = 1833.693; % Mars radius in nautical miles.
elseif strcmp(units, 'ft') == 1
    r = 2110.177*5280; % Mars radius in feet.
end

spacing = 10; % 10° spacing.

%% Draw Longitude Lines

[long1, lat1] = meshgrid(-180:spacing:180, linspace(-90, 90, 300));

[x1, y1, z1] = sph2cart(deg2rad(long1), deg2rad(lat1), r);

plot3(x1, y1, z1, 'r', 'LineWidth', 1.0)

hold on

%% Draw Latitude Lines

[lat2, long2] = meshgrid(-90:spacing:90, linspace(-180, 180, 300));

[x2, y2, z2] = sph2cart(deg2rad(long2), deg2rad(lat2), r);

plot3(x2, y2, z2, 'r', 'LineWidth', 1.0)

axis equal

end