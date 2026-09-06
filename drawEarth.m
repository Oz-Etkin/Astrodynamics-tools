function [] = drawEarth(units, color)

% Draws latitude and longitude lines in 3D, in specified units and color.

if strcmp(units, 'sm') == 1
    r = 3963.105; % Earth radius in statute miles.
elseif strcmp(units, 'km') == 1
    r = 6378; % Earth radius in km.
elseif strcmp(units, 'm') == 1
    r = 6378*1000; % Earth radius in m.
elseif strcmp(units, 'nm') == 1
    r = 3443.844; % Earth radius in nautical miles.
elseif strcmp(units, 'ft') == 1
    r = 3963.105*5280; % Earth radius in feet.
end

spacing = 10; % 10° spacing.

%% Draw Longitude Lines

[long1, lat1] = meshgrid(-180:spacing:180, linspace(-90, 90, 300));

[x1, y1, z1] = sph2cart(deg2rad(long1), deg2rad(lat1), r);

plot3(x1, y1, z1, color, 'LineWidth', 1.0)

hold on

%% Draw Latitude Lines

[lat2, long2] = meshgrid(-90:spacing:90, linspace(-180, 180, 300));

[x2, y2, z2] = sph2cart(deg2rad(long2), deg2rad(lat2), r);

plot3(x2, y2, z2, color, 'LineWidth', 1.0)

axis equal

end