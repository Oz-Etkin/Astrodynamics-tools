function [] = drawOrbit(a, e, I, AoP, LAN)

% Draws and orbit in ECI frame given semi-major axis, eccentricity,
% inclination, argument of periapsis, and longitude of ascending node.
% Angles in degrees, distances in km.
% Dependent on COEs2state and drawEarth.
% -- MATLAB function by Oz Etkin

n = 1;
i = zeros(360, 1);
j = zeros(360, 1);
k = zeros(360, 1);
for f = 0:360
    [r_vec, v_vec] = COEs2state(a, e, I, AoP, LAN, f);
    i(n) = r_vec(1);
    j(n) = r_vec(2);
    k(n) = r_vec(3);
    n = n + 1;
end


plot3(i, j, k, 'y', 'LineWidth', 2.0)
xlabel('I (km)');
ylabel('J (km)');
zlabel('K (km)');
axis equal
hold on
drawEarth('km', 'w')
theme dark

end