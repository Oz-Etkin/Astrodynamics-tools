function [L, G, H, l, g, h] = COE2Delaunay(a, e, I, AoP, LAN, M, mu)

% Converts classical orbital elements to Delaunay.

l = M;
g = AoP;
h = LAN;
L = sqrt(mu*a);
G = L*sqrt(1 - e^2);
H = G*cosd(I);

end