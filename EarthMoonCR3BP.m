function statedot = EarthMoonCR3BP(t, state)

% Computes satellite trajectory in Earth-Moon CR3BP rotating frame when
% integrated numerically using ode45 or similar.
% -- MATLAB function by Oz Etkin.

x = state(1);
y = state(2);
z = state(3);
xd = state(4);
yd = state(5);
zd = state(6);

m1 = 6.0748E+24;
m2 = 7.34767E+22;

mu = m2/(m1 + m2);

x1 = -mu;
x2 = 1 - mu;

D = norm([x y z]' - [x1 0 0]');
R = norm([x y z]' - [x2 0 0]');

xdd = 2*yd + x - (1 - mu)*(x + mu)/D^3 - mu*(x + mu - 1)/R^3;

ydd = -2*xd + y - (1 - mu)*y/D^3 - mu*y/R^3;

zdd = -(1 - mu)*z/D^3 - mu*z/R^3;

statedot = [xd yd zd xdd ydd zdd]';

end