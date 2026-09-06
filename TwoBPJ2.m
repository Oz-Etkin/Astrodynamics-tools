function statedot = TwoBPJ2(t, state)

% Computes satellite trajectory with J2 perturbation when integrated
% numerically using ode45 or similar.
% -- MATLAB function by Oz Etkin.

Re = 6378.137;
mu = 398600.4418;
J2 = 0.001083;

x = state(1);
y = state(2);
z = state(3);
xdot = state(4);
ydot = state(5);
zdot = state(6);

r_vec = [x y z]';
r = norm(r_vec);

aJ2 = -3/2*J2*(mu/r^2)*(Re/r)^2*[(1 - 5*(z/r)^2)*x/r (1 - 5*(z/r)^2)*y/r (3 - 5*(z/r)^2)*z/r]';

a = -mu/r^3*r_vec + aJ2;

statedot = [xdot ydot zdot a(1) a(2) a(3)]';

end
