function statedot = EarthSunCR3BP(t, state)

% Computes satellite trajectory in Earth-Sun CR3BP rotating frame when
% integrated numerically using ode45 or similar. Also computes state
% transition matrix.
% -- MATLAB function by Oz Etkin.

x = state(1);
y = state(2);
z = state(3);
xd = state(4);
yd = state(5);
zd = state(6);

Phi = state(7:42);

Phi = reshape(Phi, [6 6]);

m1 = 1.989E+30;
m2 = 6.0748E+24;

mu = m2/(m1 + m2);

x1 = -mu;
x2 = 1 - mu;

D = norm([x y z]' - [x1 0 0]');
R = norm([x y z]' - [x2 0 0]');

xdd = 2*yd + x - (1 - mu)*(x + mu)/D^3 - mu*(x + mu - 1)/R^3;

ydd = -2*xd + y - (1 - mu)*y/D^3 - mu*y/R^3;

zdd = -(1 - mu)*z/D^3 - mu*z/R^3;
% A  =[                                                                                                                                                                                                                                                                                                                0,                                                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                                               0,  1, 0, 0;
%                                                                                                                                                                                                                                                                                                                 0,                                                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                                               0,  0, 1, 0;
%                                                                                                                                                                                                                                                                                                                 0,                                                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                                               0,  0, 0, 1;
% (mu - 1)/(abs(mu + x)^2 + abs(y)^2 + abs(z)^2)^(3/2) - mu/(abs(mu + x - 1)^2 + abs(y)^2 + abs(z)^2)^(3/2) - (3*abs(mu + x)*sign(mu + x)*(mu + x)*(mu - 1))/(abs(mu + x)^2 + abs(y)^2 + abs(z)^2)^(5/2) + (3*mu*abs(mu + x - 1)*sign(mu + x - 1)*(mu + x - 1))/(abs(mu + x - 1)^2 + abs(y)^2 + abs(z)^2)^(5/2) + 1,                                                                                               (3*mu*abs(y)*sign(y)*(mu + x - 1))/(abs(mu + x - 1)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (3*abs(y)*sign(y)*(mu + x)*(mu - 1))/(abs(mu + x)^2 + abs(y)^2 + abs(z)^2)^(5/2),                                                                                           (3*mu*abs(z)*sign(z)*(mu + x - 1))/(abs(mu + x - 1)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (3*abs(z)*sign(z)*(mu + x)*(mu - 1))/(abs(mu + x)^2 + abs(y)^2 + abs(z)^2)^(5/2),  0, 2, 0;
%                                                                                                                                   (3*mu*y*abs(mu + x - 1)*sign(mu + x - 1))/(abs(mu + x - 1)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (3*y*abs(mu + x)*sign(mu + x)*(mu - 1))/(abs(mu + x)^2 + abs(y)^2 + abs(z)^2)^(5/2), (mu - 1)/(abs(mu + x)^2 + abs(y)^2 + abs(z)^2)^(3/2) - mu/(abs(mu + x - 1)^2 + abs(y)^2 + abs(z)^2)^(3/2) + (3*mu*y*abs(y)*sign(y))/(abs(mu + x - 1)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (3*y*abs(y)*sign(y)*(mu - 1))/(abs(mu + x)^2 + abs(y)^2 + abs(z)^2)^(5/2) + 1,                                                                                                             (3*mu*y*abs(z)*sign(z))/(abs(mu + x - 1)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (3*y*abs(z)*sign(z)*(mu - 1))/(abs(mu + x)^2 + abs(y)^2 + abs(z)^2)^(5/2), -2, 0, 0;
%                                                                                                                                  (3*mu*z*abs(mu + x - 1)*sign(mu + x - 1))/(abs(mu + x - 1)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (3*z*abs(mu + x)*sign(mu + x)*(mu - 1))/(abs(mu + x)^2 + abs(y)^2 + abs(z)^2)^(5/2),                                                                                                                 (3*mu*z*abs(y)*sign(y))/(abs(mu + x - 1)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (3*z*abs(y)*sign(y)*(mu - 1))/(abs(mu + x)^2 + abs(y)^2 + abs(z)^2)^(5/2), (mu - 1)/(abs(mu + x)^2 + abs(y)^2 + abs(z)^2)^(3/2) - mu/(abs(mu + x - 1)^2 + abs(y)^2 + abs(z)^2)^(3/2) + (3*mu*z*abs(z)*sign(z))/(abs(mu + x - 1)^2 + abs(y)^2 + abs(z)^2)^(5/2) - (3*z*abs(z)*sign(z)*(mu - 1))/(abs(mu + x)^2 + abs(y)^2 + abs(z)^2)^(5/2),  0, 0, 0];
% 

Uxx = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (3*mu*(2*mu + 2*x - 2)^2)/(4*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(2*mu + 2*x)^2*(mu - 1))/(4*((mu + x)^2 + y^2 + z^2)^(5/2)) + 1;
Uxy = (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2));
Uxz = (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);
Uyy = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) + 1;
Uyz = (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);
Uzz = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*z^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*z^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2);

K = -[Uxx Uxy Uxz;
      Uxy Uyy Uyz;
      Uxz Uyz Uzz];
C = -[0 2 0;
     -2 0 0;
      0 0 0];

A = [zeros(3), eye(3);
    -K -C];

Phidot = A*Phi;

Phidot = reshape(Phidot, [1 36]);

statedot = [xd yd zd xdd ydd zdd Phidot]';


end