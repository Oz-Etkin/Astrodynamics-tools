function [r_vec, v_vec] = COEs2state(a, e, I, AoP, LAN, f)

% Given the classical orbital elements, find the state vector using the 
% routine from Space Flight Dynamics by Craig A. Kluever p. 66-73.
% All angles in degrees, distances in km.
% -- MATLAB function by Oz Etkin

mu = 398600.4418; % Gravitational parameter of Earth, km^3/s^2

p = a*(1 - e^2); %Semi-latus rectum

r = p/(1 + e*cosd(f)); %Current r

h = sqrt(mu*p); %Angular momentum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STATE VECTOR IN PERIFOCAL COORDINATES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_pqw = r*cosd(f);

y_pqw = r*sind(f);

rPQW = [x_pqw y_pqw 0]';

xdot_pqw = -mu/h*sind(f);

ydot_pqw = mu/h*(e + cosd(f));

vPQW = [xdot_pqw ydot_pqw 0]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONVERT FROM PERIFOCAL TO ECI %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C_AoP = [cosd(-AoP) -sind(-AoP) 0;
         sind(-AoP)  cosd(-AoP) 0;
         0           0         1];

A_I = [1 0         0;
       0 cosd(-I) -sind(-I);
       0 sind(-I)  cosd(-I)];

C_LAN = [cosd(-LAN) -sind(-LAN) 0;
         sind(-LAN)  cosd(-LAN) 0;
         0           0          1];

R = (C_AoP*A_I*C_LAN)';

r_vec = R*rPQW;

v_vec = R*vPQW;

end