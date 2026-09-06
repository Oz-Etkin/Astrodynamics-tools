function [a, e, I, AoP, LAN, f] = state2COEs(r_vec, v_vec)

% Given position and velocity vectors, compute the classical orbital
% elements a, e, I, AoP, LAN, and f (all angles in degrees) using the 
% routine described in Space Flight Dynamics by Craig A. Kluever p. 60-63.
% -- MATLAB script by Oz Etkin

mu = 398600.4418;

%%%%%%%%%%%%%%%%%%%
%% UNPACK STATES %%
%%%%%%%%%%%%%%%%%%%

r = norm(r_vec); %Position vector

v = norm(v_vec); %Velocity vector

epsilon = v^2/2 - mu/r; %Orbital energy

h_vec = cross(r_vec, v_vec); %Angular momentum vector

h = norm(h_vec); %Angular momentum scalar

%%%%%%%%%%%%%%%%%%%%%%%
%% ECI BASIS VECTORS %%
%%%%%%%%%%%%%%%%%%%%%%%

Ihat = [1 0 0]';
Jhat = [0 1 0]';
Khat = [0 0 1]';

nhat = cross(Khat, h_vec);

n = norm(nhat);

%%%%%%%%%%%%%%%%%%%%%
%% SEMI-MAJOR AXIS %%
%%%%%%%%%%%%%%%%%%%%%

a = -mu/(2*epsilon);

%%%%%%%%%%%%%%%%%%
%% ECCENTRICITY %%
%%%%%%%%%%%%%%%%%%

e_vec = 1/mu*((v^2 - mu/r)*r_vec - dot(r_vec, v_vec)*v_vec);

e = norm(e_vec);

%%%%%%%%%%%%%%%%%
%% INCLINATION %%
%%%%%%%%%%%%%%%%%

I = acosd(dot(h_vec, Khat)/h);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% ARGUMENT OF PERIGEE %%
%%%%%%%%%%%%%%%%%%%%%%%%%

AoP = acosd(dot(nhat, e_vec)/(n*e));
if e_vec(3) < 0
    AoP = 360 - AoP;
else 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LONGITUDE OF ASCENDING NODE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cLAN = dot(Ihat, nhat)/n;
sLAN = dot(Jhat, nhat)/n;

LAN  = atan2d(sLAN, cLAN);


%%%%%%%%%%%%%%%%%%
%% TRUE ANOMALY %%
%%%%%%%%%%%%%%%%%%

f = acosd(dot(e_vec, r_vec)/(e*r));

if dot(r_vec, v_vec) < 0
    f = 360 - f;
else
end

end