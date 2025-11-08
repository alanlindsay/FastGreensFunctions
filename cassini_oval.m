function [r, rp, rpp] = cassini_oval(t, a, b)
% CASSINI_PARAM  Parametric Cassini oval with exact derivatives.
%
%   [r, rp, rpp] = cassini_param(t, a, b)
%
%   Inputs:
%       t : parameter array (0 ≤ t ≤ 2π)
%       a : half distance between foci
%       b : Cassini constant
%
%   Outputs:
%       r   = [x; y]   positions
%       rp  = [x'; y'] first derivatives
%       rpp = [x''; y''] second derivatives
%
%   Equation in polar form:
%       r^4 - 2a^2 r^2 cos(2t) = b^4 - a^4

t = t(:)'; % row vector

h = b^4 - a^4 * sin(2*t).^2;
chi = a^2 *cos(2*t) + sqrt(h);

% --- r(t) ---
%R = sqrt(a^2 * cos(2*t) + sqrt(b^4 - a^4 * sin(2*t).^2) );  % r(t)

R = sqrt(chi);

% --- r'(t) ---

hp =  -2 * a^4 * sin(4*t); % derivative of inner sqrt term
hpp = -8 * a^4 * cos(4*t);

chip =  -2*a^2 * sin(2*t) + 0.5*hp./sqrt(h);
chipp = -4*a^2 * cos(2*t) + 0.5*hpp./sqrt(h) - 0.25*(hp.^2)./(sqrt(h).^3);

dR_dt = 0.5*chip./sqrt(chi); % derivative of r = sqrt(chi);

% --- r''(t) ---
% Derivative of dR_dt analytically

d2R_dt2 = 0.5*chipp./sqrt(chi) - 0.25*(chip.^2)./(sqrt(chi).^3);

% --- Convert to Cartesian and derivatives ---
x = R .* cos(t);
y = R .* sin(t);

xp  = dR_dt .* cos(t) - R .* sin(t);
yp  = dR_dt .* sin(t) + R .* cos(t);

xpp = d2R_dt2 .* cos(t) - 2 * dR_dt .* sin(t) - R .* cos(t);
ypp = d2R_dt2 .* sin(t) + 2 * dR_dt .* cos(t) - R .* sin(t);

r   = [x; y];
rp  = [xp; yp];
rpp = [xpp; ypp];
end