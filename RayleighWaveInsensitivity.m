%% Copyright
% This work is supplied under the common development and distribution license
% (CDDL) Version 1.0
% 
%% Creation details
% Author: Dr. James Andrews
% Employer: University of Birmingham
% Date: 01/04/2019
%
%% Disclaimer
% The author and the employer is not liable for any damages arising in contract,
% tort or otherwise from the use of or inability to use this script 
% or any material contained in it, or from any action or decision taken as a 
% result of using this script.
%
% This script comprises the author's view; it does not constitute legal or 
% professional advice. f you need specific advice, please seek a professional 
% who is licensed or knowledgeable in that area. 
%
%% Moisture insensitivity
% Fig 4

% Material parameters
h = [10^(-2) 3 10^(-2) 100]*10^(-3);  % Thickness of each layer (m)
Evals = [10^2 5*10^2 10^3 5*10^3 10^4];
for j=1:length(Evals)
  E = [Evals(j) 32 10^7 221*10^3]*10^3;      % Young's modulus (Pa)
  nu = [0.3 0.4 0.31 0.494];            % Poisson's ratio (-)
  rho = [1000 1000 1000 1800];          % Density (kg/m^3)
  % Derived parameters
  cp = sqrt((E.*nu./((1 + nu).*(1 - 2*nu)) + 2*E./(2 + 2*nu))./rho);  % P-wave speed
  cs = sqrt(E./(2*rho.*(1+nu)));              % S-wave speed

  % Pacinian corspucle dimensions
  a = 0.1 * 10^(-3);  % Semi-minor axis (m)
  b = 0.12 * 10^(-3); % Semi-major axis (m)

  % Frequency range
  f = 10^2*[1:1:10];

  % Second Lame constant
  mu = E.*nu./((1+nu).*(1-2*nu));

  cr = fdispersion_rayleigh(cs./cp, 0, rho, h, mu)*cs(4);
end;