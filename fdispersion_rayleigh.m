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
%% Function dispersion relation

function fdispersion_rayleigh=fdispersion_rayleigh(alpha, beta, rho, h, mu, A)
  %% alpha: cs/cp 
  %% beta:  flag 0 -> terrestrial animal, flag 1 -> (semi-)aquatic animal
  %% rho:   density vector consisting of density of each layer (kg/m^3)
  %% h:     thickness vector consisting of thickness of each layer (m)
  %% mu:    second Lame constant for each layer
  %% Returns c/cs for the three layer Rayleigh wave
  elayers = 0;
  for i=1:3
    elayers = h(i)*(1-sqrt(1-alpha(i+1)^2)/sqrt(1-alpha(i+1)^2))*(rho(i)/rho(4)-mu(i)/mu(4))+elayers;
  end
  % Newton-Raphson
  x = 6.59092450;
  tol = 10^(-6);
  err = 1;
  while err>tol    
    x0 = x;
    F = (2-x^2)^4-16*(1-alpha(4)^2*x^2)*(1-x^2)+0.522*beta*2*pi*rho(1)*h(1)*x/rho(4)+4*x^2*(1-x^2)*elayers;
    dFdx = -8*(-x^2 + 2)^3*x + 32*alpha(4)^2*x*(-x^2 + 1) + 32*(1 - alpha(4)^2*x^2)*x + 1.044*beta*pi*rho(1)*h(1)/rho(4) + 8*x*(-x^2 + 1)*elayers - 8*x^3*elayers;
    x = x-F/dFdx;
    err = abs((x-x0)/x);
  end
  fdispersion_rayleigh = x;
  
endfunction
