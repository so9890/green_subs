function f = calibration_production(x, MOM, list, paramss, polhelp)
% function determines pn, pg, and delta in equilibrium! 
% percieves Af, Ag, An in 2015-2019 as parameters but from here back
% out Af0, Ag0, An0

% 1) read in parameters and variables
% parameters
eppse=paramss(list.paramsdir=='eppse');
eppsy=paramss(list.paramsdir=='eppsy');
tauf =polhelp(list.poldir=='tauf');

% variables
pn = exp(x(list.prod=='pn'));
pg = exp(x(list.prod=='pg'));
deltay = 1/(1+exp(x(list.prod=='deltay')));


% 2) auxiliary equations: final output and demand final sector
Y   = MOM.Y; 

pf  = MOM.FG^(-1/eppse)*pg-tauf; % optimality energy
pe  = ((pf+tauf).^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); % definition prices
E   = MOM.EpeY*Y/pe; 
N   = (pe/pn)^eppsy*(1-deltay)/deltay*E; % optimality final good
F   = E* (1+(1/MOM.FG)^((eppse-1)/eppse))^(-eppse/(eppse-1)); % uses optimality energy producers
G   = F/MOM.FG; 
Yout = (deltay^(1/eppsy)*E^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N^((eppsy-1)/eppsy))^((eppsy)/(eppsy-1)); 

% 3) equations
f(1) = Yout-Y;                  % determines share of energy in final production
f(2) = (pf+tauf)*F+pn*N+pg*G-Y; % market clearing demand
f(3)=  1-(deltay*pe^(1-eppsy)+(1-deltay)*pn^(1-eppsy))^(1/(1-eppsy)); % final good as numeraire, optimality final good producers => price of final good

end