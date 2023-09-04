function f=laissez_faire(x, params, list, pol, laggs, indic, Emlim,t)
% Model
% equilibrium for one period!
% takes policy as given

% starts to solve the model in 2020-2024;
% i.e. Aj_laggs  refer to 2015-2019 period
% A_lag is except for the initial condition is 
%- read in policy and parameters
read_in_params;
taul    = pol(:,list.pol=='taul');
%taus    = pol(:,list.pol=='taus'); % earmarking: determined by carbon tax revenues
tauresg = pol(:,list.pol=='tauresg');
taurese = pol(:,list.pol=='taurese');

if indic.limit_LF==0 % otherwise, tauf is choice variable
    tauf    = pol(:,list.pol=='tauf');
else
    tauf=x(list.choice=='tauf');
end

%- initial condition
Af_lag=laggs(list.init=='Af0');
An_lag=laggs(list.init=='An0');
Ag_lag=laggs(list.init=='Ag0');

% choice variables
%- transform variables directly instead of in code

H      = upbarH/(1+exp(x(list.choice=='H')));
w      = exp(x(list.choice=='w'));
Ln     = exp(x(list.choice=='Ln'));
Lg     = exp(x(list.choice=='Lg'));
Lf     = exp(x(list.choice=='Lf'));
gammal = (x(list.choice=='gammal')).^2;

C      = exp(x(list.choice=='C'));
F      = exp(x(list.choice=='F'));
G      = exp(x(list.choice=='G'));

Af     = exp(x(list.choice=='Af'));
Ag     = exp(x(list.choice=='Ag'));
An     = exp(x(list.choice=='An'));
sff    = exp(x(list.choice=='sff'));
sn     = exp(x(list.choice=='sn'));
sg     = exp(x(list.choice=='sg'));
ws     = exp(x(list.choice=='ws'));

pg     = exp(x(list.choice=='pg'));
pn     = exp(x(list.choice=='pn'));
pee    = exp(x(list.choice=='pee'));
pf     = exp(x(list.choice=='pf'));
Trans      = exp(x(list.choice=='Trans')); 

%% - read in auxiliary equations
A_lag   = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);

muu   = C.^(-thetaa); % same equation in case thetaa == 1
   
E     = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));

N     = (1-deltay)/deltay.*(pee./pn)^(eppsy).*E; % demand N final good producers 
Y     = (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy).*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % production function Y 

if indic.notaul==2 % then earmarking on machines of green sector
    xg    = pg.*alphag.*G+tauf.*F; % from firm demand for green machines and gov budget
    taus  = tauf.*F./xg; % from government budget to use carbon tax revenues as green technology adoption
else
    taus =zeros(size(F));
end
%% model equations
q=0;

%1- household optimality (muu auxiliary variable determined above)
%-- budget
q=q+1;
f(q)= -C+(1-taul).*w.*H+Trans;
%-- FOC labor (muu determined from FOC consumption) 
q=q+1;
f(q)= chii*H.^sigmaa-(muu.*(1-taul).*w-gammal);

%- Production
%- output fossil
q=q+1;
f(q) = (alphaf.*pf).^(alphaf./(1-alphaf)).*Af.*Lf -F; 
%5- output neutral
q=q+1;
f(q) = N-An.*Ln.*(pn.*alphan).^(alphan./(1-alphan)); 
%6- output green
q=q+1;
f(q)=  G-Ag.*Lg.*(pg.*alphag./(1-taus)).^(alphag./(1-alphag));

%- Research
% demand scientists

q=q+1;
f(q)= ws.*(1-taurese) - (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*F.*(1-alphaf).*Af_lag)./(rhof^etaa.*Af); 
q=q+1;
f(q)= ws.*(1-tauresg).*(1-taurese).*(1-taus) - (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G.*(1-alphag).*Ag_lag)./(rhog^etaa.*Ag);
q=q+1;
f(q)= ws - (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N.*(1-alphan).*An_lag)./(rhon^etaa.*An);

%- LOM technology
q=q+1;
f(q) = An-An_lag.*(1+gammaa.*(sn./rhon).^etaa.*(A_lag./An_lag).^phii);
q=q+1;
f(q) = Af-Af_lag.*(1+gammaa.*(sff./rhof).^etaa*(A_lag./Af_lag).^phii);
q=q+1;
f(q) = Ag-Ag_lag.*(1+gammaa.*(sg./rhog).^etaa*(A_lag./Ag_lag).^phii);

%- labor demand
q=q+1; 
f(q) = w  - (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*(pf).^(1/(1-alphaf)).*Af; % labour demand fossil
q=q+1;
f(q) = Ln -pn.*(1-alphan).*N./w; % labour demand neutral 
q=q+1;
f(q) = Lg - pg.*(1-alphag).*G./w; % independent of taus! indirectly respected through G

%- Definition/ Optimal Prices
%- optimality energy producers
q=q+1;
f(q) = (pf+tauf).*F.^(1/eppse)- (G).^(1/eppse).*pg; 
q=q+1;
f(q) = pee - ((pf+tauf).^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); 
q=q+1;
f(q) =  1-(deltay.*pee.^(1-eppsy)+(1-deltay).*pn.^(1-eppsy)).^(1/(1-eppsy));

%- market clearing 
q=q+1;
f(q) = Lf+Lg+Ln-H;
q=q+1;
f(q)= gammal.*(H-upbarH);

q=q+1;
f(q)= sff+sg+sn-S; % determines wage in neutral sector

% balanced budget government
q=q+1;
if indic.notaul~=2 % it is env. tax revenues are not redistributed lump sum
    f(q)= -Trans+taul.*w.*H+tauf.*F;
else
    f(q)= -Trans+taul.*w.*H;
end

if indic.limit_LF==1
    q=q+1;
    if t==1 % base year period, tauf =0
        f(q)= tauf;
    else
        f(q)=omegaa*F-deltaa-Emlim;
    end
end

% fprintf('number equations: %d; number variables %d', q, length(list.choice));
end