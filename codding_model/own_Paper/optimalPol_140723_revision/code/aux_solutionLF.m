function LF_t=aux_solutionLF( SLF,pol, laggs, list, symms, indexx, params, indic, Emlim, t)

% output
% LF_t: column vector of simulated variables in period t

%- params
read_in_params;
taul    = pol(:,list.pol=='taul');
tauresg = pol(:,list.pol=='tauresg');
taurese = pol(:,list.pol=='taurese');

if indic.limit_LF==0 % otherwise, tauf is choice variable
    tauf    = pol(:,list.pol=='tauf');
else 
    tauf=SLF.tauf;
end

% read in vars
gammal=SLF.gammal;

Lg=SLF.Lg;
Ln=SLF.Ln;
Lf=SLF.Lf;
w=SLF.w;
H=SLF.H;

F=SLF.F;
G=SLF.G;
C=SLF.C;
Af=SLF.Af;
Ag=SLF.Ag;
An=SLF.An;

sff=SLF.sff;
sg=SLF.sg;
sn=SLF.sn;
ws=SLF.ws;

pg=SLF.pg;
pn=SLF.pn;
pee=SLF.pee;
pf=SLF.pf;
Trans=SLF.Trans;

%- auxiliary equations
E  = (SLF.F^((eppse-1)/eppse)+SLF.G^((eppse-1)/eppse))^(eppse/(eppse-1)); 
N  =  (1-deltay)/deltay.*(SLF.pee./SLF.pn)^(eppsy).*E; % demand N final good producers 
Y  = (deltay^(1/eppsy)*E^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N^((eppsy-1)/eppsy))^(eppsy/(eppsy-1));


xn = SLF.pn*alphan*N;
xf = SLF.pf*alphaf*SLF.F;

if indic.notaul==2 % then earmarking on machines of green sector
    xg    = SLF.pg.*alphag.*SLF.G+tauf.*SLF.F; % from firm demand for green machines and gov budget
    taus  = tauf.*F./xg; % from government budget to use carbon tax revenues as green technology adoption
else
    xg    = SLF.pg.*alphag.*SLF.G; 
    taus  = zeros(size(F));
end

A   = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);

muu      = C.^(-thetaa); % same equation in case thetaa == 1
w     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*(pf).^(1/(1-alphaf)).*Af; 
Emnet     = omegaa*F-deltaa; % net emissions

% utility
if thetaa~=1
    Utilcon = (C.^(1-thetaa))./(1-thetaa);
elseif thetaa==1
    Utilcon = log(C);
end
Utillab = chii.*(H.^(1+sigmaa))./(1+sigmaa);
SWF = Utilcon-Utillab;

% test market clearing

Cincome=Y-xn-xf-xg;

diff=C-Cincome;

if max(abs(diff))>1e-6
    error('market clearing does not hold')
else
    fprintf('goods market cleared!')
end

% test variables read in properly
xx=eval(symms.choice);
guess_trans=trans_guess(indexx('LF'), xx, params, list.params);       
f=laissez_faire(guess_trans, params, list, pol, laggs, indic, Emlim,t);



if (max(abs(f)))>1e-7
    fprintf('f only solved at less than 1e-7, max abs f = %f',max(abs(f)) )
else
    fprintf('saved variables are correct!')
end

% save stuff
LF_t= eval(symms.allvars)';

end
 