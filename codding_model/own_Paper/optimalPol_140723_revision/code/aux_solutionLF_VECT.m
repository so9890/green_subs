function LF_SIM=aux_solutionLF_VECT(x, list, symms,varrs, params, T , indic)

read_in_params;

if indic.limit_LF==0
    tauf=varrs(list.allvars=='tauf', :)'; 
else
    tauf=(x((find(list.choice=='tauf')-1)*T+1:(find(list.choice=='tauf'))*T));
end

gammal = x((find(list.choice=='gammal')-1)*T+1:(find(list.choice=='gammal'))*T).^2;

 Ln    = exp(x((find(list.choice=='Ln')-1)*T+1:(find(list.choice=='Ln'))*T));
 Lg    = exp(x((find(list.choice=='Lg')-1)*T+1:(find(list.choice=='Lg'))*T));
 Lf    = exp(x((find(list.choice=='Lf')-1)*T+1:(find(list.choice=='Lf'))*T));
 w     = exp(x((find(list.choice=='w')-1)*T+1:(find(list.choice=='w'))*T));
 H     = upbarH./(1+exp(x((find(list.choice=='H')-1)*T+1:find(list.choice=='H')*T)));

C      = exp(x((find(list.choice=='C')-1)*T+1:(find(list.choice=='C'))*T));

 F      = exp(x((find(list.choice=='F')-1)*T+1:(find(list.choice=='F'))*T));
 G      = exp(x((find(list.choice=='G')-1)*T+1:(find(list.choice=='G'))*T));
 Af     = exp(x((find(list.choice=='Af')-1)*T+1:(find(list.choice=='Af'))*T));
 Ag     = exp(x((find(list.choice=='Ag')-1)*T+1:(find(list.choice=='Ag'))*T));
 An     = exp(x((find(list.choice=='An')-1)*T+1:(find(list.choice=='An'))*T));
sff     = (x((find(list.choice=='sff')-1)*T+1:(find(list.choice=='sff'))*T)).^2;
sg     = (x((find(list.choice=='sg')-1)*T+1:(find(list.choice=='sg'))*T)).^2;
sn     = (x((find(list.choice=='sn')-1)*T+1:(find(list.choice=='sn'))*T)).^2;
ws     = (x((find(list.choice=='ws')-1)*T+1:(find(list.choice=='ws'))*T)).^2;


 pg     = exp(x((find(list.choice=='pg')-1)*T+1:(find(list.choice=='pg'))*T));
 pn     = exp(x((find(list.choice=='pn')-1)*T+1:(find(list.choice=='pn'))*T));
 pee     = exp(x((find(list.choice=='pee')-1)*T+1:(find(list.choice=='pee'))*T));
 pf     = exp(x((find(list.choice=='pf')-1)*T+1:(find(list.choice=='pf'))*T));
Trans = (x((find(list.choice=='Trans')-1)*T+1:(find(list.choice=='Trans'))*T));
taul=varrs(list.allvars=='taul', :)';
taus=varrs(list.allvars=='taus', :)';
tauresg=varrs(list.allvars=='tauresg', :)';
taurese=varrs(list.allvars=='taurese', :)';
% auxiliary variables 
muu      = C.^(-thetaa); % same equation in case thetaa == 1
S=(sg+sff+sn);
  

E  = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1)); 
N  =  (1-deltay)/deltay.*(pee./pn).^(eppsy).*E; % demand N final good producers 
Y  = (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy).*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1));
xn = pn.*alphan.*N;
xg = pg.*alphag.*G;
xf =  pf.*alphaf.*F;

A  = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);

Emnet  = omegaa.*F-deltaa; % net emissions

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

if max(abs(diff))>1e-7
    error('market clearing does not hold')
else
    fprintf('goods market cleared!')
end

% % test variables read in properly
%  xx=eval(symms.choice);
%  xx=xx(:);
%  guess_trans=trans_guess(indexx('LF'), xx, params, list.params);
%  f=laissez_faire_VECT(guess_trans, params, list, varrs, laggs,T, indic, Ems);
% 
%  if (max(abs(f)))>1e-10
%      fprintf('f only solved at less than 1e-10')
%  else
%      fprintf('saved variables are correct!')
%  end

% save stuff
LF_SIM= eval(symms.allvars);
end
 