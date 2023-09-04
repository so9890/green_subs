function f=laissez_faire_VECT(x, params, list, varrs, laggs,T, indic, Ems)

% Policy version with GOV=tauf pf F
% can handle version with and without taul

% called by script 'test_results.m'
% takes optimal policy results as input

% starts to solve the model in 2020-2024;
% i.e. Aj_laggs  refer to 2015-2019 period
% A_lag is except for the initial condition is 
%- read in policy and parameters
read_in_params;

if indic.limit_LF==0
    tauf=varrs(list.allvars=='tauf', :)'; 
else
    tauf=(x((find(list.choice=='tauf')-1)*T+1:(find(list.choice=='tauf'))*T));
end

taul    = varrs(list.allvars=='taul', :)'; 
taus    = varrs(list.allvars=='taus', :)'; 
tauresg = varrs(list.allvars=='tauresg', :)'; 
taurese = varrs(list.allvars=='taurese', :)'; 


% choice variables
 w     = exp(x((find(list.choice=='w')-1)*T+1:(find(list.choice=='w'))*T));
 H     = upbarH./(1+exp(x((find(list.choice=='H')-1)*T+1:find(list.choice=='H')*T)));
 Lf    = exp(x((find(list.choice=='Lf')-1)*T+1:(find(list.choice=='Lf'))*T));
 Lg    = exp(x((find(list.choice=='Lg')-1)*T+1:(find(list.choice=='Lg'))*T));
 Ln    = exp(x((find(list.choice=='Ln')-1)*T+1:(find(list.choice=='Ln'))*T));

 C      = exp(x((find(list.choice=='C')-1)*T+1:(find(list.choice=='C'))*T));
 F      = exp(x((find(list.choice=='F')-1)*T+1:(find(list.choice=='F'))*T));
 G      = exp(x((find(list.choice=='G')-1)*T+1:(find(list.choice=='G'))*T));
 Af     = exp(x((find(list.choice=='Af')-1)*T+1:(find(list.choice=='Af'))*T));
 Ag     = exp(x((find(list.choice=='Ag')-1)*T+1:(find(list.choice=='Ag'))*T));
 An     = exp(x((find(list.choice=='An')-1)*T+1:(find(list.choice=='An'))*T));
 
        sff    = (x((find(list.choice=='sff')-1)*T+1:(find(list.choice=='sff'))*T)).^2;
        sg     = (x((find(list.choice=='sg')-1)*T+1:(find(list.choice=='sg'))*T)).^2;
        sn     = (x((find(list.choice=='sn')-1)*T+1:(find(list.choice=='sn'))*T)).^2;
        ws     = (x((find(list.choice=='ws')-1)*T+1:(find(list.choice=='ws'))*T)).^2;
   
 gammal = x((find(list.choice=='gammal')-1)*T+1:(find(list.choice=='gammal'))*T).^2;
 pg     = exp(x((find(list.choice=='pg')-1)*T+1:(find(list.choice=='pg'))*T));
 pn     = exp(x((find(list.choice=='pn')-1)*T+1:(find(list.choice=='pn'))*T));
 pee     = exp(x((find(list.choice=='pee')-1)*T+1:(find(list.choice=='pee'))*T));
 pf     = exp(x((find(list.choice=='pf')-1)*T+1:(find(list.choice=='pf'))*T));
 Trans     = exp(x((find(list.choice=='Trans')-1)*T+1:(find(list.choice=='Trans'))*T));

 %% - read in auxiliary equations
%- initial condition
Af_lag=[laggs(list.init=='Af0'); Af(1:T-1)];
An_lag=[laggs(list.init=='An0'); An(1:T-1)];
Ag_lag=[laggs(list.init=='Ag0'); Ag(1:T-1)];
A_lag   = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);


muu      = C.^(-thetaa); % same equation in case thetaa == 1
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));
N       =  (1-deltay)/deltay.*(pee./pn).^(eppsy).*E; % demand N final good producers 
Y       =  (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy).*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % production function Y 

% xn      = (alphan*pn).^(1/(1-alphan)).*Ln*An;
% xf      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*Lf*Af;
% xg      = (alphag*pg).^(1/(1-alphag)).*Lg*Ag;

%% model equations
q=0;

%1- household optimality (muu auxiliary variable determined above)

q=q+1;
f(q:T)= chii*H.^sigmaa-(muu.*(1-taul).*w-gammal);
%3- budget
q=q+1;
f((q-1)*T+1:T*q) = -C+(1-taul).*w.*H+Trans;

%4- output fossil
q=q+1;
f((q-1)*T+1:T*q) = (alphaf.*pf).^(alphaf./(1-alphaf)).*Af.*Lf -F; 

%5- output neutral
q=q+1;
f((q-1)*T+1:T*q) = N-An.*Ln.*(pn.*alphan).^(alphan./(1-alphan)); 

%6- output green
q=q+1;
f((q-1)*T+1:T*q)=  G-Ag.*Lg.*(pg.*alphag./(1-taus)).^(alphag./(1-alphag));

%7- demand green scientists

q=q+1;
f((q-1)*T+1:T*q)= ws.*(1-taurese) - (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*F.*(1-alphaf).*Af_lag)./(rhof^etaa.*Af); 
%8
q=q+1;
f((q-1)*T+1:T*q)= ws.*(1-tauresg).*(1-taurese).*(1-taus) - (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G.*(1-alphag).*Ag_lag)./(rhog^etaa.*Ag);
%9
q=q+1;
f((q-1)*T+1:T*q)= ws - (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N.*(1-alphan).*An_lag)./(rhon^etaa.*An);
 
%10- LOM technology
q=q+1;
f((q-1)*T+1:T*q) = An-An_lag.*(1+gammaa.*(sn./rhon).^etaa.*(A_lag./An_lag).^phii);
%11
q=q+1;
f((q-1)*T+1:T*q) = Af-Af_lag.*(1+gammaa.*(sff./rhof).^etaa.*(A_lag./Af_lag).^phii);
%12
q=q+1;
f((q-1)*T+1:T*q) = Ag-Ag_lag.*(1+gammaa.*(sg./rhog).^etaa.*(A_lag./Ag_lag).^phii);
 
q=q+1; % labour demand fossil
f((q-1)*T+1:T*q) =  w  - (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*(pf).^(1/(1-alphaf)).*Af; % labour demand fossil
q=q+1;
f((q-1)*T+1:T*q) = Ln -pn.*(1-alphan).*N./w; % labour demand neutral 
q=q+1;
f((q-1)*T+1:T*q) = Lg - pg.*(1-alphag).*G./w;
% prices and wages
%19- optimality energy producers
q=q+1;
f((q-1)*T+1:T*q) = (pf+tauf).*F.^(1/eppse)- (G).^(1/eppse).*pg; 

%- definitions prices
%20
q=q+1;
f((q-1)*T+1:T*q) = pee - ((pf+tauf).^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); %definition
%21
q=q+1;
f((q-1)*T+1:T*q) =  1-(deltay.*pee.^(1-eppsy)+(1-deltay).*pn.^(1-eppsy)).^(1/(1-eppsy));

%22- market clearing (consumption good=> numeraire)
q=q+1;
f((q-1)*T+1:T*q) = Lf+Lg+Ln-H;  
q=q+1;
f((q-1)*T+1:T*q)= sff+sg+sn-S; % determines wage in neutral sector
q=q+1;
f((q-1)*T+1:T*q)= gammal.*(H-upbarH);


% income schedule budget clearing
q=q+1;
if indic.notaul<2
    f((q-1)*T+1:T*q)= -Trans+taul.*w.*H+tauf.*F;
else
    f((q-1)*T+1:T*q)= -Trans+taul.*w.*H;    
end
if indic.limit_LF==1
    q=q+1;
    f((q-1)*T+1:T*q)=omegaa*F-deltaa-Ems';
end

%fprintf('number equations: %d; number variables %d', q, length(list.choice));
end