function [Sparams, Spol, params, pol,...
    x0LF, SL, SP, SR, Sall, ...
    Sinit201014, init201014 , Sinit201519, init201519]=calibration_results(symms,trProd, trLab, trRes, paramss, list, poll, MOM, indic) 

% parameters
read_in_pars_calib;

%soltions
cell_par=arrayfun(@char, symms.calib, 'uniform', 0);
SL=cell2struct(num2cell(trLab), cell_par, 2);
cell_par=arrayfun(@char, symms.prod, 'uniform', 0);
SP=cell2struct(num2cell(trProd), cell_par, 2);
cell_par=arrayfun(@char, symms.calibRes, 'uniform', 0);
SR=cell2struct(num2cell(trRes), cell_par, 2);

%- parameters
deltay=SP.deltay;
omegaa=SP.omegaa;
chii = SL.chii;
gammaa=SR.gammaa;

%-- save results
params = eval(subs(symms.params));
taurese = 0;
tauresg = 0;
taus = 0;
pol    = eval(symms.pol);

cell_par=arrayfun(@char, symms.params, 'uniform', 0);
Sparams=cell2struct(num2cell(params), cell_par, 2);
cell_par=arrayfun(@char, symms.pol, 'uniform', 0);
Spol=cell2struct(num2cell(pol), cell_par, 2);

%-Allocation
%-- labor
Trans = SL.Trans;
w = SL.w;
H = MOM.targethour;
chii = SL.chii;
gammal = SL.gammal;

%-- production and consumption
[C, Lnw, Lgw, Lfw, pf, F, pn, pg, pee, E, Y, N, G, xn, xg, xf, ...
            AfLf, AgLg, AnLn, omegaa, deltay]=results_production(list, trProd, MOM, ...
            paramss, poll, 'calib');
Ln=Lnw./w;
Lg=Lgw./w;
Lf=Lfw./w;
muu = C^(-thetaa);

%-- research
Af   = AfLf/Lf; 
Ag   = AgLg/Lg; 
An   = AnLn/Ln; 

%--- to save as initial values
        Af0=Af; An0=An; Ag0=Ag;
        init201519 = eval(symms.init);
        clearvars Af0 Ag0 An0

sff = SR.sff;
sg = SR.sg;
sn = SR.sn;
ws= SR.ws;

Af0 = SR.Af_lag; % 2010-2014
Ag0 = SR.Ag_lag;
An0 = SR.An_lag;

A  = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);
A0  = (rhof*Af0+rhon*An0+rhog*Ag0)/(rhof+rhon+rhog);
Agtest= Ag0*(1+gammaa*(sg/rhog)^etaa*(A0/Ag0)^phii); 
 
 if abs(Ag-Agtest)>1e-10
     error('growth rate off')
 else
     fprintf('growth rate works!!')
 end

 %-- remaining variables
GovRev  = MOM.Debt;
Emnet   = omegaa*F-deltaa; % net emissions

if thetaa~=1
 Utilcon = (C.^(1-thetaa))./(1-thetaa);
elseif thetaa==1
 Utilcon = log(C);
end

Utillab = chii.*H.^(1+sigmaa)./(1+sigmaa);
SWF = Utilcon-Utillab;

%-- save variables and initial values
cell_par=arrayfun(@char, symms.allvars, 'uniform', 0);
Sall=cell2struct(num2cell(eval(symms.allvars)), cell_par, 2);

init201014 = eval(symms.init);
cell_par=arrayfun(@char, symms.init, 'uniform', 0);
Sinit201014=cell2struct(num2cell(init201014), cell_par, 2);
Sinit201519=cell2struct(num2cell(init201519), cell_par, 2);

x0LF= eval(symms.choice);

% tests
Cincome= (1-taul)*H*w+Trans;

if abs(Cincome-C)>1e-10
    error('goods market does not clear!')
else
    fprintf('goods market clears!')
end

if abs(H-Ln-Lf-Lg)>1e-10
    error('labor market does not clear!')
else
    fprintf('Labor market clears!')
end
end