function  [ceq]= calibRem_nows_fsolve_GOOD_nonsep(x, MOM, list, trProd, paramss, poll, Af, An, Ag, lambdaa)
% this function backs out missing functions:
% 

read_in_pars_calib

% calibration to 2019 (lag = 2010-2014)
Af_lag  = exp(x(list.calib3=='Af_lag'));
Ag_lag  = exp(x(list.calib3=='Ag_lag'));
An_lag  = exp(x(list.calib3=='An_lag'));
sff     = exp(x(list.calib3=='sff'));
sg      = exp(x(list.calib3=='sg'));
sn      = exp(x(list.calib3=='sn'));
ws      = exp(x(list.calib3=='ws'));
% sigmaas = exp(x(list.calib3=='sigmaas'));
chiis   = exp(x(list.calib3=='chiis'));
gammaa  = exp(x(list.calib3=='gammaa'));
% rhon  = exp(x(list.calib3=='rhon'));
% rhog  = exp(x(list.calib3=='rhog'));
% rhof  = exp(x(list.calib3=='rhof'));


% auxiliary
% rhog=rhon/MOM.rhong;
% rhof=rhon/MOM.rhonf;
A =  (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);
A_lag  = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);

% required vars from previous calibration steps
[ C, ~, ~, ~, pf, FF, pn, pg, ~, ~, ~, N, G,...
    ~, ~, ~, ~, ~, ~, ~, ~]=resProd(list, trProd, MOM , paramss, poll, 'calib'); 
muu=C^(-thetaa);
% gammaa = MOM.growth/((sn/rhon)^etaa*(A_lag/An_lag)^phii); % n is growth rate in non-energy technology: n0=An'/An-1
%  q=q+1;
% ceq(q)= MOM.S-(ws/(chiis.*C.^thetaa)).^(1/sigmaas); % hours supply by scientists to find sigmaas
% chiis= ws/(C^thetaa*MOM.S^sigmaas);

%% target equations

q=0;
% % target gammaa
   q=q+1;
   ceq(q)= (A/A_lag-1)-MOM.growth; % targeting 5 year growth rate
%  q=q+1;
%  ceq(q)= gammaa - MOM.growth/((sn/rhon)^etaa*(A_lag/An_lag)^phii);

% ceq(q)=
% q=q+1;
% ceq(q)= gammaa - MOM.growth/((sn/rhon)^etaa*(A_lag/An_lag)^phii);

% LOM => lagged technology
q=q+1;
ceq(q) = Af- Af_lag*(1+gammaa*(sff/rhof)^etaa*(A_lag/Af_lag)^phii);
q=q+1;
ceq(q) = Ag- Ag_lag*(1+gammaa*(sg/rhog)^etaa*(A_lag/Ag_lag)^phii);
 q=q+1;
 ceq(q) = An- An_lag*(1+gammaa*(sn/rhon)^etaa*(A_lag/An_lag)^phii);
% scientist demand
q=q+1;
ceq(q)= ws - (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*(1-tauf).*FF.*(1-alphaf).*Af_lag)./(rhof^etaa.*Af); 
%8
q=q+1;
ceq(q)= ws - (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G.*(1-alphag).*Ag_lag)./(rhog^etaa.*(1-taus)*Ag);
%9
q=q+1;
ceq(q)= ws - (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N.*(1-alphan).*An_lag)./(rhon^etaa.*An);

q=q+1;
ceq(q)= sff+sg+sn-MOM.targethour;
q=q+1;   
ceq(q)= chiis.*MOM.targethour.^(sigmaas+taul)-muu.*lambdaa.*(1-taul).*ws.^(1-taul);
            
%ceq(q)= MOM.targethour-(muu*ws/(chiis)).^(1/sigmaas);  % equal disutility as for other labour => pins down ws
%  q=q+1;
%  ceq(q)= MOM.rhon-rhon;  % equal disutility as for other labour => pins down ws

% fprintf('number equations %d, number variables %d', q, length(x));
end