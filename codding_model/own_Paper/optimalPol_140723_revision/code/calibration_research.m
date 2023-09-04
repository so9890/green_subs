function  [ceq]= calibration_research(x, MOM, list, trProd, paramss, poll, Af, An, Ag)

% this function calibrates the research side of the model
% Af_lag, An_lag, Ag_lag
% gammaa
%... and matches the allocation of scientists and scientists wage rate in
% eqbm:
% sff, sn, sg
% ws

read_in_pars_calib;

% calibration to 2019 (lag = 2010-2014)
Af_lag  = exp(x(list.calibRes=='Af_lag'));
Ag_lag  = exp(x(list.calibRes=='Ag_lag'));
An_lag  = exp(x(list.calibRes=='An_lag'));
sff     = exp(x(list.calibRes=='sff'));
sg      = exp(x(list.calibRes=='sg'));
sn      = exp(x(list.calibRes=='sn'));
ws      = exp(x(list.calibRes=='ws'));
gammaa  = exp(x(list.calibRes=='gammaa'));



% auxiliary
A =  (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);
A_lag  = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);

% required vars from previous calibration steps
[ C, ~, ~, ~, pf, FF, pn, pg, ~, ~, ~, N, G,...
    ~, ~, ~, ~, ~, ~, ~, ~]=results_production(list, trProd, MOM , paramss, poll, 'calib'); 
muu=C^(-thetaa);

%% target equations

q=0;
% % target gammaa
   q=q+1;
   ceq(q)= gammaa -3.96;%(A/A_lag-1)-MOM.growth; % targeting 5 year growth rate
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
ceq(q)= ws - (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*FF.*(1-alphaf).*Af_lag)./(rhof^etaa.*Af); 
%8
q=q+1;
ceq(q)= ws - (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G.*(1-alphag).*Ag_lag)./(rhog^etaa.*Ag);
%9
q=q+1;
ceq(q)= ws - (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N.*(1-alphan).*An_lag)./(rhon^etaa.*An);

q=q+1;
ceq(q)= sff+sg+sn-S;
% q=q+1;   
% ceq(q)= chiis.*MOM.targethour.^(sigmaas+taul)-muu.*lambdaa.*(1-taul).*ws.^(1-taul);
            
%ceq(q)= MOM.targethour-(muu*ws/(chiis)).^(1/sigmaas);  % equal disutility as for other labour => pins down ws
%  q=q+1;
%  ceq(q)= MOM.rhon-rhon;  % equal disutility as for other labour => pins down ws

% fprintf('number equations %d, number variables %d', q, length(x));
end