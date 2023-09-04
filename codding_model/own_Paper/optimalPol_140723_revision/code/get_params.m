function [params, Sparams,  pol, init201014, init201519, list, symms, Ems,  Sall, x0LF, MOM , indexx , StatsEms]...
            =get_params( T, indic, lengthh)

% calibration without emission limit
 indic.limit_LF=0; 
   
% function to read in parameter values and to calibrate direct parameters
% calls on calibration_matching and calibration_emissions 

% input
% T       : length of optimization periods
% indic   : indicators to choose model version
% lengthh : length of model period

% output
% params:       numeric vector of calibrated parameters and initial conditions
% pol:          numeric vector of policy
% targets:      numeric vector of emission targets
% init:         numeric vector of initial conditions
% Sall:         structure of all model variables in baseyear
% Ems:          numeric vector of emission targets
% x0LF:         numeric vector of LF solution in baseyear ordered as in list.choice (used for LF solution) 
% 
%% symbolic vector and list
syms sigmaa...      % 1/sigmaa = Frisch elasticity of labour
     thetaa...      % courvature consumption utility
     betaa...       % discount factor 5 years
     chii ...       % disutility labour
     upbarH...      % time endowment 
     alphaf ...     % machine share fossil
     alphan ...     % machine share neutral
     alphag ...     % machine share green
     eppsy ...      % elasticity final good production
     eppse ...      % elasticity energy production
     deltay ...     % weight on energy in final production
     gammaa ...     % quality research
     etaa ...       % returns to scale scientists
     rhof ...       % number tasks scientists fossil
     rhon ...       % number tasks scientists neutral
     rhog ...       % number tasks scientists green
     phii ...       % spill over curvature
     Ems ...        % vector of net emission targets
     deltaa ...     % regeneration rate nature
     omegaa ...     % emission share of dirty output
     Af0 ...        % initial technology level fossil
     Ag0 ...        % initial technology level green
     An0 ...        % initial technology level neutral
     S ...          % number of scientists
     real 
 
syms taul ...       % income tax progressivity
     tauf ...       % sales tax fossil
     taus ...       % green sector subsidies (earmarking)
     taurese ...    % research subsidies energy
     tauresg ...    % additional research subsidies green
     real
 
symms.params = [sigmaa, thetaa, betaa, chii, upbarH, alphaf, alphan, alphag, S,...
                eppsy, eppse, deltay,gammaa, etaa, rhof, rhon, rhog, phii, deltaa, omegaa];   
list.params  = string(symms.params);

symms.init   = [Af0, An0, Ag0];
list.init    = string(symms.init);

symms.pol     = [taul, tauf, taus, taurese, tauresg];
list.pol      = string(symms.pol);

% parameters directly calibrated
symms.paramsdir = [sigmaa, thetaa, betaa, upbarH, alphaf, alphan, alphag,S,...
                   eppsy, eppse, etaa, phii,  rhof, rhon, rhog, deltaa];   
list.paramsdir  = string(symms.paramsdir);

symms.poldir     = [taul, tauf];
list.poldir      = string(symms.poldir);
%% Calibration 
sigmaa   = 1/0.75;      % from Chetty et al 

%if indic.util== 0
thetaa   = 1; % log utility
% else
%     if indic.Bop==1 % income effect dominates! => labor supply less responsive
%         thetaa=(0.2+1/sigmaa)./((1-0.2)/sigmaa);
%     else
%         thetaa=0.4;
%     end
% end

betaa    = (.985)^5;  % Barrage, but here for 5 years
upbarH   = 1;
S        = 0.01;
eppse    = 1.5;            % Fried
eppsy    = 0.05;           % Fried
alphaf   = 1-0.28;         % Fried: fossil has a higher labour share!
alphag   = 1-0.09;         % Fried
alphan   = 1-0.64;         % Fried
 
gammaa   = 3.96;           % Fried
%if indic.spillovers==0
etaa     = 0.79;% Fried returns to research
%else
%    etaa     = 1.2; % positive spillovers=> to accomodate zero scientists in competitive eqbm 
%end
rhof     = 0.01; % Fried
rhon     = 1;    % Fried
rhog     = 0.01; % Fried
phii     = 0.5;  % Fried: 0.5; knowledge spillovers! 

%- policies
taul    = 0.24; % linear tax! As in Barrage 2020
tauf    = 0; 

%% - direct calibration 
%-- get moments
MOM = calibration_moments();
% MOM.S = 0.01; % from fried: Supply scientists in base year
MOM.growth = (1.017795)^5 -1; %5 year grwoth rate from OECD over initial period
MOM.Debt=0; %GovRev;  % as share of GDP
MOM.GDP1519MILLION=101950887.298532; %sum GDP over 2015-2019 expressedn in 2019 preisen
%102140824.682949; %98280493; 

%% - emissions
[deltaa, Ems, MOM, StatsEms]= calibration_emissions(T, lengthh, MOM); 
Ems  = StatsEms.Emslimit_constantEmsRat_Budget; % use equal allocation of remaining carbon budget

% -omegaa follows in main calibration
%% save directly calibrated variables
parsHelp = eval(symms.paramsdir);
polhelp= eval(symms.poldir);
%%
[x0LF, SL, SP, SR, Sall, Sinit201014, init201014 , Sinit201519, init201519, Sparams, Spol, params, pol, symms, MOM, indexx, list]...
 = calibration_matching(MOM, symms, list, parsHelp, polhelp, indic);

%- TFP

Sparams.TFPF0= init201014(list.init=='Af0')^(1-Sparams.alphaf);
Sparams.TFPG0= init201014(list.init=='Ag0')^(1-Sparams.alphag);
Sparams.TFPN0= init201014(list.init=='An0')^(1-Sparams.alphan);


end