indic.util =0; % ==0 log utility, otherwise as in Boppart
indic.Bop=0; % indicator ==1 then uses version as discussed in Boppart: 
                 % income effect stronger than substitution effect and
                 % thetaa > 1
indic.elasE = 0; % ==0 then standard as in Fried, % == 1 then eppsee=10
indic.sigmaWorker =0; % ==0 then standard calibration with chetty sigmaa=1/0.75; ==1 then higher frish elasticity: sigma smaller! sigmaa= 1/1.5;
                      % == 2 then smaller frisch elasticity sigmaa= 1/0.5; 
indic.targetWhat =0; % uses baseline emission limit (equal allocation of remaining carbon budget)
                      
%indic.sep =0; % ==0 one joint market (in calibration very low fossil and green scientists to satisfy wage clearing 
              % ==1 3 separate markets 
              % ==2 if partial equbm; relative to joint market
              % ==3 energy market joint and non-energy market separate
%indic.Sun = 2; %==2 then scientists are taxed too!  ; ==0 then scientsist form part of household; ==1 then scientists are separate households             
indic.target =0; % ==1 if uses emission target
indic.know_spill =0; % ==0 knowledge spillovers as in Fried =0.5; ==1 then without;
                        % ==2 then smaller: phii=0.25; ==3 then bigger: phii=0.75
%indic.sizeequ=0; %==1 then research sectors have same size => is there still a higher progressive tax when there are spillovers?
%indic.spillovers =0; % ==1 then there are positive spillover effects of scientists within sectors! 
%indic.noskill = 0; % == 1 if no skill calibration of model
indic.notaul=0; % Indicator of policy
                % ==0 baseline: lump sum, taul is an option 
                % ==1 lump sum trans, no taul
                % ==2 use env tax revenues as research subsidies (earmarking)
indic.limit_LF =0; % ==1 then tauf is determined by meeting limit in each period
                   %  set by a planner who knows how economy works but each
                   %  period; not dynamic! (in optimal policy taking dynamics into account)
%indic.taus =0; %==0 then no subsedy on green 
indic.subsres = 0; % == 0 no lump-sum financed research subsidies,
                   % == 1 with lump-sum financed research subsidies
indic.xgrowth=0;
%indic.extern=0; % extern==0 when uses no externality in utility
% but ensure no externality when target is used 
%if indic.target==1
%    indic.extern=0;
%end
indic.count_techgap=0; % if ==1 then uses technology gap as in Fried
%indic.subs = 0; %==1 eppsy>1 (energy and neutral good are substitutes)
indic.PV = 1; % ==1 if continuation value is added to planners problem
%indic.PVwork =0; %==0 then disutility of work is not in 
%indic.emsbase=0; % ==0 then uses emission limits as calculated
%indic.zero = 0;  % ==1 then version without growth
%indic.GOV=0; % ==0 then no gov revenues
indic.taul0=0; %==0 then calibrated value for taul; ==1 then 0
%indic.labshareequ =0; %==1 then equal capital shares in fossil green and non-energy sectors
indic

percon = 0;  % periods nonconstrained before 50\% constrained
