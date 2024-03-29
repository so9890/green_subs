%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The role of fiscal %%%
%%%%% policies in the env %%
%%%%% policy %%%%%%%%%%%%%%%LF_SIM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Author: Sonja Dobkowitz
% building on Lint Barrage's code ReStud2019
% Version: July 2023
clear
cd '/home/sonja/Documents/DocumentsSonja/projects/Feasible_OptimalClimatePol/codding_model/own_Paper/optimalPol_140723_revision/code'
%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 1: Select Scenario        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = 12;  % Direct optimization period time horizon: 2020-2080
         % one period = 5 years

lengthh = 5; % number of years per period         
read_indicators;% read in indicator values; change version here

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 2: Parameters        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfile(sprintf('../output/params_2407.mat'))
    fprintf('loading parameter values')
    load(sprintf('../output/params_2407'),...
        'params', 'Sparams', 'polCALIB', 'init201014', 'init201519', 'list', 'symms', 'Ems', 'Sall', 'x0LF', 'MOM', 'indexx', 'StatsEms')
else
    fprintf('calibrating model')
    [params, Sparams,  polCALIB,  init201014, init201519, list, symms, Ems,  Sall, x0LF, MOM, indexx, StatsEms]...
        =get_params( T, indic, lengthh);
        % function 1) sets direct parameters, 
        %          2) calibrates model to indirect params.
    save(sprintf('../output/params_2407')) 
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 3: BAU Simulation        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in this section I simulate the economy (BAU) starting from 2015-2019
% order of variables in LF_SIM as in list.allvars
indic.notaul=0; % baseline policy with lump-sum transfers and labor tax
for scen="Limit" %[ "BAU", "LF", "Limit"]
    fprintf("running version %s", scen)
    if scen=="BAU"
        POL=polCALIB;
        indic.limit_LF=0;
    elseif scen == "LF"
        POL=[0,0, 0, 0, 0];
        indic.limit_LF=0;
    elseif scen=="Limit" % with earmarking does not solve
        indic.limit_LF=1;
        POL= polCALIB; % use calibrated labor tax
    end
        % if xgr==0
        %     if indic.count_techgap==0
[LF_SIM, pol, FVAL] = solve_LF(T, list, POL, params, symms, x0LF, init201014, indexx, indic, Ems);
helper.LF_SIM=LF_SIM;
[COMP]=solve_LF_VECT(T, list,  params,symms, init201519, helper, indic, Ems, MOM);
%- save results

save(sprintf('../output/%s_phii%d_sigma%d_Bop%d_util%d.mat',...
            scen, indic.know_spill, indic.sigmaWorker, indic.Bop), 'COMP');
end
% else
                 % indic.limit_LF=1;
                 %% Counterfactual productivity levels
                %- version with counterfactual technology gap
                
                 % An0=init201014(list.init=='An0');
                 % Ag0=0.9*An0;
                 % Af0=Ag0/0.4; 
                 % initcount= eval(symms.init); % vector or counterfactual technology 
                 % iin=load('../output/init_techgap.mat');

                 % iin=load('../output/init_techgap.mat');
                 % [LF_SIM, pol, FVAL] = solve_LF(T, list, POL, params, symms, x0LF, iin.initcount, indexx, indic, Ems);
                 %  helper.LF_SIM=LF_SIM;
                 % [COMP]=solve_LF_VECT(T, list, params,symms, iin.init1519count, helper, indic, Ems, MOM);
            % end
        %- Transforming tauf in per ton of carbon in 2014-19 us dollars
         tauf=COMP(:,list.allvars=='tauf');
         tauf_CO2=tauf./Sparams.omegaa;
         tauf_perton2019 = tauf_CO2*(MOM.GDP1519MILLION*1e6)./(1e9); % denominator to go from gigaton to ton in 2019 prices
         TAUF(:,nnt+1)=tauf_perton2019*1.12; % to have it in 2022 prices

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 4: Sociel Planner allocation                             %%%
% Timing: starting from 2020-2025  as initial period                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indic.know_spill=0;
indic.util=0;
indic.Bop=0;
count=25;
Tinit=T;
indic.sigmaWorker=0;
indic.targetWhat =0; % uses baseline emission limit (equal allocation of remaining carbon budget)

for xgr=0
    indic.xgrowth=xgr;
       for nknk=[0]
            indic.noknow_spill=nknk;
%             if ~isfile(sprintf('SP_target_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, params(list.params=='etaa')))
%                 indic.target=1;
%                 fprintf('solving Social planner solution with target, noskill%d', indic.noskill);
       for tar=[0,1]
            indic.target=tar;
            indic     
            if indic.count_techgap==0
           %SP_solve_extT(list, symms, params, count, init201519, indic, Tinit, Emss, MOM, percon);
                SP_solve(list, symms, params, Sparams, x0LF, init201014, init201519, indexx, indic, T, Ems, MOM, percon);
            else
            %    SP_solve(list, symms, params, Sparams, x0LF, iin.initcount, iin.init1519count, indexx, indic, T, Emss, MOM, percon);

            end
       end
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 5: Solve for Optimal Allocation       %%%
% Timing: starting from 2020-2025 the gov. chooses      %%
% the optimal allocation                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indic.elasE =0;
indic.sigmaWorker=0;
indic.Bop=0;
indic.util=0; % have both to be changed
indic.count_techgap=0;
indic.limit_LF=0; % no need to test this
indic.testT =0; % do not test value of T but only run with T=12
indic.targetWhat = 0; %==0 then baseline, ==1 then equal shares
indic.subsres = 1; %=1 with lump-sum financed research subsidies

for tr =[1]
    indic.target=tr;
for xgr=0
    indic.xgrowth=xgr;
for nknk=[0]
    indic.know_spill=nknk; %==0 then fried,, ==3 then 3/4

 for nnt=[3] % policy regime: 
                % ==0 baseline: lump-sum transfers, taul is an option 
                % ==1 lump-sum trans, no taul
                % ==2 use env tax revenues as research subsidies
                %     (earmarking; feasible policy); 
                %     taul is an option
                % == 3 as 2 but taul is not an option; taul=0
     indic.notaul=nnt;
     indic
 if indic.count_techgap==0
     OPT_solve(list, symms, params, x0LF, init201519, indexx, indic, T, Ems, MOM, percon);
 else
     OPT_solve(list, symms, params, x0LF, iin.init1519count, indexx, indic, T, Ems, MOM, percon);
 end
 end
end
end
end

%%
%-- extend optimality for count

indic.taus  = 0; % with ==0 no taus possible!
indic.sep =0;
indic.extern=0;
indic.GOV=0; % ==0 then no gov revenues
indic.sizeequ=0; 
indic.util=0;
indic.Bop=0;
indic.sigmaWorker=0;
indic.noknow_spill=3;
indic.Sun=2;
indic.targetWhat=1;
indic.limit_LF=0; % no need to test this
indic.testT =0; % do not test value of T
indic.taulFixed=0;
indic
count=30;% addiitonal periods

if indic.targetWhat==0
    Emss=Ems;
else
    Emss=StatsEms.Emslimit_constantEmsRat_Budget;
end

for tr =[1,0]
    indic.target=tr;
for xgr=0
    indic.xgrowth=xgr;
for nsk=[0]
    indic.noskill=nsk;
 for nnt=[4]
     indic.notaul=nnt;
     indic
[symms, list, opt_all]= OPT_solve_sep_ExtT(list, symms, params, x0LF, init201519, indexx, indic, T, Emss, MOM, percon, count);
 end
end
end
end
%%
%-- extend optimality for count
% for 0,3,4 => code with taul generated directly from non-taul alternative
indic.taus  = 0; % with ==0 no taus possible!
indic.sep =0;
indic.extern=0;
indic.util=0;
indic.Bop=00;
indic.GOV=0; % ==0 then no gov revenues
indic.sizeequ=0; 
indic.noknow_spill=3;
indic.limit_LF=0; % no need to test this
indic.testT =0; % do not test value of T
indic.Sun=2;
indic.targetWhat=1;
indic.taulFixed=0;

indic
count=30;% addiitonal periods

if indic.targetWhat==0
    Emss=Ems;
else
    Emss=StatsEms.Emslimit_constantEmsRat_Budget;
end
for tr =[1,0]
    indic.target=tr;
for xgr=0
    indic.xgrowth=xgr;
for nsk=0
    indic.noskill=nsk;
 for nnt=4
     indic.notaul=nnt;
     indic
[symms, list, opt_all]= OPT_solve_sep_ExtT_direct(list, symms, params, x0LF, init201519, indexx, indic, T, Emss, MOM, percon, count);
 end
end
end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 6: Competitive equi 
%%%      counterfactual policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for tf=[1,5] %==4 (set optimal pol without no spil in benchmark model); ==2 then uses policy as in helper.opt_all; only benchmark taul, tauf=0 in other models
         %==5 uses optimal taul from without target as given, derives tauf
         %==6 taul fixed, tauf from joint optimal with target
         %(with limit_LF==1)
indic.tauf=tf; % ==1 uses version with optimal taul=0 but tauf=1; ==0 uses version with tauf=0 optimal but taul =1
indic.PV=1;
indic.notaul=4;
indic.limit_LF=0; % for simulation, not a policy calculation
indic.noknow_spill=3; % counterfactuals so far only without knowledge spillovers
indic.oldCalib=0; % ==0 then uses new calbration
indic.Bop=0;
indic.Sun=2;
indic.targetWhat=0;
indic.target=1;
T=12;
count=30;

if indic.targetWhat==0
    Emss=Ems;
else
    Emss=StatsEms.Emslimit_constantEmsRat_Budget;
end

for xgr=0
    indic.xgrowth=xgr;
    for nsk=0
        indic.noskill=nsk;
% load benchmark policy
if tf>=2 && tf<5 % read in benchmark model wrt skill xgr
    if tf~=4
        if indic.target==1
            helper=load(sprintf('OPT_target_2112_emnet%d_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
            count, indic.targetWhat, indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
  
  %          helper=load(sprintf('OPT_target_plus%d_0501_emnet%d_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
   %         count, indic.targetWhat, indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
  
            %helper=load(sprintf('OPT_target_2112_emnet%d_Sun%d_spillover%d_knspil%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.targetWhat, indic.Sun, indic.spillovers, indic.noknow_spill, indic.noskill, indic.notaul, indic.sep, indic.xgrowth, indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
        elseif indic.target ==0
            helper=load(sprintf('OPT_notarget_plus%d_0501_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
               count, indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul,indic.sep, indic.extern, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa'))); 
        end
    elseif tf==4
         helper=load(sprintf('OPT_target_plus30_0509_spillover%d_knspil1_taus0_noskill0_notaul%d_sep%d_xgrowth0_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.spillovers, indic.notaul, indic.sep, indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
    end
elseif tf<2 % load in same model wsrt skill xgr
    if indic.oldCalib==1
        helper=load(sprintf('OPT_target_plus30_0509_spillover%d_knspil%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.spillovers,indic.noknow_spill,indic.noskill, indic.notaul, indic.sep, indic.xgrowth , indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
    elseif indic.oldCalib==0  
        if indic.target==1
            helper=load(sprintf('OPT_target_2112_emnet%d_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
            indic.targetWhat, indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
%    helper=load(sprintf('OPT_target_plus%d_0501_emnet%d_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
%             count, indic.targetWhat, indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
%   
            %helper=load(sprintf('OPT_target_2112_emnet%d_Sun%d_spillover%d_knspil%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.targetWhat, indic.Sun, indic.spillovers, indic.noknow_spill, indic.noskill, indic.notaul, indic.sep, indic.xgrowth, indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
        elseif indic.target ==0
            helper=load(sprintf('OPT_notarget_plus%d_0501_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
               count, indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul,indic.sep, indic.extern, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa'))); 
        end
    end
elseif tf>=5 
    % use optimal taul from version without tauf
%        helper=load(sprintf('OPT_notarget_2112_Sun%d_spillover0_knspil%d_taus0_noskill%d_notaul4_sep%d_extern0_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.Sun, indic.noknow_spill, indic.noskill, indic.sep, indic.xgrowth, indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
       helper=load(sprintf('OPT_notarget_plus%d_0501_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul4_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
               count, indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill,indic.sep, indic.extern, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa'))); 
   
end
         if  indic.tauf>1 && indic.tauf<5 % get better starting values!plotts.nsk
             if indic.tauf<2
                 T=11;
             end
             LF_SIM=helper.opt_all;       
                for ll=list.choice(1:end-1)
%                     fprintf('%s', ll)
                    x0LF(list.choice==ll)=LF_SIM(1, list.allvars==ll);
                end
                 if indic.tauf>=2
                     poll= [LF_SIM(:,list.allvars=='taul'), LF_SIM(:,list.allvars=='taus'),LF_SIM(:,list.allvars=='tauf'), LF_SIM(:,list.allvars=='lambdaa')];
                 elseif indic.tauf==0
                     poll= [LF_SIM(:,list.allvars=='taul'), LF_SIM(:,list.allvars=='taus'),zeros(size(LF_SIM(:,list.allvars=='tauf'))), LF_SIM(:,list.allvars=='lambdaa')];
                 elseif indic.tauf==1
                     poll= [zeros(size(LF_SIM(:,list.allvars=='taul'))), LF_SIM(:,list.allvars=='taus'),LF_SIM(:,list.allvars=='tauf'), LF_SIM(:,list.allvars=='lambdaa')];
                     
                 end
             if indic.xgrowth==0
                 [LF_SIM, pol, FVAL] = solve_LF_nows(T, list, poll, params, symms, x0LF, init201014, indexx, indic, Sall, Emss);
             else
                 [LF_SIM, pol, FVAL, indexx] = solve_LF_nows_xgrowth(T, list, poll, params, Sparams,  symms, x0LF, init201014, indexx, indic, Sall);
             end
             helper.LF_SIM=LF_SIM'; % but not correct policy!

             if indic.tauf<2
                helper.LF_SIM=[helper.LF_SIM;helper.LF_SIM(end,:)];
             end
                 
             for pp=list.pol(list.pol~='lambdaa')
                 helper.LF_SIM(:, list.allvars==pp)= helper.opt_all(:, list.allvars==pp);
             end
         else
           helper.LF_SIM=helper.opt_all;

         end
         T=12; 
        [LF_COUNT]=compequ(T, list, params, init201519, symms, helper.LF_SIM,indic, Emss, MOM);
        % helper.opt_all: as initial values and to deduce policy 
        if indic.oldCalib==1
           save(sprintf('COMPEquN_SIM_0501_taufopt%d_knspil%d_spillover%d_notaul%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat', indic.tauf,indic.noknow_spill,  indic.spillovers, indic.notaul, indic.noskill, indic.sep, indic.xgrowth, indic.PV, params(list.params=='etaa')),'LF_COUNT', 'Sparams');
        elseif indic.oldCalib==0
            if indic.target==1
                save(sprintf('COMPEquN_SIM_0501_taufopt%d_newCalib_target_emnet%d_Sun%d_knspil%d_spillover%d_notaul%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat', indic.tauf,indic.targetWhat, indic.Sun, indic.noknow_spill,  indic.spillovers, indic.notaul, indic.noskill, indic.sep, indic.xgrowth, indic.PV, params(list.params=='etaa')),'LF_COUNT', 'Sparams');
            else
               save(sprintf('COMPEquN_SIM_0501_taufopt%d_newCalib_notarget_emnet%d_Sun%d_knspil%d_spillover%d_notaul%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat', indic.tauf,indic.targetWhat, indic.Sun, indic.noknow_spill,  indic.spillovers, indic.notaul, indic.noskill, indic.sep, indic.xgrowth, indic.PV, params(list.params=='etaa')),'LF_COUNT', 'Sparams');
            end
        end
end
end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 6: PLOTS       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
etaa=params(list.params=='etaa');
weightext=0.01;
indic

% choose sort of plots to be plotted
plotts.regime_gov=  4; % = equals policy version to be plotted

plotts.table=       0;
plotts.cev  =       0; 
plotts.analyta =    0;
plotts.limit=       0; %==1 if plots emission target
plotts.robust=      0;
plotts.Ems_limit=   1; 
plotts.Ems_limit_Decomp=   0;
plotts.phi_sens =   0;
plotts.phiSens_newcalib_TvsNoT =0;
plotts.phi_effBauOPt_noBAu_newcalib =0;
plotts.phi_eff_newcalib     =0;
plotts.taulFixedtaufJoint_newcalib_polPer_Tauf=0;
plotts.sens_other=  0;
plotts.phi_newcalib        = 0;
plotts.phi_LF_newcalib     = 0;
plotts.phi_effLF_newcalib  = 0;
plotts.phi_effLFOPt_newcalib = 0;
plotts.phi_effBau_newcalib =0;
plotts.phi_effBauOPt_newcalib=0;
plotts.phi_newcalib_noeff  = 0;
plotts.comp_OPTPer_NCalib  = 0;
plotts.comp_OPT_NCAlib     = 0;
plotts.count_devs_both_NC  = 0;
plotts.phi_newcalib_TvsNoT = 0;
plotts.taulFixed_newcalib  = 0;
plotts.taulFixed_newcalib_pol  = 0;
plotts.taulFixed_newcalib_polPer = 0;
plotts.taulFixedtaufJoint_newcalib_polPer=0;
 plotts.single_pol_NC      =0;

plotts.countcomp=   0;
plotts.countcomp2=  0;
plotts.countcomp3=  0;
plotts.extern=      0;
plotts.compEff_mod_dev1         = 0;
plotts.count_taul_nsk_LF        = 0;
plotts.count_taul_xgr_LF        = 0;
plotts.count_taul_xgr_lev       = 0;
plotts.count_tauflev            = 0; % counterfactual with only tauf in laissez faire
plotts.count_taullev            = 0; % counterfactual with only taul in laissez faire
plotts.count_tauflev_Ben        = 0; % laissez faire, only optimal tauf and benchmark policy
plotts.count_tauflev_Ben_noLF   = 0;
plotts.compnsk_xgr              = 0;
plotts.compnsk_xgr1             = 0;

plotts.compnsk_xgr_dev          = 0;
plotts.compnsk_xgr_dev1         = 0;
plotts.count_modlev             = 0; 

plotts.count_modlev_eff         = 0;
plotts.single_pol               = 0;
plotts.singov                   = 0;

plotts.notaul                   = 0; % policy comparisons; this one needs to be switched on to get complete table
plotts.bau                      = 0; % do plot bau comparison
plotts.lf                       = 0; % comparison to laissez faire allocation 

plotts.comptarg                 = 0; % comparison with and without target
plotts.compeff                  = 0; % efficient versus optimal benchmark and non-benchmark
plotts.compeff3                 = 0; % sp versus optimal benchmark
plotts.compeff3_NC              = 0; % sp versus optimal benchmark

plotts.comp_LFOPT               = 0; % laissez faire and optimal with and without taul
plotts.compeff1=    0; %1; only social planner
plotts.compeff2=    0; %1; efficient and non benchmark
plotts.comp_OPT=    0; % laissez faire and optimal with and without taul
plotts.comp_OPT_NK= 0; % laissez faire and optimal with and without taul
plotts.comp_Bench_CountNK =0; % policy from model without knowledge spillovers in benchmark model
plotts.per_BAUt0 =  0;
plotts.per_effopt0= 0;
plotts.per_effoptd= 0;
plotts.per_baud =   0;
plotts.per_LFd  =   0; % dynamic lf as benchmark
plotts.per_LFd_NC = 0;
plotts.per_LFd_nt=  0; % dynamic lf as benchmark plus no income tax
plotts.per_LFd_ne_nt=0; % dynamic lf as benchmark plus no income tax

plotts.per_LFt0  =  0; % 2020  lf as benchmark
plotts.per_optd =   0;

plotts.tauf_comp=0;
plotts.compREd=0;

indic.noknow_spill=0;
indic.Sun=2;
indic.targetWhat=1;
for xgr =0
    for nsk=0
        for nknk=3
plotts.xgr = xgr; % main version to be used for plots
plotts.nsk = nsk;
plotts.sizeequ =0; % important for comparison of 
plotts.GOV =0;
plotts.nknk=nknk; % in the benchmark allocation there are kn spillovers
indic.count_techgap=0;
plotts
%%
%     plottsSP_PolRegimes(list, T, etaa, weightext,indic, params, Ems, plotts, percon);
plottsSP_tidiedUp(list, T-1, etaa, weightext,indic, params, Ems, plotts, percon, MOM); 
        end
    end
end
%% Optimal policy results
etaa=params(list.params=='etaa');
weightext=0.01;
indic

% choose sort of plots to be plotted
plotts.ems =        1;
plotts.ems_goals =  0;


plotts.table=       0;
plotts.cev  =       0; 
plotts.analyta =    0;
plotts.limit=       0; %==1 if plots emission target
plotts.robust=      0;
plotts.countcomp=   0;
plotts.countcomp2=  0;
plotts.countcomp3=  0;
plotts.extern=      0;
plotts.compEff_mod_dev1         = 0;
plotts.count_taul_nsk_LF        = 0;
plotts.count_taul_xgr_LF        = 0;
plotts.count_taul_xgr_lev       = 0;
plotts.count_tauflev            = 0; % counterfactual with only tauf in laissez faire
plotts.count_taullev            = 0; % counterfactual with only taul in laissez faire
plotts.count_tauflev_Ben        = 0; % laissez faire, only optimal tauf and benchmark policy
plotts.count_tauflev_Ben_noLF   = 0;
plotts.compnsk_xgr              = 0;
plotts.compnsk_xgr1             = 0;

plotts.compnsk_xgr_dev          = 0;
plotts.compnsk_xgr_dev1         = 0;
plotts.count_modlev             = 0; 
plotts.count_devs               = 0;
plotts.count_devs_fromcto       = 0;
plotts.count_devs_both          = 0;
plotts.count_modlev_eff         = 0;
plotts.single_pol               = 0;     
plotts.singov                   = 0;

plotts.notaul                   = 0; % policy comparisons; this one needs to be switched on to get complete table
plotts.bau                      = 0; % do plot bau comparison
plotts.lf                       = 0; % comparison to laissez faire allocation in levels

plotts.comptarg                 = 0; % comparison with and without target
plotts.compeff                  = 0; % efficient versus optimal benchmark and non-benchmark
plotts.compeff3                 = 0; % sp versus optimal benchmark
plotts.compeff4                 = 0; % sp versus optimal benchmark
plotts.comp_LFOPT               = 0; % laissez faire and optimal with and without taul
plotts.compeff1=    0; %1; only social planner
plotts.compeff2=    0; %1; efficient and non benchmark
plotts.comp_OPT=    0; % laissez faire and optimal with and without taul
plotts.comp_OPTPer= 0; % comparison in percent with and without taul
plotts.comp_OPT_NK= 0; % laissez faire and optimal with and without taul
plotts.comp_Bench_CountNK =0; % policy from model without knowledge spillovers in benchmark model
plotts.per_BAUt0 =  0;
plotts.per_effopt0= 0;
plotts.per_effoptd= 0;
plotts.per_baud =   0;
plotts.per_LFd  =   0; % dynamic lf as benchmark
plotts.per_LFd_nt=  0; % dynamic lf as benchmark plus no income tax
plotts.per_LFd_ne_nt=0; % dynamic lf as benchmark plus no income tax

plotts.per_LFt0  =  0; % 2020  lf as benchmark
plotts.per_optd  =  0;

plotts.tauf_comp=0;
plotts.compREd=0;
for rr= [4]
    plotts.regime_gov=  rr; % = equals policy version to be plotted

for xgr =0
    for nsk=0
        for nknk=0
            T=12;
plotts.xgr = xgr; % main version to be used for plots
plotts.nsk = nsk;
plotts.sizeequ =0; % important for comparison of 
plotts.GOV =0;
plotts.extern =0;
indic.noknow_spill=nknk; % in the benchmark allocation there are kn spillovers
indic.slides =1;
plotts
%%
%     plottsSP_PolRegimes(list, T, etaa, weightext,indic, params, Ems, plotts, percon);

plotts_extT(list, T, etaa, weightext,indic, params, Ems, plotts, percon, MOM, StatsEms)
        end
    end
end
end
%% tables
% constructed in plotsSP file

tt=load('Table_SWF_July22_sep1_noskill0_etaa0.79_xgrowth0_extern0.mat');
addpath('tools')
table2latex(tt.TableSWF_PV)
%% Analytical model files
% solves analytical model and compares policy regimes
solve_easy;

%% Simulate economy after last optimisation period using competitive economy files
% under the assumption of constant policy
% GOAL: check if emission target gets violated

indic.notaul =3;

indic.noskill=0;
indic.taus   =0;
indic.sep    =1;
indic.spillovers = 0;

for xgr =0
    indic.xgrowth=xgr;
    res_opt=load(sprintf('OPT_target_spillover%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_etaa%.2f.mat', indic.spillovers, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.xgrowth, params(list.params=='etaa')), 'opt_all', 'Sparams');
    % simulate in laissez faire for peroiods after 2075 with assumption of
    % constant policy (after 2080 because laissez faire code starts from first period)

    %read in policy
    taul = res_opt.opt_all(T, list.sepallvars=='taul');
    tauf = res_opt.opt_all(T, list.sepallvars=='tauf');
    taus = res_opt.opt_all(T, list.sepallvars=='taus');
    lambdaa = res_opt.opt_all(T, list.sepallvars=='lambdaa');
    polLFcon = eval(symms.pol);
    %- initial 2080 technology levels
    % => first simulation period is 2085
    Af0 = res_opt.opt_all(T, list.sepallvars=='Af');
    Ag0 = res_opt.opt_all(T, list.sepallvars=='Ag');
    An0 = res_opt.opt_all(T, list.sepallvars=='An');
    initcon=eval(symms.init);

    %- initial values as in x0LF; indicator: list.choice
    x0LFcon=x0LF;
    for ii =list.choice
        if ismember(ii, ["ws" "gammas"])
            x0LFcon(list.choice==ii) = res_opt.opt_all(T, list.sepallvars==sprintf('%sg', ii ));
        else
            x0LFcon(list.choice==ii) = res_opt.opt_all(T, list.sepallvars==ii);
        end
    end
    %- choose periods for which to simulate
    Tcon = 12; 

    if indic.xgrowth==0
        [LF_SIM, polLF, FVAL] =solve_LF_nows_continued(Tcon, list, polLFcon, params, Sparams,  symms, x0LFcon, initcon, indexx, indic, Sall);
        helper.LF_SIM=LF_SIM;
        indic.xgrowth=0;
        [LF_SIM]=solve_LF_VECT(T, list, params,symms, initcon, helper, indic);
    else
        [LF_SIM, polLF, FVAL, indexx] = solve_LF_nows_xgrowth_continued(Tcon, list, polLFcon, params, Sparams,  symms, x0LFcon, initcon, indexx, indic, Sall);
         LF_SIM = LF_SIM';   
    end

    save(sprintf('LF_CON_spillover%d_noskill%d_sep%d_xgrowth%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.xgrowth, Sparams.etaa), 'Tcon', 'LF_SIM', 'Sparams', 'polLFcon');
% evaluate how economy evolves
Fcon= LF_SIM(:,list.sepallvars=='F');
Gcon= LF_SIM(:,list.sepallvars=='G');
GFFcon = Gcon./Fcon;
Ftargetcon = (Ems(end)+Sparams.deltaa)/Sparams.omegaa;
    if sum(Fcon >Ftargetcon)
        fprintf('The emission target gets violated in at least one period after intervention; exog growth %d', indic.xgrowth);
    end    
end

concon= load(sprintf('LF_CON_spillover%d_noskill%d_sep%d_xgrowth%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.xgrowth, Sparams.etaa),'LF_SIM');
%=> given the model's assumptions there needs to be continued intervention
%   to meet the emission limit
%%
for BN=1
    indic.BN=BN;
    for inn=0:1
        indic.ineq=inn;
        indic
        plottsSP(list, T, etaa, weightext,indic, params, Ems, plotts);
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% find good initial starting point for target opt %%
%%% as a function of taul
taul= 0.5;
pf= 1;
exx = polExp(pf, params, polCALIB, list, taul, T, Ems, indexx);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Symbolic approach to solve Ramsey problem %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indic.target=1;
indic.noskill=0;
% indic.sep=0;
indic.target=1; 
%- with etaa=1
% params(list.params=='etaa')=1;
% Sparams.etaa=params(list.params=='etaa');

%1) get objective function 
if indic.sep==1
    [OB_RAM, list, symms, Ftarget]= model_ram_sep( list, params, T, init201519, indic, Ems, symms);
else
    [OB_RAM, list, symms, Ftarget]= model_ram( list, params, T, init201519, indic, Ems, symms);
end
%- x is a symbolic vector of choice variables! 

%2) take derivatives and write resulting equations as function
if indic.target==1
    [indexx, model, list]=symmodel_eq(OB_RAM, symms.optALL, params,  Ftarget, 'Ram_Model_target_1905_KTS', list, indic, indexx);
else
    if indic.sep==0
        [indexx, model]=symmodel_eq(OB_RAM, symms.optALL, params,  Ftarget, 'Ram_Model_notarget_1905_KTS', list, indic, indexx);
    else
        [indexx, model]=symmodel_eq_sep(OB_RAM, symms.optALL, params,  Ftarget, 'Ram_Model_notarget_sep_1905', list, indic, indexx);
    end
end

%3) solve model using fsolve
if indic.sep==1
    RAM = solve_sym_sep(symms, list, Ftarget, indic);
else
    RAM = solve_sym(symms, list, Ftarget, indic);
end
