function []=plotts_fun(list, T, etaa, weightext,indic, params, Ems, plotts, percon, MOM)

% this script plots results

date="30Aug23";
if ~isfile(sprintf('../figures/all_%s', date ))
    mkdir(sprintf('../figures/all_%s', date));
end

%- color for graphs
orrange= [0.8500 0.3250 0.0980];
grrey = [0.6 0.6 0.6];

%- variables
syms H Y F E N Emnet G pg pn pf pee tauf taul taus taurese tauresg Trans w ws C Lg Lf Ln xn xg xf sn sff sg SWF Af Ag An A real
syms analyTaul PV CEVv CEVvPV CEVvDy Tauf AgAf sffsg sgsff GFF EY CY LgLf pgpftf pepn gAg gAf gAn gAagg Utilcon Utillab real
symms.plotsvarsProd =[Y N E G F];
symms.plotsvarsHH =[H C SWF Emnet]; 
symms.plotsvarsRes =[sn sff sg Af Ag An A];  
symms.plotsvarsProdIn =[xn xg xf Ln Lg Lf];  

symms.plotsvarsPri =[pg pn pf pee w ws];  
symms.plotsvarsPol =[taus tauf taul Trans];  
symms.plotsvarsAdd = [analyTaul PV Tauf AgAf sffsg sgsff GFF EY CY LgLf pgpftf pepn gAagg gAg gAf gAn Utilcon Utillab ];
% already exists: symms.addgov
symms.comp=[ CEVv CEVvDy CEVvPV ]; % for comparison of policy interventions, 


listt.plotsvarsProd=string(symms.plotsvarsProd);
listt.plotsvarsProdIn=string(symms.plotsvarsProdIn);
listt.plotsvarsHH=string(symms.plotsvarsHH);
listt.plotsvarsRes=string(symms.plotsvarsRes);
listt.plotsvarsPol=string(symms.plotsvarsPol);
listt.plotsvarsPri=string(symms.plotsvarsPri);
listt.plotsvarsAdd=string(symms.plotsvarsAdd);
list.comp=string(symms.comp);
list.growthrates =string([gAg gAf gAn gAagg]);

lisst = containers.Map({'Prod', 'ProdIn','Res', 'HH', 'Pol', 'Pri', 'Add'}, {listt.plotsvarsProd, listt.plotsvarsProdIn, ...
    listt.plotsvarsRes,listt.plotsvarsHH,listt.plotsvarsPol, listt.plotsvarsPri, listt.plotsvarsAdd});
 
%- variables which to plot in one graph plus legend
lissComp = containers.Map({'Labour', 'Science', 'Growth', 'LabourInp'}, {string([H]), string([sff sg sn]), string([gAf gAg  gAn gAagg]), string([Lf Lg])});
legg = containers.Map({'Labour', 'Science', 'Growth', 'LabourInp'}, {["Hours"], ["fossil", "green", "non-energy", "total"], ["fossil", "green", "non-energy", "aggregate"], ["fossil", "green"]});

%- new list 
% varlist_polcomp=[list.allvars, list.addgov, string(symms.plotsvarsAdd)];

varlist=[list.allvars, string(symms.plotsvarsAdd)];
%% read in results

%-initialise cells to store matrices: one for each combination of xgr and nsk
OTHERPOL={}; % cell to save containers of different policy versions
% OTHERPOL_xgr={}; 
% OTHERPOL_nsk={};
% OTHERPOL_xgr_nsk={}; 

% baseline results 

%- sp solution independent of policy

helper= load(sprintf('../output/SP_target_2407_emnet%d_phii%d_sigma%d_Bop%d_util%d_PV%d.mat',...
            indic.targetWhat, indic.know_spill, indic.sigmaWorker, indic.Bop, indic.util), 'sp_all', 'obs');
sp_t=helper.sp_all';

helper=load(sprintf('../output/SP_notarget_2407_phii%d_sigma%d_Bop%d_util%d_PV%d.mat',...
             indic.know_spill, indic.sigmaWorker, indic.Bop, indic.util), 'sp_all', 'obs');        
sp_not=helper.sp_all';

%- other results
for i=[0,1,2,3] % loop over policy versions
  
    helper = load(sprintf('../output/LF_phii%d_sigma%d_Bop%d_util%d.mat', ...
              indic.know_spill, indic.sigmaWorker, indic.Bop));
    LF = helper.COMP';
    
    helper=load(sprintf('../output/OPT_notarget_2407_notaul%d_subsres%d_phii%d_sigma%d_Bop%d_util%d_PV%d.mat',...
            i, indic.subsres, indic.know_spill, indic.sigmaWorker, indic.Bop, indic.util), 'opt_all', 'obs');
    opt_not=helper.opt_all';
    
    helper=load(sprintf('../output/OPT_target_2407_emnet%d_notaul%d_subsres%d_phii%d_sigma%d_Bop%d_util%d_PV%d.mat',...
        indic.targetWhat, i, indic.subsres, indic.know_spill, indic.sigmaWorker, indic.Bop, indic.util), 'opt_all', 'obs');

    opt_t=helper.opt_all';

    RES = containers.Map({'LF', 'SP_T', 'SP_NOT' ,'OPT_T', 'OPT_NOT'},...
                            { LF, sp_t, sp_not, opt_t, opt_not});
    %- add additional variables
    %if xgr==0 && nsk==0
    OTHERPOL{i+1}=add_vars(RES, list, params, indic, list.allvars, symms, MOM);
    % elseif xgr==0 && nsk==1
    %     OTHERPOL_nsk{i+1}=add_vars(RES, list, params, indic, list.allvars, symms, MOM);
    % elseif xgr==1 && nsk==0
    %     OTHERPOL_xgr{i+1}=add_vars(RES, list, params, indic, list.allvars, symms, MOM);
    % elseif xgr==1 && nsk==1
    %     OTHERPOL_xgr_nsk{i+1}=add_vars(RES, list, params, indic, list.allvars, symms, MOM);
    % end
end

%% counterfac technology

% helper=load('SP_target_1008_countec_spillover0_knspil0_noskill0_sep0_xgrowth0_zero1_PV0_sizeequ7.900000e-01_etaa.mat');
% sp_t=helper.sp_all';
% helper=load('SP_notarget_1008_countec_spillover0_knspil0_noskill0_sep0_extern0_xgrowth0_PV1_sizeequ0_etaa0.79.mat');
% sp_not=helper.sp_all';
% RES_count=containers.Map({'SP_T', 'SP_NOT' },...
%                                 {  sp_t, sp_not});
% RES_count=add_vars(RES_count, list, params, indic, list.allvars, symms, MOM);

%% counterfac emission limit
% for ii=[4,5]
%      helper=load(sprintf('OPT_target_0509_emnet1_spillover%d_knspil%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat', indic.spillovers,plotts.nknk, plotts.nsk, ii, indic.sep, plotts.xgr,indic.PV, plotts.sizeequ, plotts.GOV, etaa));
%      opt_t_notaus=helper.opt_all';
%      RES_Ems=containers.Map({'OPT_T_NoTaus', },{opt_t_notaus});
%      RES_Ems=add_vars(RES_Ems, list, params, indic, list.allvars, symms, MOM);
%      OTHERPOL_Ems{ii+1}=RES_Ems;
% end


%% Tables
if plotts.table==1 % 
  %- choose respective containers cell
     OTHERPOLL= OTHERPOL;
  %- discount vector
      betaa=params(list.params=='betaa');

      % version with T=11
      disc=repmat(betaa, 1,T);
      expp=0:T-1;
      vec_discount= disc.^expp;

     %- Table
     TableSWF_PV=table(keys(RES)',zeros(length(keys(RES)),1), zeros(length(keys(RES)),1),zeros(length(keys(RES)),1),zeros(length(keys(RES)),1));
     TableSWF_PV.Properties.VariableNames={'Allocation','lump-sum trans with taul', 'lump-sum trans no taul',...
                                            'green subsidies with taul', 'green subsidies no taul'}; 
                                            % row names equivalent/ordered
                                            % as in indic.notaul numeration

    %- all results
    for i =keys(RES) % rows: policy versions
                     % columns: equilibrium (LF, Optimal, Social Planner with and without target)
         ii=string(i);
         count=0; % to keep track of which container, i.e. policy,  is used
        for ccc=OTHERPOLL
            pp=ccc{1}; % to get content of cell (allocation under a certain policy
            count=count+1; % policies
            allvars=pp(ii);
            TableSWF_PV{TableSWF_PV.Allocation==ii,count+1}=vec_discount*allvars(find(varlist=='SWF'),:)'+indic.PV*allvars(find(varlist=='PV'),1);
        end
    end
    %%
    %     OTHERPOL_nsk{i+1}=add_vars(RES, list, params, indic, list.allvars, symms, MOM);
    % elseif xgr==1 && nsk==0
    %     OTHERPOL_xgr{i+1}=add_vars(RES, list, params, indic, list.allvars, symms, MOM);
    % elseif xgr==1 && nsk==1
    %     OTHERPOL_xgr_nsk{i+1}=add_vars(RES, list, params, indic, list.allvars, symms, MOM);
    % end
     save(sprintf('../figures/all_%s/Table_SWF_subsres%d_phii%d_sigma%d_Bop%d.mat',...
         date, indic.subsres,indic.know_spill, indic.sigmaWorker, indic.Bop), 'TableSWF_PV');
end


%% Pick main policy version for plots
%if plotts.xgr ==0 && plotts.nsk==0
OTHERPOLL= OTHERPOL;
% elseif plotts.xgr ==1 && plotts.nsk==0
%     OTHERPOLL= OTHERPOL_xgr;
% elseif plotts.xgr ==0 && plotts.nsk==1
%     OTHERPOLL= OTHERPOL_nsk;
% elseif plotts.xgr ==1 && plotts.nsk==1
%     OTHERPOLL= OTHERPOL_xgr_nsk;
% end

%% table CEV
if plotts.cev==1
    %- calculate CEV for a pair of policy regimes each
    for gg =[0,2]
        plotts.regime_gov=gg;
    if plotts.regime_gov==0 % lump-sum transfers
        h1= OTHERPOLL{1}; % taul can be used
        h2= OTHERPOLL{2}; % taul cannot be used
    elseif plotts.regime_gov==2 % green subsidies
        h1= OTHERPOLL{3}; % taul can be used
        h2= OTHERPOLL{4}; % taul cannot be used
    end
        
%     [COMP, COMPTable] = comp_CEV(h1('OPT_T_NoTaus'),h2('OPT_T_NoTaus') , varlist, varlist, symms, list, params, T, indic);   
    [COMP, COMPTable] = comp_CEV(h1('OPT_T'),h2('OPT_T') , varlist, varlist, symms, list, params, T, indic);   

    save(sprintf('Table_CEV_%s_regime%d_emnet%d_subsres%d_phii%d_sigma%d_Bop%d.mat',...
        date, plotts.regime_gov,indic.subsres,indic.know_spill, indic.sigmaWorker, indic.Bop), 'COMPTable');

    end
end
%% Plots
%- axes
time = 1:T;
txx=1:2:T; % reducing indices
%- using start year of beginning period 
Year =transpose(year(['2020'; '2025';'2030'; '2035';'2040'; '2045';'2050'; '2055'; '2060';'2065';'2070';'2075'],'yyyy'));
Year10 =transpose(year(['2020';'2030'; '2040'; '2050';'2060';'2070'],'yyyy'));

%'LF', 'BAU', 'SP_T', 'SP_NoT' 'OPT_T_NOTaul', 'OPT_T_WithTaul',
%'OPT_NoT_NOTaul', 'OPT_NoT_WithTaul', 'Count_TaulFixed' count_tauf_TaulFixed_TaufJoint


%% All figures single
if plotts.single_pol_NC==1
    
    fprintf('plotting single graphs')

    %- loop over economy versions
    for nnt=[ "base"]
        if nnt =="base"
            RES=RES_NCalib;
        else
            RES= RES_NCalib_RS;
        end
    for s=["T", "NoT"]
        ss=string(s);
        wt=RES(sprintf("OPT_%s_WithTaul", ss));
        
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
            if ll=="HH" && varr=="Emnet"
                main=plot(time,wt(find(varlist==varr),1:T),time(percon+1:end),Ems(1:T), 'LineWidth', 1.1);  
                set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
                lgd=legend('net emissions' , 'net emission limit',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');

            else
                main=plot(time,wt(find(varlist==varr),1:T), 'LineWidth', 1.1);   
                set(main, {'LineStyle'},{'-'}, {'color'}, {'k'} )   

            end
           xticks(txx)
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)

           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            if varr=="dTaulAv" || varr=="dTaulAvS"
                ytickformat('%.1f')
            elseif varr =="taul"
                ytickformat('%.3f')
            end
            xticklabels(Year10)
            path=sprintf('figures/all_%s/Single_NC_%s_%s_emnet%d_Sun%d_regime%s_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_extern%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.png',...
               date, ss, varr, indic.targetWhat, indic.Sun, nnt, indic.spillovers,plotts.nknk, plotts.nsk, indic.sep,plotts.xgr,indic.extern, indic.PV, plotts.sizeequ, plotts.GOV, etaa);
         
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end

%% levels LF with and without target 
if plotts.levels_LF_all==1
    fprintf('plott levels with lf')

    for poll=[2]
    RES    = OTHERPOLL{poll+1}; % taul available
    RES_nt = OTHERPOLL{poll+2}; % no taul available
    LFall= RES("LF");
    for s= "T" %["T", "NOT"]
        ss=string(s);
        nt= RES_nt(sprintf("OPT_%s", ss));
        wt=RES(sprintf("OPT_%s", ss));
        eff=RES(sprintf("SP_%s",ss) );
            
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(LFall(find(varlist==varr),1:T)), time,(nt(find(varlist==varr),1:T)), time,(wt(find(varlist==varr),1:T)), time,(eff(find(varlist==varr),1:T)), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{':'; '--';'-'; '--'}, {'color'}, {'k';'k';'k';grrey} )   
            if lgdind==1
               lgd=legend('laissez-faire','no labor tax', 'benchmark', 'first best',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           % if ismember(varr, list.growthrates)
           %      xlim([1, time(end-1)])
           % else             
           xlim([1, time(end)])
           % end
           
            ax=gca;
            ax.FontSize=13;
            if varr=="Tauf" || varr=="GFF"
                ytickformat('%.0f')
            elseif varr =="H"
                ytickformat('%.3f')
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            path=sprintf('../figures/all_%s/Levels_eff2pol_LF_pol%d_%s_%s_emnet%d_subsres%d_knspil%d_sigma%d_Bop%d_util%d_lgd%d.png',...
                date,poll, ss, varr, indic.targetWhat, indic.subsres, indic.know_spill,indic.sigmaWorker, indic.Bop, indic.util, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
        end
    end
    end
end



%% plotts with and without taul
if plotts.pol_comp==1
    fprintf('plott with and without taul')

    for poll=[0] % lump sum trans or subsidies
        RES    = OTHERPOLL{poll+1}; % taul available
        RES_nt = OTHERPOLL{poll+2}; % no taul available
      for s =["NOT", "T"]
        ss=string(s);

        nt = RES_nt(sprintf("OPT_%s",ss));
        wt = RES(sprintf("OPT_%s", ss));
        
    for lgdind=0:1
    for l ="Pol" %keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(nt(find(varlist==varr),1:T)), time,(wt(find(varlist==varr),1:T)),  'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'--';'-'}, {'color'}, {'k';'k'} )   
            if lgdind==1
               lgd=legend('no labor tax', 'benchmark', 'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            if varr=="Tauf"
                ytickformat('%.0f')
            elseif varr =="H"
                ytickformat('%.3f')
            elseif  varr=="taul"
                ytickformat('%.2f')
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            path=sprintf('../figures/all_%s/Optimal_comp_pol%d_%s_%s_emnet%d_subsres%d_knspil%d_sigma%d_Bop%d_util%d_lgd%d.png',...
                date,poll, ss, varr, indic.targetWhat, indic.subsres, indic.know_spill,indic.sigmaWorker, indic.Bop, indic.util, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
    end
    end
    end
end

%% plotts with and without taul
if plotts.pol_regime_comp==1
    fprintf('plott taul with lump-sum trans vs. subsidies')

        RES_lt    = OTHERPOLL{1}; % lump-sum trans with taul
        RES_subs  = OTHERPOLL{3}; % subs with taul
        RES_subs_nt  = OTHERPOLL{4}; % no taul available

        for s =[ "T", "NOT"]
        ss=string(s);

        lt = RES_lt(sprintf("OPT_%s",ss));
        subs = RES_subs(sprintf("OPT_%s", ss));
        subs_nt = RES_subs_nt(sprintf("OPT_%s", ss));

    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(lt(find(varlist==varr),1:T)), time,(subs(find(varlist==varr),1:T)),time,(subs_nt(find(varlist==varr),1:T)),  'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'--';'-'; '--'}, {'color'}, {'k';'k'; 'b'} )   
            if lgdind==1
               lgd=legend('lump-sum transfers', 'green subsidies', 'green subsidies, no labor tax',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            if varr=="Tauf"
                ytickformat('%.0f')
            elseif varr =="dTaulAv"
                ytickformat('%.1f')
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            path=sprintf('../figures/all_%s/Comp_LSTrans_Subs_%s_%s_emnet%d_subsres%d_knspil%d_sigma%d_Bop%d_util%d_lgd%d.png',...
                date, ss, varr, indic.targetWhat, indic.subsres, indic.know_spill,indic.sigmaWorker, indic.Bop, indic.util, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
    end
    end
end

%% no efficient
if plotts.phi_newcalib_TvsNoT==1
    fprintf('plott new calib no T vs NoT')

for nnt=["base"]
    nnt=string(nnt);
    if nnt=="base"
        RES=RES_NCalib;
    else
        RES=RES_NCalib_RS;
    end
    wt = RES(sprintf("OPT_T_WithTaul"));
    not = RES(sprintf("OPT_NoT_WithTaul"));

    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(not(find(varlist==varr),1:T)), time,(wt(find(varlist==varr),1:T)),  'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'--';'-'}, {'color'}, {'b';'k'} )   
            
           xticks(txx)
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)

           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            if lgdind==1
               lgd=legend('no emission limit', 'with emission limit', 'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
           
            ax=gca;
            ax.FontSize=13;
            if varr=="Tauf"
                ytickformat('%.0f')
            elseif varr =="dTaulAv"
                ytickformat('%.1f')
            elseif varr =="taul"
                ytickformat('%.3f')
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            path=sprintf('figures/all_%s/NewCalib_pol_TvsNoT_%s_%s_emnet%d_Sun%d_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date,varr,nnt, indic.targetWhat, indic.Sun, indic.spillovers, plotts.nknk, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end

%% no efficient
if plotts.phiSens_newcalib_TvsNoT==1

    fprintf('plott new calib no T vs NoT')
        nknk0= RES_Sens("nknk0_T");
        nknk05= RES_Sens("nknk05_T");
        nknk025= RES_Sens("nknk025_T");
        nknk075=RES_Sens("nknk075_T"); %RES_NCalib("OPT_T_WithTaul");
        nknk0_nt= RES_Sens("nknk0_NoT");
        nknk05_nt= RES_Sens("nknk05_NoT");
        nknk025_nt= RES_Sens("nknk025_NoT");
        nknk075_nt=RES_Sens("nknk075_NoT"); %RES_NCalib("OPT_T_WithTaul");

    perdif0=100*(nknk0-nknk0_nt)./nknk0_nt;
    perdif025=100*(nknk025-nknk025_nt)./nknk025_nt;
    perdif05=100*(nknk05-nknk05_nt)./nknk05_nt;
    perdif075=100*(nknk075-nknk075_nt)./nknk075_nt;

    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(perdif075(find(varlist==varr),1:T)), time,(perdif05(find(varlist==varr),1:T)),...
                 time,(perdif025(find(varlist==varr),1:T)),  time,(perdif0(find(varlist==varr),1:T)),  time,zeros(size(perdif0(find(varlist==varr),1:T))), 'LineWidth', 1.1);   
             set(main, {'LineStyle'},{'--';'-'; '--'; ':'; '--'}, {'color'}, {'b';'k';orrange; grrey; grrey} )   
       
            
           xticks(txx)
           xlim([1, time(end-1)])
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)
           
           if lgdind==1
               lgd=legend('$\phi=0.75$', '$\phi=0.5$', '$\phi=0.25$','$\phi=0$',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            ax=gca;
            ax.FontSize=13;
            if varr=="Tauf"
                ytickformat('%.0f')
            elseif varr =="dTaulAv"
                ytickformat('%.1f')
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            path=sprintf('figures/all_%s/NewCalib_Sens_TvsNoT_%s_emnet%d_Sun%d_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date,varr, indic.targetWhat, indic.Sun, indic.spillovers, plotts.nknk, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end

%% comparison only carbon tax (reg 5) to only tauf from combined
if plotts.count_devs_both_NC==1
    
    fprintf('plott counterfactual deviation from carbon tax only pol => role of adjustment in tauf')

    %- read in variable container of chosen regime
         for s=["NoT"]
             ss=string(s);
        allvarsntaul= RES_NCalib(sprintf("OPT_%s_NOTaul", ss));
        allvarscount=RES_count_NCalib(sprintf("CountOnlyTauf_%s",ss)); % version with only tauf
        perdiftauf= 100*(allvarscount-allvarsntaul)./allvarsntaul;
        
         
        allvars= RES_NCalib(sprintf("OPT_%s_WithTaul",ss));
        %allvarscount=RES_count("CountOnlyTauf"); % version with only tauf
        perdiftaul= 100*(allvars-allvarscount)./allvarsntaul; % relative to carbon-tax-only policy
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

             main=plot(time,perdiftauf(find(varlist==varr),1:T), time,perdiftaul(find(varlist==varr),1:T), time,zeros(size(perdiftauf(find(varlist==varr),1:T))), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-'; '--'; '--'}, {'color'}, {'k';orrange; grrey} )   
            
            
           xticks(txx)
           xlim([1, time(end-1)])
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)

            ax=gca;
            ax.FontSize=13;
%             if varr== "sgsff"
%                 ytickformat('%.1f')
%             elseif varr == "Hagg"
%                 ytickformat('%.2f')
%                 ylim([-0.12, 0.061])
%             else
%                  ytickformat('%.2f')
%             end
            
            if lgdind==1
               lgd=legend('effect $\tau_F$' , 'effect $\tau_\iota$', 'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');
            end
            xticklabels(Year10)

            path=sprintf('figures/all_%s/CountTAUF_Both_Opt_NewCalib_%s_emnet%d_Sun%d_target_%s_nsk%d_xgr%d_knspil%d_regime%d_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date,ss, indic.targetWhat, indic.Sun, varr ,plotts.nsk, plotts.xgr, plotts.nknk, plotts.regime_gov, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
    end
    end
end

%% optimal with and without taul in levels 
if plotts.comp_OPTPer_NCalib==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting percentage opt and no taul no eff') 
    
   for s=["T", "NoT"]
       ss=string(s);
        allvarsnt= RES_NCalib(sprintf("OPT_%s_NOTaul", ss));
        allvars=RES_NCalib(sprintf("OPT_%s_WithTaul", ss));
        Perdif=100*(allvars-allvarsnt)./allvarsnt;
    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot( time,(Perdif(find(varlist==varr),1:T)),time,zeros(size(time)));            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; grrey},{'LineWidth'}, {1.1; 1} )   
           xticks(txx)
           xlim([1, time(end-1)])
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)
           
            ax=gca;
            ax.FontSize=13;
%             if varr=="SWF" || varr== "sn"
%                ytickformat('%.0f')
%             elseif varr=="sff" || varr=="sg" ||  varr == "sgsff" || varr =="GFF" || varr =="sffsg"
%                ytickformat('%.1f')
%             else
%                ytickformat('%.2f')
%             end
            xticklabels(Year10)
%            if lgdind==1
%               lgd=legend('with income tax', 'without income tax', 'Interpreter', 'latex');
%               set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
%            end
        path=sprintf('figures/all_%s/%s_COMPtaulPerNewCalib_%s_regime%d_emnet%d_Sun%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.png',date, varr,ss, plotts.regime_gov,indic.targetWhat, indic.Sun, indic.spillovers, plotts.nknk, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa);
        
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
    end
end

%% comparison parameter alternatives
if plotts.phi_sens==1
    
    fprintf('plott sensitivity phi')

        nknk0= RES_Sens("nknk0_T");
        nknk05= RES_Sens("nknk05_T");
        nknk025= RES_Sens("nknk025_T");
        nknk075=RES_Sens("nknk075_T"); %RES_NCalib("OPT_T_WithTaul");

    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(nknk075(find(varlist==varr),1:T)), time,(nknk05(find(varlist==varr),1:T)), time,(nknk025(find(varlist==varr),1:T)),time,nknk0(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'--';'-'; '--'; ':'}, {'color'}, {'b';'k';orrange; grrey} )   
            if lgdind==1
               lgd=legend('$\phi=0.75$', '$\phi=0.5$', '$\phi=0.25$','$\phi=0$',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            if varr=="Tauf"
                ytickformat('%.0f')
            elseif varr =="dTaulAv"
                ytickformat('%.1f')
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            path=sprintf('figures/all_%s/Phi_SensN_%s_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr, indic.spillovers,indic.noknow_spill, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end
%% comparison parameter alternatives
if plotts.sens_other==1
    
    fprintf('plott sensitivity other')

        RES=OTHERPOL{plotts.regime_gov+1}; 
        ben= RES("OPT_T_NoTaus");
        sig=RES_sigma("OPT_T_NoTaus");
        bop=RES_bop("OPT_T_NoTaus");
        ela=RES_ela("OPT_T_NoTaus");
        tec = RES_tec("OPT_T_NoTaus");
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(ben(find(varlist==varr),1:T)), time,(sig(find(varlist==varr),1:T)), ...
                time,(bop(find(varlist==varr),1:T)),time,ela(find(varlist==varr),1:T),  'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'; '--'; '--'}, {'color'}, {'k';'b';orrange; grrey} )   
            if lgdind==1
               lgd=legend('baseline', '$\sigma=2/3$', '$\theta=2$','$\varepsilon_e=10$',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           xlim([1, time(end-1)])

           
            ax=gca;
            ax.FontSize=13;
            if varr=="Tauf"
                ytickformat('%.0f')
            elseif varr =="dTaulAv"
                ytickformat('%.0f')
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            path=sprintf('figures/all_%s/SensOther_%s_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr, indic.spillovers,indic.noknow_spill, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end
%% effect of taul in model without xgr and with xgr
% compare effect of only taul in benchmark and in exogenous growth model
if plotts.count_taul_nsk_LF==1
    
    fprintf('plott only taul nsk model devLF')

        RES=OTHERPOL{plotts.regime_gov+1}; 
        RESnsk=OTHERPOL_nsk{plotts.regime_gov+1}; 

        allvarslf= RES("LF");
        allvarstlf= RESnsk("LF");
        allvars= RES_count_SAMETAUL_onlytaul("test");
        allvarst=RES_count_SAMETAUL_onlytaul('nsk');
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(allvars(find(varlist==varr),1:T)), time,allvarslf(find(varlist==varr),1:T),...
                time,allvarst(find(varlist==varr),1:T), time, allvarstlf(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';':'; '--'; ':'}, {'color'}, {'k'; orrange; 'b';grrey} )   
            if lgdind==1
               lgd=legend('benchmark' , 'benchmark, LF', 'homogeneous skill', 'homogeneous skill, LF',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/CountNskTaulLF_target_%s_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end

%% effect of taul in model without xgr and with xgr
% compare effect of only taul in benchmark and in exogenous growth model
if plotts.count_taul_xgr_LF==1
    
    fprintf('plott only taul xgr model devLF')

        RES=OTHERPOL{plotts.regime_gov+1}; 
        RESxgr=OTHERPOL_xgr{plotts.regime_gov+1}; 

        allvarslf= RES("LF");
        allvarstlf= RESxgr("LF");
        allvars= RES_count_SAMETAUL_onlytaul("test");
        allvarst=RES_count_SAMETAUL_onlytaul('xgr');
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(allvars(find(varlist==varr),1:T)), time,allvarslf(find(varlist==varr),1:T),...
                time,allvarst(find(varlist==varr),1:T), time, allvarstlf(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';':'; '--'; '--'}, {'color'}, {'k';orrange; 'b';grrey} )   
            if lgdind==1
               lgd=legend('benchmark' , 'benchmark, LF', 'exogenous growth', 'exogenous growth, LF',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/CountXgrTaulLF_target_%s_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end


%% Comparison Emission limits
if plotts.Ems_limit==1
    
    fprintf('plott levels with other emission limit')

        RES=OTHERPOL{plotts.regime_gov+1}; 
        RES_ems=OTHERPOL_Ems{plotts.regime_gov+1}; 

        allvars= RES("OPT_T_NoTaus");
        allvarsems= RES_ems("OPT_T_NoTaus");
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(allvars(find(varlist==varr),1:T)), time,allvarsems(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k';orrange} )   
            if lgdind==1
               lgd=legend('benchmark' , 'equal \% reduction' ,  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            if varr=="Tauf"
                ytickformat('%.0f')
            elseif varr =="dTaulAv"
                ytickformat('%.1f')
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            path=sprintf('figures/all_%s/Ems_Sens_%s_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr, indic.spillovers,plotts.nknk, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end

%% policy decomposition with less strict emission limit; comparison to benchmark decomposition
if plotts.Ems_limit_Decomp==1
    
    fprintf('plott decomposition ems limit')

        RES=OTHERPOL{plotts.regime_gov+1};
        RES_ems=OTHERPOL_Ems{plotts.regime_gov+1}; 

        if plotts.regime_gov==4
            RESnt=OTHERPOL{plotts.regime_gov+1+1}; 
            RES_emsnt=OTHERPOL_Ems{plotts.regime_gov+1+1};
        end
        allvars= RES("OPT_T_NoTaus");
        allvarsnt= RESnt("OPT_T_NoTaus");
        allvarsems= RES_ems("OPT_T_NoTaus");
        allvarsemsnt= RES_emsnt("OPT_T_NoTaus");

        Perdif=100*(allvars-allvarsnt)./allvarsnt;
        Perdifems=100*(allvarsems-allvarsemsnt)./allvarsemsnt;
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(Perdif(find(varlist==varr),1:T)), time,Perdifems(find(varlist==varr),1:T),time,zeros(1,T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-'; '--'; '--'}, {'color'}, {'k';orrange; grrey} )   
            if lgdind==1
               lgd=legend('benchmark' , 'equal \% reduction' ,  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            if varr=="sffsg"
                ytickformat('%.0f')
            elseif varr=="Tauf"
                ytickformat('%.1f')                
            else
                ytickformat('%.0f')
            end
            xticklabels(Year10)
            path=sprintf('figures/all_%s/EmsDecomp_CTO_Sens_%s_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr, indic.spillovers,plotts.nknk, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end
%% effect of taul in model without xgr and with xgr
% compare effect of only taul in benchmark and in exogenous growth model
if plotts.count_taul_xgr_lev==1
    
    fprintf('plott only taul xgr model level')

        allvars= RES_count_SAMETAUL_onlytaul("test");
        allvarst=RES_count_SAMETAUL_onlytaul('xgr');
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,allvars(find(varlist==varr),1:T),time,allvarst(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'}, {'color'}, {'k'; 'b'} )   
            if lgdind==1
               lgd=legend('benchmark' , 'exogenous growth',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/CountXgrTaul_target_%s_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end

%% effect tauf relative to laissez faire
% in levels
if plotts.count_tauflev==1
    
    fprintf('plott only tauf model level')

    %- loop over economy versions
    for k=keys(RES_count_onlytauf)
        kk=string(k);
            %- read in variable container of chosen regime
            if kk=="test"
                RES=OTHERPOL{plotts.regime_gov+1}; 
            elseif kk=="nsk"
                RES=OTHERPOL_nsk{plotts.regime_gov+1}; 
            elseif kk=="xgr"
                RES=OTHERPOL_xgr{plotts.regime_gov+1}; 
            elseif kk=="xgr_nsk"
                RES=OTHERPOL_xgr_nsk{plotts.regime_gov+1}; 
            end
        allvars= RES("LF");
        allvarst= RES_count_onlytauf(kk);
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,allvars(find(varlist==varr),1:T),time,allvarst(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'}, {'color'}, {'k'; 'b'} )   
            if lgdind==1
               lgd=legend('laissez-faire' , 'only environmental tax',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/CountTauf_mod%s_target_%s_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, kk, varr, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end
%% effect tauf relative to laissez faire
% in levels
if plotts.count_tauflev_Ben==1
    
    fprintf('plott only tauf model level and benchmark pol')

    %- loop over economy versions
    for k=keys(RES_count_onlytauf)
        kk=string(k);
            %- read in variable container of chosen regime
            if kk=="test"
                RES=OTHERPOL{plotts.regime_gov+1}; 
            elseif kk=="nsk"
                RES=OTHERPOL_nsk{plotts.regime_gov+1}; 
            elseif kk=="xgr"
                RES=OTHERPOL_xgr{plotts.regime_gov+1}; 
            elseif kk=="xgr_nsk"
                RES=OTHERPOL_xgr_nsk{plotts.regime_gov+1}; 
            end
        allvarsLF= RES("LF");
        allvars=RES("OPT_T_NoTaus");
        allvarst= RES_count_onlytauf(kk);
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,allvarsLF(find(varlist==varr),1:T),time,allvarst(find(varlist==varr),1:T),time,allvars(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-.';'--'; '-'}, {'color'}, {grrey; 'b'; 'k'} )   
            if lgdind==1
               lgd=legend('laissez-faire' , 'only environmental tax', 'both taxes optimal', 'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/CountTauf_Ben_mod%s_target_%s_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, kk, varr, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end
%% effect tauf relative to laissez faire
% in levels
if plotts.count_tauflev_Ben_noLF==1
    
    fprintf('plott only tauf model level and benchmark pol no lf')

    %- loop over economy versions
    for k=keys(RES_count_onlytauf)
        kk=string(k);
            %- read in variable container of chosen regime
            if kk=="test"
                RES=OTHERPOL{plotts.regime_gov+1}; 
            elseif kk=="nsk"
                RES=OTHERPOL_nsk{plotts.regime_gov+1}; 
            elseif kk=="xgr"
                RES=OTHERPOL_xgr{plotts.regime_gov+1}; 
            elseif kk=="xgr_nsk"
                RES=OTHERPOL_xgr_nsk{plotts.regime_gov+1}; 
            end
        allvars=RES("OPT_T_NoTaus");
        allvarst= RES_count_onlytauf(kk);
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,allvarst(find(varlist==varr),1:T),time,allvars(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'--'; '-'}, {'color'}, {'b'; 'k'} )   
            if lgdind==1
               lgd=legend( 'only environmental tax', 'both taxes optimal', 'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/CountTauf_Ben_noLF_mod%s_target_%s_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, kk, varr, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end
%% effect tauf relative to laissez faire
% in levels
if plotts.count_taullev==1
    
    fprintf('plott only taul model level')
%- loop over economy versions
    for k=keys(RES_count_onlytaul)
      kk=string(k);
            %- read in variable container of chosen regime
            if kk=="test"
                RES=OTHERPOL{plotts.regime_gov+1}; 
            elseif kk=="nsk"
                RES=OTHERPOL_nsk{plotts.regime_gov+1}; 
            elseif kk=="xgr"
                RES=OTHERPOL_xgr{plotts.regime_gov+1}; 
            elseif kk=="xgr_nsk"
                RES=OTHERPOL_xgr_nsk{plotts.regime_gov+1}; 
            end
        allvars= RES("LF");
        allvarst=RES_count_onlytaul(kk);

    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
     
            main=plot(time,allvars(find(varlist==varr),1:T),time,allvarst(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'}, {'color'}, {'k'; 'b'} )   
            if lgdind==1
               lgd=legend('laissez-faire' , 'only income tax',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/CountTaul_mod%s_target_%s_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, kk, varr, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end



%% Comparison model versions in one graph
% in levels
if plotts.count_modlev==1
    
    fprintf('plott counterfactual model level')

    %- read in variable container of chosen regime
    RES=OTHERPOL{plotts.regime_gov+1};
%     RESnsk=OTHERPOL_nsk{plotts.regime_gov+1};
%     RESxgr=OTHERPOL_xgr{plotts.regime_gov+1};  
    %- loop over economy versions
        allvars= RES("OPT_T_NoTaus");
        allvarsnsk=RES_count("nsk");
        allvarsxgr=RES_count("xgr");
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,allvars(find(varlist==varr),1:T),time,allvarsxgr(find(varlist==varr),1:T) ,time,allvarsnsk(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'; ':'}, {'color'}, {'k'; 'b'; orrange} )   
            if lgdind==1
               lgd=legend('benchmark' , 'exogenous growth', 'homogeneous skills',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/CountMod1_target_%s_regime%d_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr , plotts.regime_gov, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end

%% Comparison model versions: in deviation from efficient 1 graph
% in levels
if plotts.count_modlev_eff
    
    fprintf('plotting comp nsk xgr percent from efficient 1 counterfact')

    %- read in variable container of chosen regime
        RES=OTHERPOL{plotts.regime_gov+1};
        RESnsk=OTHERPOL_nsk{plotts.regime_gov+1};
        RESxgr=OTHERPOL_xgr{plotts.regime_gov+1};
        
        allvars= RES("OPT_T_NoTaus");
        allvarseff=RES("SP_T");
        allvarsnsk=RES_count("nsk");
        allvarsnskeff=RESnsk("SP_T");
        allvarsxgr=RES_count("xgr"); % counterfactual
        allvarsxgreff=RESxgr("SP_T");
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,100*(allvars(find(varlist==varr),1:T)-allvarseff(find(varlist==varr),1:T))./allvarseff(find(varlist==varr),1:T),time,100*(allvarsxgr(find(varlist==varr),1:T)-allvarsxgreff(find(varlist==varr),1:T))./allvarsxgreff(find(varlist==varr),1:T),...
                        time,(allvarsnsk(find(varlist==varr),1:T)-allvarsnskeff(find(varlist==varr),1:T))./allvarsnskeff(find(varlist==varr),1:T)*100, 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'; ':'}, {'color'}, {'k'; 'b'; orrange} )   
            if lgdind==1
                lgd=legend('benchmark' , 'exogenous growth', 'homogeneous skills',  'Interpreter', 'latex');
                
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/Per1_CountMod_%s_regime%d_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr , plotts.regime_gov, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end


%% Comparison model versions in one graph
% in levels
if plotts.compnsk_xgr1==1
    
    fprintf('plotting comp nsk xgr')

    %- read in variable container of chosen regime
    RES=OTHERPOL{plotts.regime_gov+1};
    RESnsk=OTHERPOL_nsk{plotts.regime_gov+1};
    RESxgr=OTHERPOL_xgr{plotts.regime_gov+1};  
    %- loop over economy versions
    for i = ["OPT_T_NoTaus" "OPT_NOT_NoTaus"]% only plotting polcies separately
        ii=string(i);
        allvars= RES(ii);
        allvarsnsk=RESnsk(ii);
        allvarsxgr=RESxgr(ii);
        
    fprintf('plotting %s',ii );
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,allvars(find(varlist==varr),1:T),time,allvarsxgr(find(varlist==varr),1:T) ,time,allvarsnsk(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'; ':'}, {'color'}, {'k'; 'b'; orrange} )   
            if lgdind==1
               lgd=legend('benchmark' , 'exogenous growth', 'homogeneous skills',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/CompMod1_%s_%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, ii,varr , plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end
    
%% Comparison model versions
% in levels
if plotts.compnsk_xgr==1
    
    fprintf('plotting comp nsk xgr')

    %- read in variable container of chosen regime
    RES=OTHERPOL{plotts.regime_gov+1};
    for kk="xgr_nsk" %["nsk" "xgr" ]
        if kk =="nsk" % no skill
            RESalt=OTHERPOL_nsk{plotts.regime_gov+1};
        elseif kk =="xgr"
            RESalt=OTHERPOL_xgr{plotts.regime_gov+1};
        elseif kk=="xgr_nsk"
            RESalt=OTHERPOL_xgr_nsk{plotts.regime_gov+1};
        end
  
    %- loop over economy versions
    for i = ["OPT_T_NoTaus" "OPT_NOT_NoTaus"]% only plotting polcies separately
        ii=string(i);
        allvars= RES(ii);
        allvarsalt=RESalt(ii);
        
    fprintf('plotting %s',ii );
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,allvars(find(varlist==varr),1:T),time,allvarsalt(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'}, {'color'}, {'k'; 'b'} )   
            if lgdind==1
                if kk == "nsk"
                    lgd=legend('benchmark' , 'homogeneous skills',  'Interpreter', 'latex');
                elseif kk=="xgr"
                    lgd=legend('benchmark' , 'exogenous growth',  'Interpreter', 'latex');
                elseif kk=="xgr_nsk"
                    lgd=legend('benchmark' , ['exogenous growth,' newline  'homogeneous skills'],  'Interpreter', 'latex');
                end
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/CompMod_%s_%s_%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, kk, ii,varr , plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
    end % loop over model versions
end
%% Comparison model versions: in deviation from efficient 1 graph
% in levels
if plotts.compEff_mod_dev1==1
    
    fprintf('plotting comp models efficient 1')

    %- read in variable container of chosen regime
        RES=OTHERPOL{plotts.regime_gov+1};
        RESnsk=OTHERPOL_nsk{plotts.regime_gov+1};
        RESxgr=OTHERPOL_xgr{plotts.regime_gov+1};

  
    %- loop over economy versions
    for ind = 1:2
        opt=["OPT_T_NoTaus" "OPT_NOT_NoTaus"];% only plotting polcies separately
        eff=["SP_T" "SP_NOT"];% only plotting polcies separately

        ii=string(opt(ind));
        ef=string(eff(ind));
        
        allvarseff=RES(ef);
        allvarsnskeff=RESnsk(ef);
        allvarsxgreff=RESxgr(ef);
        
    fprintf('plotting %s',ii );
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,allvarseff(find(varlist==varr),1:T),time,allvarsxgreff(find(varlist==varr),1:T),...
                        time,allvarsnskeff(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'; ':'}, {'color'}, {'k'; 'b'; orrange} )   
            if lgdind==1
                lgd=legend('benchmark' , 'exogenous growth', 'homogeneous skills',  'Interpreter', 'latex');
                
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/EFF_CompMod_%s_%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, ii,varr , plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end

%% Comparison model versions: in deviation from efficient
% in levels
if plotts.compnsk_xgr_dev==1
    
    fprintf('plotting comp nsk xgr percent from efficient')

    %- read in variable container of chosen regime
    RES=OTHERPOL{plotts.regime_gov+1};
    for kk= "xgr_nsk" %["nsk" "xgr"]
        if kk =="nsk" % no skill
            RESalt=OTHERPOL_nsk{plotts.regime_gov+1};
        elseif kk =="xgr"
            RESalt=OTHERPOL_xgr{plotts.regime_gov+1};
        elseif kk=="xgr_nsk"
            RESalt=OTHERPOL_xgr_nsk{plotts.regime_gov+1};
        end
  
    %- loop over economy versions
    for ind = 1:2
        opt=["OPT_T_NoTaus" "OPT_NOT_NoTaus"];% only plotting polcies separately
        eff=["SP_T" "SP_NOT"];% only plotting polcies separately

        ii=string(opt(ind));
        ef=string(eff(ind));
        
        allvars= RES(ii);
        allvarseff=RES(ef);
        allvarsalt=RESalt(ii);
        allvarsalteff=RESalt(ef);
        
    fprintf('plotting %s',ii );
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,100*(allvars(find(varlist==varr),1:T)-allvarseff(find(varlist==varr),1:T))./allvarseff(find(varlist==varr),1:T),time,100*(allvarsalt(find(varlist==varr),1:T)-allvarsalteff(find(varlist==varr),1:T))./allvarsalteff(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'}, {'color'}, {'k'; 'b'} )   
            if lgdind==1
                if kk == "nsk"
                    lgd=legend('benchmark' , 'homogeneous skills',  'Interpreter', 'latex');
                elseif kk=="xgr"
                    lgd=legend('benchmark' , 'exogenous growth',  'Interpreter', 'latex');
                elseif kk=="xgr_nsk"
                    lgd=legend('benchmark' , ['exogenous growth,' newline  'homogeneous skills'],  'Interpreter', 'latex');
                end
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/Per_CompMod_%s_%s_%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, kk, ii,varr , plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
    end % loop over model versions
end
%% Comparison model versions: in deviation from efficient 1 graph
% in levels
if plotts.compnsk_xgr_dev1==1
    
    fprintf('plotting comp nsk xgr percent from efficient 1')

    %- read in variable container of chosen regime
        RES=OTHERPOL{plotts.regime_gov+1};
        RESnsk=OTHERPOL_nsk{plotts.regime_gov+1};
        RESxgr=OTHERPOL_xgr{plotts.regime_gov+1};

  
    %- loop over economy versions
    for ind = 1:2
        opt=["OPT_T_NoTaus" "OPT_NOT_NoTaus"];% only plotting polcies separately
        eff=["SP_T" "SP_NOT"];% only plotting polcies separately

        ii=string(opt(ind));
        ef=string(eff(ind));
        
        allvars= RES(ii);
        allvarseff=RES(ef);
        allvarsnsk=RESnsk(ii);
        allvarsnskeff=RESnsk(ef);
        allvarsxgr=RESxgr(ii);
        allvarsxgreff=RESxgr(ef);
        
    fprintf('plotting %s',ii );
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,100*(allvars(find(varlist==varr),1:T)-allvarseff(find(varlist==varr),1:T))./allvarseff(find(varlist==varr),1:T),time,100*(allvarsxgr(find(varlist==varr),1:T)-allvarsxgreff(find(varlist==varr),1:T))./allvarsxgreff(find(varlist==varr),1:T),...
                        time,(allvarsnsk(find(varlist==varr),1:T)-allvarsnskeff(find(varlist==varr),1:T))./allvarsnskeff(find(varlist==varr),1:T)*100, 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'; ':'}, {'color'}, {'k'; 'b'; orrange} )   
            if lgdind==1
                lgd=legend('benchmark' , 'exogenous growth', 'homogeneous skills',  'Interpreter', 'latex');
                
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/Per1_CompMod_%s_%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, ii,varr , plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end

%% All figures single
if plotts.single_pol==1
    
    fprintf('plotting single graphs')

    %- read in variable container of chosen regime
    RES=OTHERPOLL{plotts.regime_gov+1};
%     RES=ccc;
  

    %- loop over economy versions
    for i = ["OPT_T_NoTaus" "OPT_NOT_NoTaus"]% only plotting polcies separately
        ii=string(i);
        allvars= RES(ii);
    fprintf('plotting %s',ii );
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
            if ll=="HH" && varr=="Emnet"
                main=plot(time,allvars(find(varlist==varr),1:T),time(percon+1:end),Ems(1:T), 'LineWidth', 1.1);  
                set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
                lgd=legend('net emissions' , 'net emission limit',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');

            else
                main=plot(time,allvars(find(varlist==varr),1:T), 'LineWidth', 1.1);   
                set(main, {'LineStyle'},{'-'}, {'color'}, {'k'} )   

            end
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if indic.count_techgap==0
                path=sprintf('figures/all_%s/Single_%s_%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_extern%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.png',date,  ii,varr , plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr,indic.extern, indic.PV, plotts.sizeequ, plotts.GOV,   etaa);
            else
                path=sprintf('figures/all_%s/Single_%s_%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_extern%d_countec_PV%d_sizeequ%d_GOV%d_etaa%.2f.png',date,  ii, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.extern, indic.PV,plotts.sizeequ, plotts.GOV, etaa);
            end
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end
%% figures single overlayed
if plotts.singov==1
    fprintf('plotting single overlayed graphs')
    RES=OTHERPOLL{plotts.regime_gov+1};
    varl=varlist;
    for lgdind=0:1
    for i =["OPT_T_NoTaus" "OPT_NOT_NoTaus"]
        ii=string(i);
        allvars =RES(ii);
    fprintf('plotting %s',ii );
    for l =keys(lissComp) % loop over variable groups
        ll=string(l);
        plotvars=lissComp(ll); % here plotvars is a group of variable names which are to be plotted in the same graph

        gcf=figure('Visible','off');

        if length(plotvars)==2
            main=plot(time,allvars(find(varl==plotvars(1)),1:T), time,allvars(find(varl==plotvars(2)),1:T), 'LineWidth', 1.1);    % plot vectors!        
            set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; 'k'} ) 
    %    elseif ll=="LabourInp"
    %           main=plot(time,allvars(find(varlist==plotvars(1)),:), time,allvars(find(varlist==plotvars(2)),:),time,allvars(find(varlist==plotvars(3)),:), 'LineWidth', 1.1);    % plot vectors!        
    %           set(main, {'LineStyle'},{'-'; '--'; ':'}, {'color'}, {'k'; 'k'; 'k'} )   
       elseif ll=="Growth"
              main=plot(time(1:T),allvars(find(varl==plotvars(1)),1:T), time(1:T),allvars(find(varl==plotvars(2)),1:T),...
              time(1:T),allvars(find(varl==plotvars(3)),1:T),time(1:T),allvars(find(varl==plotvars(4)),1:T),'LineWidth', 1.1);    % plot vectors!        
              set(main, {'LineStyle'},{'-'; '--'; ':'; '--'}, {'color'}, {'k'; 'k'; 'k'; grrey} )   

        elseif ll=="Science"
              main=plot(time,allvars(find(varl==plotvars(1)),1:T), time,allvars(find(varl==plotvars(2)),1:T),...
                  time,allvars(find(varl==plotvars(3)),1:T), time,allvars(find(varl==plotvars(4)),1:T),'LineWidth', 1.1);    % plot vectors!        
               set(main, {'LineStyle'},{'-'; '--'; ':'; '--'}, {'color'}, {'k'; 'k'; 'k'; grrey} )   
        end
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
         if lgdind==1
            pp=legg(ll);
            if length(pp)==2
                lgd=legend(sprintf('%s',pp(1)) ,sprintf('%s',pp(2)),  'Interpreter', 'latex');
            elseif length(pp)==3
                lgd=legend(sprintf('%s',pp(1)) ,sprintf('%s',pp(2)),sprintf('%s',pp(3)),  'Interpreter', 'latex');

            elseif length(pp)==4
                lgd=legend(sprintf('%s',pp(1)) ,sprintf('%s',pp(2)),sprintf('%s',pp(3)), sprintf('%s',pp(4)),  'Interpreter', 'latex');
            end
            set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
         end
        if indic.count_techgap==0
            path=sprintf('figures/all_%s/SingleJointTOT_regime%d_%s_%s_spillover%d_noskill%d_sep%d_xgrowth%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, plotts.regime_gov, ii,ll, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr,indic.extern, indic.PV,  etaa, lgdind);
        else
            path=sprintf('figures/all_%s/SingleJointTOT_regime%d_%s_%s_spillover%d_noskill%d_sep%d_xgrowth%d_countec_extern%d_PV%d_etaa%.2f_lgd%d.png',date,  plotts.regime_gov, ii,ll, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.extern, indic.PV, etaa, lgdind);
        end
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
    end
    end
    end
end

 %% comparison POLICY scenarios
if plotts.notaul==1
    fprintf('plotting comparison across policies') 
    bb=1:length(OTHERPOLL);
    bb=bb(bb~=plotts.regime_gov+1); % drop benchmark policy
    for nt =  bb
        pp = OTHERPOLL{nt};
        count = nt-1; % subtract 1 to get indicator of notaul
        
    for i =["OPT_T_NoTaus" "OPT_NOT_NoTaus"]

     ii=string(i);
     varl=varlist; 
     %- benchmark policy
     helpp=OTHERPOLL{plotts.regime_gov+1};
     allvars= helpp(ii);
     %- comparison
     allvarsnt=pp(ii); 
     
    %% 
    fprintf('plotting %s',ii );
    for lgdind=0:1
    for l =keys(lisst)
        ll=string(l);
      
        plotvars=lisst(ll);
      
        for v=1:length(plotvars)
           gcf=figure('Visible','off');

               varr=string(plotvars(v));
               main=plot(time,allvars(find(varl==varr),1:T), time,allvarsnt(find(varlist==varr),1:T), 'LineWidth', 1.1);

               set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; 'b'} )   
               xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
               ax=gca;
               ax.FontSize=13;
               ytickformat('%.2f')
               xticklabels(Year10)

            if lgdind==1
                if count==0 % integrated policy
                    lgd=legend('benchmark', 'integrated policy, with income tax', 'Interpreter', 'latex');
                elseif count == 1 % integrated without income tax
                     lgd=legend('benchmark', 'integrated policy, no income tax', 'Interpreter', 'latex');
                elseif count == 2 % notaul =2 gov=tauf*F*pf, without income tax
                     lgd=legend('with income tax', 'without income tax', 'Interpreter', 'latex');
                elseif count == 3 % notaul=3
                     lgd=legend('benchmark', 'no redistribution, with income tax', 'Interpreter', 'latex');
                elseif count == 4 % Tls and taul
                     lgd=legend('benchmark', 'lump-sum transfers, with income tax', 'Interpreter', 'latex');
                elseif count == 5 % Tls no taul
                     lgd=legend('benchmark', 'lump-sum transfers, no income tax', 'Interpreter', 'latex');
                end
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            path=sprintf('figures/all_%s/comp_benchregime%d_notaul%d_%s_%s_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png', date, plotts.regime_gov, count, ii, varr, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV, etaa, lgdind);
            exportgraphics(gcf,path,'Resolution', 400)
            close gcf
           end % variables in group
    end % variable group
    end % legend
    end
    end
end

%% comparison with and without target    
%- string to loop over 
if plotts.comptarg==1
    fprintf('plotting comparison target graphs')
    ssr= string({'SP_T', 'SP_NOT' });%,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'

    for i =[1] % to filter benchmark: policy with target
        ii=ssr(i);
        %- read in data
        t=string(ssr(i));
        nt=string(ssr(i+1));
        if indic.count_techgap==0
                RES=OTHERPOLL{plotts.regime_gov+1};

                 allvars= RES(t);
                 allvarsnot=RES(nt); 
        else
                 allvars =RES_count(t);
                 allvarsnot=RES_count(nt); 
        end

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
            main=plot(time,allvars(find(varlist==varr),1:T), time,allvarsnot(find(varlist==varr),1:T), 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; 'k'} )   
           xticks(txx)
           xlim([1, time(end-1)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('wih emission limit', 'no emission limit', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 18,'Orientation', 'vertical');
           end
           if indic.count_techgap==0
                path=sprintf('figures/all_%s/%s_TargetComp%s_regime%d_knspil%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png', date, varr, ii, plotts.regime_gov,indic.noknow_spill, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV, etaa, lgdind);
           else
                path=sprintf('figures/all_%s/%s_TargetComp%s_countec_knspil0_spillover0_noskill0_sep0_xgrowth0_PV1_etaa%.2f_lgd%d.png', date, varr, ii, etaa, lgdind);
           end
            exportgraphics(gcf,path,'Resolution', 400)
       close gcf
        end
        end
    end
    end      
end

%% comparison social planner and optimal policy benchmark and comparison policy
% graph incorporates with and without laissez faire allocation 
if plotts.compeff==1
    fprintf('plotting comparison efficient-optimal and non benchmark graphs')   

    %- read in container of results
    RES=OTHERPOLL{plotts.regime_gov+1};
    bb=1:length(OTHERPOLL);
    bb=bb(bb~=plotts.regime_gov+1); % drop benchmark policy
    
    for withlff=0
        lff=RES('LF');
        eff= string({'SP_T', 'SP_NOT'});
        opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});

    for i =1:length(eff)

        ie=eff(i);
        io=opt(i);
        varl=varlist;

         %- benchmark policy
         allvars= RES(io);
         allvarseff=RES(ie); 

     %- comparison         
     for nt = 3% bb % loop over policy scenarios but benchmark
        RES_help=OTHERPOLL{nt};
        count=nt-1;
        allvarsnotaul =RES_help(io);

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
           if withlff==1
               main=plot(time, lff(find(varl==varr),1:T),  time,allvars(find(varl==varr),1:T), time,allvarsnotaul(find(varl==varr),1:T), time,allvarseff(find(varl==varr),1:T));            
               set(main,{'LineWidth'}, {1; 1.2; 1.2; 1},  {'LineStyle'},{'--';'-'; '--'; ':'}, {'color'}, {grrey; 'k'; orrange; 'b'} )   
           else
               main=plot( time,allvars(find(varl==varr),1:T), time,allvarsnotaul(find(varl==varr),1:T), time,allvarseff(find(varl==varr),1:T));            
               set(main,{'LineWidth'}, { 1.2; 1.2; 1},  {'LineStyle'},{'-'; '--'; ':'}, {'color'}, {'k'; 'b'; orrange} )   
           end
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)

           if lgdind==1
               if withlff==1

                   if count ==0
                        lgd=legend('laissez-faire',  'benchmark policy',  'integrated policy, with income tax', 'efficient','Interpreter', 'latex');
                   elseif count ==1
                        lgd=legend('laissez-faire','benchmark policy', 'integrated policy, no income tax',   'efficient', 'Interpreter', 'latex');
                   elseif count ==2
                        lgd=legend('laissez-faire', 'with income tax', 'without income tax', 'efficient',  'Interpreter', 'latex'); 
                   elseif count==3
                        lgd=legend('laissez-faire', 'benchmark policy', 'no redistribution, with income tax', 'efficient',  'Interpreter', 'latex');
                   elseif count==4
                        lgd=legend('laissez-faire', 'benchmark policy', 'lump-sum transfers, with income tax', 'efficient',  'Interpreter', 'latex');
                   elseif count==5
                        lgd=legend('laissez-faire',  'benchmark policy', 'lump-sum transfers, no income tax', 'efficient',  'Interpreter', 'latex');                        
                   end
                   
                  set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
               else
                   if count ==0
                        lgd=legend(  'benchmark policy', 'integrated policy, with income tax',  'efficient', 'Interpreter', 'latex');
                   elseif count ==1
                        lgd=legend( 'benchmark policy', 'integrated policy, no income tax', 'efficient',   'Interpreter', 'latex');
                   elseif count ==2
                        lgd=legend(  'with income tax', 'without income tax',  'efficient', 'Interpreter', 'latex'); 
                   elseif count==3
                        lgd=legend(  'benchmark policy', 'no redistribution, with income tax', 'efficient',  'Interpreter', 'latex');
                   elseif count==4
                        lgd=legend( 'benchmark policy', 'lump-sum transfers, with income tax',  'efficient', 'Interpreter', 'latex');
                   elseif count==5
                        lgd=legend( 'benchmark policy', 'lump-sum transfers, no income tax', 'efficient',  'Interpreter', 'latex');                        
                   end
                   
                   
                  set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
               end
           end
        path=sprintf('figures/all_%s/%s_CompEff%s_benchregime%d_pol%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d_lff%d.png',date, varr, io, plotts.regime_gov, count, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV, etaa, lgdind, withlff);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end
    end
    end
end


%% only social planner
if plotts.compeff1==1
    fprintf('plotting efficient')   

    %- read in container of results: any fine for social planner
    
    eff= ["SP_T" "SP_NOT"];   
    for i =[1,2]

        ie=eff(i);
        if indic.count_techgap==0
            
            RES=OTHERPOLL{plotts.regime_gov+1};
            allvarseff=RES(ie); 
        else
            allvarseff =RES_count(ie);
        end

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
           
           main=plot( time,allvarseff(find(varlist==varr),1:T));            
           set(main,{'LineWidth'}, {1.2},  {'LineStyle'},{'-'}, {'color'}, {'k'} )   
           xticks(txx)
           xlim([1, time(end-1)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)

           if lgdind==1
              lgd=legend('efficient', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
           end
           if indic.count_techgap==0
                path=sprintf('figures/all_%s/%s_CompEff%s_onlyeff_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_countec%d_PV%d_etaa%.2f_lgd%d.png', date, varr, ie, indic.spillovers,indic.noknow_spill, plotts.nsk, indic.sep, plotts.xgr, indic.count_techgap, indic.PV, etaa, lgdind);
           else
                path=sprintf('figures/all_%s/%s_CompEff%s_count_onlyeff_spillover0_knspil0_noskill0_sep0_xgrowth0_PV1_etaa%.2f_lgd%d.png', date, varr, ie, etaa, lgdind);

           end
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end
end

%% comparison social planner and non-benchmark policy
if plotts.compeff2==1
    %- only efficient and no income tax
    fprintf('plotting comparison efficient-non benchmark optimal graphs')   
    
    for withlff=0
    
        eff= string({'SP_T', 'SP_NOT'});
        opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});
   
    % for withtaul=0:1
    for i =[1,2]

        ie=eff(i);
        io=opt(i);
    for nt =  1:length(OTHERPOLL) % loop over policy regimes
        count=nt-1;
        RES=OTHERPOLL{nt};
        lff=RES('LF');
        allvarsnotaul =RES(io);
        allvarseff=RES(ie); 

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
      if withlff==1
            main=plot(time, lff(find(varlist==varr),1:T), time,allvarseff(find(varlist==varr),1:T), time,allvarsnotaul(find(varlist==varr),1:T));            
           set(main, {'LineWidth'}, {1; 1.2; 1.2}, {'LineStyle'},{'--';'-'; '--'}, {'color'}, {grrey; 'k'; orrange} )   
      else
            main=plot(time,allvarseff(find(varlist==varr),1:T), time,allvarsnotaul(find(varlist==varr),1:T));            
           set(main, {'LineWidth'}, {1.2; 1.2}, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
      end
      xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
          
           if lgdind==1
               if withlff==1
                   if count ==0
                        lgd=legend('laissez-faire', 'efficient', 'integrated policy, with income tax', 'Interpreter', 'latex');
                   elseif count ==1
                        lgd=legend('laissez-faire', 'efficient',  'integrated policy, no income tax',  'Interpreter', 'latex');
                   elseif count ==2
                        lgd=legend('laissez-faire', 'efficient',  'no redistribution, no income tax', 'Interpreter', 'latex'); 
                   elseif count==3
                        lgd=legend('laissez-faire', 'efficient',  'no redistribution, with income tax', 'Interpreter', 'latex');
                   elseif count==4
                        lgd=legend('laissez-faire', 'efficient', 'lump-sum transfers, with income tax', 'Interpreter', 'latex');
                   elseif count==5
                        lgd=legend('laissez-faire', 'efficient', 'lump-sum transfers, no income tax', 'Interpreter', 'latex');                        
                   end
               else
                  if count ==0
                        lgd=legend( 'efficient', 'integrated policy, with income tax', 'Interpreter', 'latex');
                   elseif count ==1
                        lgd=legend( 'efficient',  'integrated policy, no income tax',  'Interpreter', 'latex');
                   elseif count ==2
                        lgd=legend( 'efficient', 'no redistribution, no income tax', 'Interpreter', 'latex'); 
                   elseif count==3
                        lgd=legend( 'efficient',  'no redistribution, with income tax', 'Interpreter', 'latex');
                   elseif count==4
                        lgd=legend( 'efficient',  'lump-sum transfers, with income tax', 'Interpreter', 'latex');
                   elseif count==5
                        lgd=legend( 'efficient',  'lump-sum transfers, no income tax', 'Interpreter', 'latex');                        
                   end
               end
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
           end
        path=sprintf('figures/all_%s/%s_CompEff%s_noopt_pol%d_spillover%d_noskill%d_sep%d_xgrowth%d_countec%d_PV%d_etaa%.2f_lgd%d_lff%d.png', date, varr, io, count, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr, indic.count_techgap, indic.PV, etaa, lgdind, withlff);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end
    end
    end
end
%% LF versus optimal policy in levels
if plotts.comp_LFOPT==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting levels with LF and no taul no eff') 

     if plotts.regime_gov==3
        RES=OTHERPOLL{plotts.regime_gov+1};
        RESnt =OTHERPOLL{3}; % version without taul
     end
        allvars= RES('OPT_T_NoTaus');
        revall =RES('LF');
        allvarsnt =RESnt('OPT_T_NoTaus');
        
    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot(time,(revall(find(varlist==varr),1:T)), time,(allvars(find(varlist==varr),1:T)),time,(allvarsnt(find(varlist==varr),1:T)),'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-.';'-'; '--'}, {'color'}, {grrey; 'k'; 'b'} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('laissez-faire', 'with income tax', 'without income tax', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_LF_OPT_regime%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, indic.noknow_spill, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
   end
end

%% level differences: NEW CALIBRATION with and without taul
if plotts.comp_OPT_NCAlib==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting levels opt and no taul no eff') 

        allvars= RES_NCalib('OPT_T_WithTaul');
        allvarsnt =RES_NCalib('OPT_T_NOTaul');
        
    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot( time,(allvars(find(varlist==varr),1:T)),time,(allvarsnt(find(varlist==varr),1:T)),'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; 'b'} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('with income tax', 'without income tax', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_OPT_COMPtaul_regime%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nknk, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
   end
end

%% optimal with and without taul in levels 
if plotts.comp_OPT==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting levels opt and no taul no eff') 

     if plotts.regime_gov==3
        RES=OTHERPOLL{plotts.regime_gov+1};
        RESnt =OTHERPOLL{3}; % version without taul
     end
        allvars= RES('OPT_T_NoTaus');
        allvarsnt =RESnt('OPT_T_NoTaus');
        
    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot( time,(allvars(find(varlist==varr),1:T)),time,(allvarsnt(find(varlist==varr),1:T)),'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; 'b'} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('with income tax', 'without income tax', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_OPT_COMPtaul_regime%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, indic.noknow_spill, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
   end
end

%% counterfactual optimal policy in NKS model in benchmark
if plotts.comp_Bench_CountNK==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting opt pol without kn spills in benchmark model') 

     if plotts.regime_gov==3
        RES=OTHERPOLL{plotts.regime_gov+1};
     end
        allvars= RES('OPT_T_NoTaus');
        allvarsCOUNT =RES_count_NKinBen('all');
        
    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot( time,(allvars(find(varlist==varr),1:T)),time,(allvarsCOUNT(find(varlist==varr),1:T)),'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'--'; '-'}, {'color'}, {grrey; 'k'} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('optimal policy', 'counterfactual policy', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_KNCOUNT_FullMod_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
   end
end
%% optimal with and without taul in levels 
if plotts.comp_OPT_NK==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting levels opt benchmark and no knowledge spillovers') 

        RES=OTHERPOLL{plotts.regime_gov+1};
        allvars= RES('OPT_T_NoTaus');
        allvarsNK =RES_noknspil('OPT_T_NoTaus');
        
    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot( time,(allvarsNK(find(varlist==varr),1:T)),time,(allvars(find(varlist==varr),1:T)),'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; grrey} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('no knowledge spillovers', 'benchmark', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_KN_FullMod_sizeequ%d_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_GOV%d_etaa%.2f_lgd%d.png',date, varr, plotts.sizeequ, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,plotts.GOV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
   end
end
%% comparison social planner and optimal policy benchmark
if plotts.compeff3==1
    %- only efficient and benchmark
    fprintf('plotting comparison efficient-optimal graphs')   
    RES=OTHERPOLL{plotts.regime_gov+1};

    for withlff=0
         lff=RES('LF');
    
        eff= string({'SP_T', 'SP_NOT'});
        opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});

    for i =[1,2]

        ie=eff(i);
        io=opt(i);
        
        allvars= RES(io);
        varl=varlist;
        allvarseff=RES(ie); 

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
      if withlff==1
            main=plot(time, lff(find(varl==varr),1:T), time,allvarseff(find(varl==varr),1:T), time,allvars(find(varl==varr),1:T));            
           set(main, {'LineWidth'}, {1; 1.2; 1.2}, {'LineStyle'},{'--';'-'; '--'}, {'color'}, {grrey; 'k'; orrange} )   
      else
            main=plot(time,allvarseff(find(varl==varr),1:T), time,allvars(find(varl==varr),1:T));            
           set(main, {'LineWidth'}, {1.2; 1.2}, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
      end
      xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)

           if lgdind==1
               if withlff==1
                    lgd=legend('laissez-faire', 'efficient', 'optimal policy', 'Interpreter', 'latex');
               else
%                    if varr =="tauf"
%                       lgd=legend( 'social cost of emissions', 'no income tax', 'Interpreter', 'latex');
%                    else
                     lgd=legend( 'efficient', ' optimal policy', 'Interpreter', 'latex');
%                    end
               end
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
           end
        path=sprintf('figures/all_%s/%s_CompEff%s_regime%d_opteff_spillover%d_noskill%d_sep%d_xgrowth%d_countec%d_PV%d_etaa%.2f_lgd%d_lff%d.png', date, varr, io, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr, indic.count_techgap, indic.PV, etaa, lgdind, withlff);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end
    end
end
%% comparison social planner and optimal policy benchmark
if plotts.compeff3_NC==1
    %- only efficient and benchmark
    fprintf('plotting comparison efficient-optimal graphs')   
  

    for withlff=1
    lff= RES_NCalib("LF");
    
    for s=["T", "NoT"]
        ss=string(s);
        allvarseff=RES_NCalib(sprintf("SP_%s",ss) );
        allvars=RES_NCalib(sprintf("OPT_%s_WithTaul",ss));
%    eff= string({'SP_T', 'SP_NOT'});
%         opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});
% 
%     for i =[1,2]
% 
%         ie=eff(i);
%         io=opt(i);
%         
%         allvars= RES(io);
         varl=varlist;
%         allvarseff=RES(ie); 

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
      if withlff==1
            main=plot(time, lff(find(varl==varr),1:T), time,allvarseff(find(varl==varr),1:T), time,allvars(find(varl==varr),1:T));            
           set(main, {'LineWidth'}, {1; 1.2; 1.2}, {'LineStyle'},{'--';'-'; '--'}, {'color'}, {grrey; 'k'; orrange} )   
      else
            main=plot(time,allvarseff(find(varl==varr),1:T), time,allvars(find(varl==varr),1:T));            
           set(main, {'LineWidth'}, {1.2; 1.2}, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
      end
      xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)

           if lgdind==1
               if withlff==1
                    lgd=legend('laissez-faire', 'efficient', 'optimal policy', 'Interpreter', 'latex');
               else
%                    if varr =="tauf"
%                       lgd=legend( 'social cost of emissions', 'no income tax', 'Interpreter', 'latex');
%                    else
                     lgd=legend( 'efficient', ' optimal policy', 'Interpreter', 'latex');
%                    end
               end
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
           end
        path=sprintf('figures/all_%s/NC_%s_CompEff%s_regime%d_opteff_spillover%d_noskill%d_sep%d_xgrowth%d_countec%d_PV%d_etaa%.2f_lgd%d_lff%d.png', date, varr, ss, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr, indic.count_techgap, indic.PV, etaa, lgdind, withlff);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end
    end
end
%% comparison to BAU
if plotts.bau==1
    fprintf('plotting bau graphs') 
for nt =  1:length(OTHERPOLL) % loop over policy regimes
        count=nt-1;
        RES=OTHERPOLL{nt};
        bau=RES('BAU');

    for i ={'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'}
        ii=string(i);
        allvars= RES(ii);

    %% 
    fprintf('plotting %s',ii );
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');


            varr=string(plotvars(v));
            %subplot(floor(length(plotvars)/nn)+1,nn,v)
            main=plot(time,allvars(find(varlist==varr),1:T), time,bau(find(varlist==varr),1:T), 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
               if contains(ii, 'SP')
                  lgd=legend('Social planner', 'status quo', 'Interpreter', 'latex');
               else
                  lgd=legend('Ramsey planner', 'status quo', 'Interpreter', 'latex');
               end
            set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_BAUComp%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, ii, count, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end      
end
end

%% comparison to LF
if plotts.lf==1
    fprintf('plotting LF graphs') 
for nt =  1:length(OTHERPOLL) % loop over policy regimes
        count=nt-1;
        RES=OTHERPOLL{nt};
        bau=RES('LF');

    for i ={'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'}
        ii=string(i);
        allvars= RES(ii);

    %% 
    fprintf('plotting %s',ii );
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');


            varr=string(plotvars(v));
            %subplot(floor(length(plotvars)/nn)+1,nn,v)
            main=plot(time,allvars(find(varlist==varr),1:T), time,bau(find(varlist==varr),1:T), 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
               if contains(ii, 'SP')
                  lgd=legend('Social planner', 'laissez-faire', 'Interpreter', 'latex');
               else
                  lgd=legend('Ramsey planner', 'laissez-faire', 'Interpreter', 'latex');
               end
            set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_LFComp%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, ii, count, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end      
end
end
%%
if plotts.per_BAUt0==1
    % plot graphs in percent relative to LF, efficient/optimal world
    % without tagret dynamic and with t, percentage change over time relative to BAU scenario
    % in t=0;
     fprintf('plotting percentage relative to bau eff vs opt') 

        RES=OTHERPOLL{plotts.regime_gov+1};
        bau=RES('BAU');
               
        allvars= RES('OPT_T_NoTaus');
        allvarseff= RES('SP_T');

    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            rev = bau(find(varlist==varr), 1);
            
            main=plot(time,(allvarseff(find(varlist==varr),1:T)-rev)./rev,time,(allvars(find(varlist==varr),1:T)-rev)./rev, 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('efficient', 'optimal policy', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageBAUComp_Target_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    
end

%%
if plotts.per_LFt0==1
    % plot graphs in percent relative to LF, efficient/optimal world
    % without tagret dynamic and with t, percentage change over time relative to BAU scenario
    % in t=0;
     fprintf('plotting percentage relative to bau eff vs opt') 

        RES=OTHERPOLL{plotts.regime_gov+1};
        bau=RES('BAU');
               
        allvars= RES('OPT_T_NoTaus');
        allvarseff= RES('SP_T');

    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            rev = bau(find(varlist==varr), 1);
            
            main=plot(time,(allvarseff(find(varlist==varr),1:T)-rev)./rev,time,(allvars(find(varlist==varr),1:T)-rev)./rev, 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('efficient', 'optimal policy', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageLFComp_Target_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    
end
%%
 if plotts.per_effopt0==1
    % plot graphs in percent relative to LF, efficient/optimal world
    % without tagret dynamic and with t, percentage change over time relative to BAU scenario
    % in t=0;
     fprintf('plotting percentage relative to eff vs opt first period') 

        RES=OTHERPOLL{plotts.regime_gov+1};
               
        allvars= RES('OPT_T_NoTaus');
        allvarseff= RES('SP_T');
        revall =RES('OPT_NOT_NoTaus');
        revalleff =RES('SP_NOT');

    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            revopt = revall(find(varlist==varr), 1);
             reveff = revalleff(find(varlist==varr), 1);
            
            main=plot(time,(allvarseff(find(varlist==varr),1:T)-reveff)./reveff,time,(allvars(find(varlist==varr),1:T)-revopt)./revopt, 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('efficient', 'optimal policy', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageEffOptFirstPeriod_Target_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    
 end   

%%
 if plotts.per_optd==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting percentage relative to opt dynamic') 

        RES=OTHERPOLL{plotts.regime_gov+1};
               
        allvars= RES('OPT_T_NoTaus');
        revall =RES('OPT_NOT_NoTaus');
        
    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot(time,(allvars(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T), 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'}, {'color'}, {'k'} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend( 'optimal policy', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageOptDyn_Target_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    
 end   

 %%
 % expressed in reduction relative to bau per period/ dynamic
if plotts.per_baud==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting percentage relative to bau dynamic') 

        RES=OTHERPOLL{plotts.regime_gov+1};
               
        allvars= RES('OPT_T_NoTaus');
        allvarseff= RES('SP_T');
        revall =RES('BAU');

    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot(time,(allvarseff(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T),time,(allvars(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T), 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('efficient', 'optimal policy', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageBAUDyn_Target_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    
end   
 
 %%
 % expressed in reduction relative to bau per period/ dynamic
if plotts.per_LFd==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting percentage relative to bau dynamic') 

        RES=OTHERPOLL{plotts.regime_gov+1};
               
        allvars= RES('OPT_T_NoTaus');
        allvarseff= RES('SP_T');
        revall =RES('LF');

    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot(time,(allvarseff(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T),time,(allvars(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T), 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('efficient', 'optimal policy', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageLFDyn_Target_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    
end   
 
 %% expressed in reduction relative to bau per period/ dynamic
if plotts.per_LFd_NC==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting percentage relative to bau dynamic new calib') 

    revall= RES_NCalib("LF");
    for s=["T", "NoT"]
        ss=string(s);
        allvarseff=RES_NCalib(sprintf("SP_%s",ss) );
       allvars=RES_NCalib(sprintf("OPT_%s_WithTaul",ss));


%         RES=OTHERPOLL{plotts.regime_gov+1};
%                
%         allvars= RES('OPT_T_NoTaus');
%         allvarseff= RES('SP_T');
%         revall =RES('LF');

    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot(time,(allvarseff(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T),time,(allvars(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T), 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('efficient', 'optimal policy', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/NC_%s_PercentageLFDyn_%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr,ss, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end
end   
 
%%
 % expressed in reduction relative to bau per period/ dynamic
if plotts.per_LFd_nt==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting percentage relative to LF dynamic plus no taul') 

     if plotts.regime_gov==3
        RES=OTHERPOLL{plotts.regime_gov+1};
        RESnt =OTHERPOLL{3}; 
     end
        allvars= RES('OPT_T_NoTaus');
        allvarseff= RES('SP_T');
        revall =RES('LF');
        allvarsnt =RESnt('OPT_T_NoTaus');
        
    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot(time,(allvars(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T)*100,time,(allvarsnt(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T)*100, time,100*(allvarseff(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T),'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'; ':'}, {'color'}, {'k'; 'b'; orrange} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('with income tax', 'without income tax', 'efficient', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageLFDynNT_Target_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
end
%% relative change lF dynamic no efficient
if plotts.per_LFd_ne_nt==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting percentage relative to LF dynamic plus no taul no eff') 

     if plotts.regime_gov==3
        RES=OTHERPOLL{plotts.regime_gov+1};
        RESnt =OTHERPOLL{3}; 
     end
        allvars= RES('OPT_T_NoTaus');
        allvarseff= RES('SP_T');
        revall =RES('LF');
        allvarsnt =RESnt('OPT_T_NoTaus');
        
    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot(time,(allvars(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T)*100,time,(allvarsnt(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T)*100,'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; 'b'} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('with income tax', 'without income tax', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageLFDynNT_noeff_Target_regime%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, indic.noknow_spill, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
end
end     