function [symms, list, opt_all]= OPT_solve(list, symms, params, x0LF, init201519, indexx, indic, T, Ems, MOM, percon)

% pars
read_in_params;
Ftarget =  (Ems'+deltaa)/omegaa;

if indic.notaul==9 % taking labor tax from optimal policy without target
    helper  = load(sprintf('OPT_notarget_2112_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul4_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.sep, indic.extern , indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
    taulFixed           =helper.opt_all(:,list.allvars=='taul');
else
    taulFixed=zeros(size(Ftarget));
end
% symbilic variables and lists
syms H C F G Af Ag An sff sn sg Lf Lg real

if indic.xgrowth==1
    symms.opt= [Lf Lg C F G H];
else
    symms.opt = [Lf Lg C F G H sn sff sg];
end

list.opt  = string(symms.opt); 
nn= length(list.opt); % number of variables

%- load old lists to feed in initial solution
% lisst=load("../../optimalPol_010922_revision/listallvars_old.mat");
% lissthelp=lisst.list;


%% Initial Guess %%% 
%%%%%%%%%%%%%%%%%%%%%

% laissez faie
if indic.target==0
    helper=load(sprintf('../output/OPT_notarget_2407_notaul%d_subsres%d_phii%d_sigma%d_Bop%d_util%d_PV%d.mat',...
            2,indic.subsres, indic.know_spill, indic.sigmaWorker, indic.Bop, indic.util));
else
    helper=load(sprintf('../output/OPT_target_2407_emnet0_notaul%d_subsres%d_phii%d_sigma%d_Bop%d_util%d_PV%d.mat',...
            2,indic.subsres, indic.know_spill, indic.sigmaWorker, indic.Bop, indic.util));
end

opt_all=helper.opt_all;
x0 = zeros(nn*T,1);
% listsave=list;
% list.allvars=list.allvars(list.allvars~='taurese'& list.allvars~='tauresg'& list.allvars~='taus');

%%
if indic.target==1
        kappaa= Ftarget'./opt_all(1:T,list.allvars=='F')';    
        kappaa = kappaa*(1-1e-10);
    
        x0(T*(find(list.opt=='Lf')-1)+1:T*(find(list.opt=='Lf'))) =opt_all(:,list.allvars=='Lf'); % hhf; first period in LF is baseline
        x0(T*(find(list.opt=='Lg')-1)+1:T*(find(list.opt=='Lg'))) =opt_all(:,list.allvars=='Lg'); % hhg
        x0(T*(find(list.opt=='H')-1)+1:T*(find(list.opt=='H'))) =opt_all(:,list.allvars=='H'); % as starting value use hh
   
        x0(T*(find(list.opt=='C')-1)+1:T*(find(list.opt=='C')))     =opt_all(:,list.allvars=='C');   % C
        x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F')))     =kappaa'.*opt_all(:,list.allvars=='F'); %0.999*[Ftarget(1); Ftarget(1); Ftarget];%0.8999*0.1066/0.1159*opt_all(:,list.allvars=='F');
        x0(T*(find(list.opt=='G')-1)+1:T*(find(list.opt=='G')))     =opt_all(:,list.allvars=='G');   % G
        
   if indic.xgrowth==0
        x0(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff')))  =opt_all(:,list.allvars=='sff');  % Af
        x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))   =opt_all(:,list.allvars=='sg');  % Ag
        x0(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn')))   =opt_all(:,list.allvars=='sn');  % An
        % x0(T*(find(list.opt=='wsg')-1)+1:T*(find(list.opt=='wsg')))   =opt_all(:,list.allvars=='ws');  % An
   end  

elseif indic.target==0

    x0(T*(find(list.opt=='Lf')-1)+1:T*(find(list.opt=='Lf'))) =opt_all(:,list.allvars=='Lf'); % hhf; first period in LF is baseline
    x0(T*(find(list.opt=='Lg')-1)+1:T*(find(list.opt=='Lg'))) =opt_all(:,list.allvars=='Lg'); % hhg
    x0(T*(find(list.opt=='H')-1)+1:T*(find(list.opt=='H'))) =opt_all(:,list.allvars=='H'); % as starting value use hh
    x0(T*(find(list.opt=='C')-1)+1:T*(find(list.opt=='C')))     =opt_all(:,list.allvars=='C');   % C
   
    x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F')))     =opt_all(:,list.allvars=='F');
    x0(T*(find(list.opt=='G')-1)+1:T*(find(list.opt=='G')))     =opt_all(:,list.allvars=='G');   % G
    
    if indic.xgrowth==0
        x0(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff')))   =opt_all(:,list.allvars=='sff');  % Af
        x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))   =opt_all(:,list.allvars=='sg');  % Ag
        x0(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn')))   =opt_all(:,list.allvars=='sn');  % An
    end 
end

  
%%
% Transform to unbounded variables %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- most of variables bounded by zero below
guess_trans=log(x0);

guess_trans(T*(find(list.opt=='H')-1)+1:T*(find(list.opt=='H')))=log((params(list.params=='upbarH')-x0(T*(find(list.opt=='H')-1)+1:T*(find(list.opt=='H'))))./...
     x0(T*(find(list.opt=='H')-1)+1:T*(find(list.opt=='H'))));

if indic.target==1
     guess_trans(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F')))=log((Ftarget-x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F'))))./...
     x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F'))));
end

lb=[];
ub=[];

%%
% Test Constraints and Objective Function %%%
% percon=0;
f =  objective(guess_trans, T, params, list, Ftarget, indic, init201519, percon, MOM, taulFixed);
[c, ceq] = constraints(guess_trans, T, params, init201519, list, Ems, indic, MOM, percon, taulFixed);

% to examine stuff
% ind=1:length(ceq);
%  ss=ind(abs(ceq)>1e-6);
%  tt=floor(ss/T)+1; 

%%% Optimize %%%
%%%%%%%%%%%%%%%%
%Note: Active-set algorithm is benchmark and assumed for calculations below 
%using Lagrange Multipliers. However, for convergence and accuracy, may
%first (need to) run scenario in interior-point (non-specified) and/or sqp algorithm,
%utilize results to generate initial point and then re-run active-set algorithm.

constf=@(x)constraints(x, T, params, init201519, list, Ems, indic, MOM, percon, taulFixed);
objf=@(x)objective(x, T, params, list, Ftarget, indic, init201519, percon, MOM, taulFixed);

%%
    options = optimset('algorithm','active-set','TolCon',1e-10,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
    [x,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constf,options);
    options = optimset('algorithm','active-set','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
    [x,fval,exitflag,output,lambda] = fmincon(objf,x,[],[],[],[],lb,ub,constf,options);
    
    save(sprintf('../output/2407_OPT_notaul%d_subsres%d_target%d', indic.notaul, indic.subsres, indic.target))
    %xx=load(sprintf('../output/2407_OPT_notaul%d_subsres%d_target%d', indic.notaul, indic.subsres,  indic.target));

    %%
if indic.testT==1
    % code to test number of periods for direct optimization
    Test_lengthOPt;
end
%% transform
% helper=load(sprintf('active_set_solu_notargetOPT_505_spillover%d_taus%d_possible', indic.spillovers, indic.taus))
% x=xsqp;

out_trans=exp(x);
    out_trans((find(list.opt=='H')-1)*T+1:find(list.opt=='H')*T)=upbarH./(1+exp(x((find(list.opt=='H')-1)*T+1:find(list.opt=='H')*T)));
if indic.target==1
    out_trans((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)=Ftarget./(1+exp(x((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)));
end

[xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, H, A_lag, Emnet, A,muu,...
            pn, pg, pf, pee,  wsf, wsn, wsg,  tauf, taul, taus, taurese, tauresg, Trans,...
            w, SWF, PV,PVSWF, objF]= aux_OPT(out_trans, list, params, T, init201519, indic, MOM, taulFixed);

ws = wsg; 
gammal = zeros(size(pn));
obs =[PV, PVSWF, objF]; % save measures of objective function 

opt_all=eval(symms.allvars);

%% Test if optimal allocation is solution to Competitive eqbm

helper.LF_SIM=opt_all';
indic.limit_LF=0; % for testing no constraint on tauf
%list.allvars=listsave.allvars;
test_OPT(T, list,  params,symms, init201519, helper, indic, Ems);

%% save results
   if indic.target==1
        save(sprintf('../output/OPT_target_2407_emnet%d_notaul%d_subsres%d_phii%d_sigma%d_Bop%d_util%d_PV%d.mat',...
            indic.targetWhat, indic.notaul, indic.subsres, indic.know_spill, indic.sigmaWorker, indic.Bop, indic.util), 'opt_all', 'obs')
        %fprintf('saved')
    else
        save(sprintf('../output/OPT_notarget_2407_notaul%d_subsres%d_phii%d_sigma%d_Bop%d_util%d_PV%d.mat',...
            indic.notaul, indic.subsres, indic.know_spill, indic.sigmaWorker, indic.Bop, indic.util), 'opt_all', 'obs')
        
    end
end