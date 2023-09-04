function [LF_SIM, poll, FVAL, indexx] = solve_LF(T, list, poll, params, symms, x0LF, init, indexx, indic, Ems)
% simulate economy under laissez faire
% input: 
% pol: numeric policy vector

% output
% LF_SIM: matrix of simulated results in LF, rows= time period, column =
% variables as ordered in list.allvars

%-- to save results
LF_SIM=zeros(length(list.allvars),T); 
FVAL  = zeros(T,1);

%-- initialise values
x0      = x0LF; % initial guess from calibration 
laggs   = init; % (init should refer to 2010-2014 period)
t       = 1; % number of periods: t=1: 2015-2019 => does include base year period (in matrix on first row) but dont save!

%%
%-- check size of policy matrix
[row]=size(poll);

%- meeting emission limit: extend set of choice variables, adjust indexx
if indic.limit_LF==1
    
    syms tauf real
    symms.choice = [symms.choice, 'tauf'];
    
    list.choice = string(symms.choice);
    x0=[x0, 2];
    hhelper=indexx('LF');
    hhelper.lab = [hhelper.lab, boolean(zeros(1,1))];
    hhelper.exp = [hhelper.exp, boolean(zeros(1,1))];
    hhelper.sqr = [hhelper.sqr, boolean(zeros(1,1))];
    indexx('LF')=hhelper;
end
%%
while t<=T+1 % because first iteration is base year
    fprintf('entering simulation of period %d', t);
    %-- emission limit (=0 in base year)
    if t==1 % baseline period: no emission limit, take same as in second period
        Emlim=0;
    else % 
        Emlim = Ems(t-1);
    end
    
    %--- read in policy (can be vector or static)
    if row(1)>1
        if t<=T
            pol=poll(t,:);
        else
            pol=poll(t-1,:);
        end
    else
        pol=poll;
    end
    %% - transforming variables to unbounded variables
    %-- index for transformation 
    guess_trans=trans_guess(indexx('LF'), x0, params, list.params);
    
    % test
    f=laissez_faire(guess_trans, params, list, pol, laggs, indic, Emlim, t);
    %% - solving model
     % lb=[];
     % ub=[];
     % if ~(indic.limit_LF==0 && indic.notaul==4 && indic.noknow_spill==0)
     %     objf=@(x)objectiveFAKE(x);
     %     constrf = @(x)laissez_faire_nows_fmincon_sep(x, params, list, pol, laggs, indic, Emlim, t);
     % 
     %     if indic.labshareequ==0 && (~indic.taul0==1)
     %        options = optimset('algorithm','active-set','TolCon', 1e-7,'Tolfun',1e-26,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
     %     else
     %        options = optimset('algorithm','sqp','TolCon', 1e-8,'Tolfun',1e-26,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
     %     end
     %    [sol3,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constrf,options);
     % end
  
%- other solvers
    modFF = @(x)laissez_faire(x, params, list, pol, laggs, indic, Emlim, t);
    options = optimoptions('fsolve', 'TolFun', 10e-9,'Display','iter', 'MaxFunEvals',8e6, 'MaxIter', 3e6,  'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
    [sol3, fval, exitf] = fsolve(modFF, guess_trans, options);

    % pass to standard algorithm
%      options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5,);%, );%, );%, 'Display', 'Iter', );
%     [sol, fval, exitf] = fsolve(modFF, x1, options);
% 
% if ~(indic.noskill==1 && indic.tauf==1 && indic.xgrowth==0)
% if indic.sep<=1
%     if ~(indic.sizeequ==1 && indic.GOV==0 && indic.noknow_spill==0 ) &&  (indic.labshareequ==1&& indic.GOV==0 && indic.noknow_spill==0 )
%          options = optimoptions('fsolve', 'TolFun', 10e-9, 'MaxFunEvals',8e3, 'MaxIter', 3e5);% 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
%          [sol3, fval, exitf] = fsolve(modFF, sol2, options);
%     else
%         sol3=sol2;
%     end
% else
%     sol3=sol2;
% end
% 
if exitf<=0
    error('code did not solve')
end
%- transform results to bounded variables
LF=trans_allo_out(indexx('LF'), sol3, params, list.params, indic);

%% - save results
    % this part also checks correctness of results!
    cell_par=arrayfun(@char, symms.choice, 'uniform', 0);
    SLF=cell2struct(num2cell(LF), cell_par, 2);
    if t>1
%         if indic.sep==0
%             LF_SIM(:,t-1)=aux_solutionLF( SLF, pol, laggs, list, symms, indexx, params, indic);
%         else
            LF_SIM(:,t-1)=aux_solutionLF( SLF, pol, laggs, list, symms, indexx, params, indic, Emlim,t);
%         end
        FVAL(t-1)=max(abs(fval));
    end
    %% - update for next round
    x0 = LF; % initial guess
        Af0= SLF.Af; % today's technology becomes tomorrow's lagged technology
        Ag0= SLF.Ag; 
        An0= SLF.An; 
    laggs=eval(symms.init);
    t=t+1;
end
end