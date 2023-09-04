function [symms, list, sp_all]=SP_solve(list, symms, params, Sparams, x0LF, init201014, init201519, indexx, indic, T, Ems, MOM, percon)

% pars
read_in_params;

% function to find social planner allocation 
syms H xn xf xg An Af Ag C F sff sg sn Lg Ln real

symms.sp = [Lg Ln xn xf xg An Af Ag C H F sff sg sn];
list.sp  = string(symms.sp); 
nn= length(list.sp); 

%%
%%%%%%%%%%%%%%%%%%%%%%
% Initial Guess %
%%%%%%%%%%%%%%%%%%%%%

% laissez faire
% if indic.target==0
% helper=load(sprintf('../output/LF_phii%d_sigma%d_Bop%d_util%d.mat', ...
%      indic.know_spill, indic.sigmaWorker, indic.Bop));
% else
%     helper=load(sprintf('../output/Limit_phii%d_sigma%d_Bop%d_util%d.mat', ...
%      indic.know_spill, indic.sigmaWorker, indic.Bop));
% end
% sp_all=helper.COMP';

if indic.target==1
%     helper=load(sprintf('../output/OPT_target_2407_emnet%d_notaul%d_subsres%d_phii%d_sigma%d_Bop%d_util%d_PV%d.mat',...
%        indic.targetWhat, 0, indic.subsres, indic.know_spill, indic.sigmaWorker, indic.Bop, indic.util), 'opt_all', 'obs');

    helper=load(sprintf('../output/SP_target_2407_emnet%d_phii%d_sigma%d_Bop%d_util%d_PV%d.mat',...
        indic.targetWhat, indic.know_spill, indic.sigmaWorker, indic.Bop, indic.util), 'sp_all');
    %fprintf('saved')
else
       % helper=load(sprintf('../output/OPT_notarget_2407_notaul%d_subsres%d_phii%d_sigma%d_Bop%d_util%d_PV%d.mat',...
       %      0, indic.subsres, indic.know_spill, indic.sigmaWorker, indic.Bop, indic.util), 'opt_all', 'obs');
 
   helper=load(sprintf('../output/SP_notarget_2407_phii%d_sigma%d_Bop%d_util%d_PV%d.mat',...
        indic.know_spill, indic.sigmaWorker, indic.Bop, indic.util), 'sp_all');
    
end
% sp_all = helper.opt_all;
sp_all=helper.sp_all;
%list.allvars=list.allvars(list.allvars~='taurese'& list.allvars~='tauresg'& list.allvars~='taus');
%%
if indic.target==0
    x0 = zeros(nn*T,1);
    Ftarget = 0; % placeholder

    if ~isfile(sprintf('../output/SP_notarget'))
        
        x0(T*(find(list.sp=='Lg')-1)+1:T*(find(list.sp=='Lg'))) =sp_all(1:T,list.allvars=='Lg'); % hlg 
        x0(T*(find(list.sp=='Ln')-1)+1:T*(find(list.sp=='Ln'))) =sp_all(1:T,list.allvars=='Ln'); % hlg 
        x0(T*(find(list.sp=='H')-1)+1:T*(find(list.sp=='H')))   =sp_all(1:T,list.allvars=='H');  % hh
    
        x0(T*(find(list.sp=='xf')-1)+1:T*(find(list.sp=='xf')))   =sp_all(1:T, list.allvars=='xf'); % hlf
        x0(T*(find(list.sp=='xg')-1)+1:T*(find(list.sp=='xg')))   =sp_all(1:T,list.allvars=='xg'); % hlg 
        x0(T*(find(list.sp=='xn')-1)+1:T*(find(list.sp=='xn')))   =sp_all(1:T,list.allvars=='xn'); % hlg 
        x0(T*(find(list.sp=='Af')-1)+1:T*(find(list.sp=='Af')))   =sp_all(1:T,list.allvars=='Af');  % Af
        x0(T*(find(list.sp=='Ag')-1)+1:T*(find(list.sp=='Ag')))   =sp_all(1:T,list.allvars=='Ag');  % Ag
        x0(T*(find(list.sp=='An')-1)+1:T*(find(list.sp=='An')))   =sp_all(1:T,list.allvars=='An');  % An

        x0(T*(find(list.sp=='sn')-1)+1:T*(find(list.sp=='sn')))   =sp_all(1:T,list.allvars=='sn');  % C
        x0(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg')))   =sp_all(1:T,list.allvars=='sg');  % C
        x0(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff'))) =sp_all(1:T,list.allvars=='sff');  % C
        
        x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))     =sp_all(1:T,list.allvars=='C');  % C
        x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))     =sp_all(1:T,list.allvars=='F');  % C
    
    else
       fprintf('using sp solution') 

      
       %p_all=helper.sp_all;
   
            x0(T*(find(list.sp=='Lg')-1)+1:T*(find(list.sp=='Lg')))   =sp_all(1:T, list.allvars=='Lg'); % hlg
            x0(T*(find(list.sp=='Ln')-1)+1:T*(find(list.sp=='Ln')))   =sp_all(1:T, list.allvars=='Ln'); % hlg
            x0(T*(find(list.sp=='H')-1)+1:T*(find(list.sp=='H')))     =sp_all(1:T, list.allvars=='hh');  % hh
          
            x0(T*(find(list.sp=='xf')-1)+1:T*(find(list.sp=='xf')))   =sp_all(1:T, list.allvars=='xf'); % hlf
            x0(T*(find(list.sp=='xg')-1)+1:T*(find(list.sp=='xg')))   =sp_all(1:T, list.allvars=='xg'); % hlg 
            x0(T*(find(list.sp=='xn')-1)+1:T*(find(list.sp=='xn')))   =sp_all(1:T, list.allvars=='xn'); % hlg 
            
            x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))     =sp_all(1:T, list.allvars=='C');  %             
            x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))     =sp_all(1:T, list.allvars=='F');  % 
            if indic.xgrowth==0
                x0(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg')))   =sp_all(1:T,list.allvars=='sg');  % C
                x0(T*(find(list.sp=='sn')-1)+1:T*(find(list.sp=='sn')))   =sp_all(1:T, list.allvars=='sn');  % C
                x0(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff'))) =sp_all(1:T, list.allvars=='sff');  % C

                x0(T*(find(list.sp=='Af')-1)+1:T*(find(list.sp=='Af')))   =sp_all(1:T, list.allvars=='Af');  % Af
                x0(T*(find(list.sp=='Ag')-1)+1:T*(find(list.sp=='Ag')))   =sp_all(1:T, list.allvars=='Ag');  % Ag
                x0(T*(find(list.sp=='An')-1)+1:T*(find(list.sp=='An')))   =sp_all(1:T, list.allvars=='An');  % An
            end
    end
    
elseif indic.target==1 
    
    Ftarget = (Ems+deltaa)/omegaa;
    x0 = zeros(nn*T,1);
        %fprintf('using sp solution as initial value')
        % with new emission target
        kappaa = [repmat(Ftarget(1),1, percon),Ftarget]./sp_all(1:T,list.allvars=='F')'; % ratio of targeted F to non-emission
        kappaa = kappaa*(1-1e-10);
      
        x0(T*(find(list.sp=='Lg')-1)+1:T*(find(list.sp=='Lg')))   =sp_all(1:T, list.allvars=='Lg'); % hlg
        x0(T*(find(list.sp=='Ln')-1)+1:T*(find(list.sp=='Ln')))   =sp_all(1:T, list.allvars=='Ln'); % hlg
        x0(T*(find(list.sp=='H')-1)+1:T*(find(list.sp=='H')))     =sp_all(1:T, list.allvars=='H');  % hh
    
        x0(T*(find(list.sp=='xf')-1)+1:T*(find(list.sp=='xf')))   =sp_all(1:T, list.allvars=='xf'); % hlf
        x0(T*(find(list.sp=='xg')-1)+1:T*(find(list.sp=='xg')))   =sp_all(1:T, list.allvars=='xg'); % hlg 
        x0(T*(find(list.sp=='xn')-1)+1:T*(find(list.sp=='xn')))   =sp_all(1:T, list.allvars=='xn'); % hlg
        x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))     =sp_all(1:T, list.allvars=='C');  % C

        x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))     =kappaa.*sp_all(1:T, list.allvars=='F')';  % C
        if indic.xgrowth==0
            x0(T*(find(list.sp=='Af')-1)+1:T*(find(list.sp=='Af')))   =sp_all(1:T, list.allvars=='Af');  % Af
            x0(T*(find(list.sp=='Ag')-1)+1:T*(find(list.sp=='Ag')))   =sp_all(1:T, list.allvars=='Ag');  % Ag
            x0(T*(find(list.sp=='An')-1)+1:T*(find(list.sp=='An')))   =sp_all(1:T, list.allvars=='An');  % An
            x0(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg')))   =sp_all(1:T, list.allvars=='sg');  % C
            x0(T*(find(list.sp=='sn')-1)+1:T*(find(list.sp=='sn')))   =sp_all(1:T, list.allvars=='sn');  % C
            x0(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff'))) =sp_all(1:T, list.allvars=='sff');  % C
        end

%     end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial values for An0, Ag0, Af0 refer to 2015-2019=> first period in
% Laissez faire solution; do not take init which refers to 2010-2014!
% this version here under assumption of first best policy in initial policy
% that is: tauf=0, taul=0, taus=0, lambdaa to balage budget
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initOPT= init201519; % as calibrated under BAU policy
% Transform variables to unbounded vars => requires less constraints! %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 guess_trans=log(x0);
  guess_trans(T*(find(list.sp=='H')-1)+1:T*(find(list.sp=='H')))=log((params(list.params=='upbarH')-x0(T*(find(list.sp=='H')-1)+1:T*(find(list.sp=='H'))))./...
  x0(T*(find(list.sp=='H')-1)+1:T*(find(list.sp=='H'))));
 
     % guess_trans(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff')))=sqrt(x0(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff'))));
     % guess_trans(T*(find(list.sp=='sn')-1)+1:T*(find(list.sp=='sn')))=sqrt(x0(T*(find(list.sp=='sn')-1)+1:T*(find(list.sp=='sn'))));
     % guess_trans(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg')))=sqrt(x0(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg'))));
 
if indic.target==1
    % only from the third period onwards F is contrained
    
    guess_trans(T*(find(list.sp=='F')-1)+1+percon:T*(find(list.sp=='F')))=log((Ftarget'-x0(T*(find(list.sp=='F')-1)+1+percon:T*(find(list.sp=='F'))))./...
     x0(T*(find(list.sp=='F')-1)+1+percon:T*(find(list.sp=='F'))));
end
lb=[];
ub=[];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Constraints and Objective Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f =  objectiveSP(guess_trans,T,params, list, Ftarget, indic, initOPT, percon);
[c, ceq] = constraintsSP(guess_trans, T, params, initOPT, list, Ems, indic, percon, MOM);

objfSP=@(x)objectiveSP(x,T,params, list, Ftarget, indic, initOPT, percon);
constfSP=@(x)constraintsSP(x, T, params, initOPT, list, Ems, indic, percon, MOM);

%%
options = optimset('algorithm','sqp', 'TolCon',1e-8, 'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
[x,fval,exitflag,output,lambda] = fmincon(objfSP,guess_trans,[],[],[],[],lb,ub,constfSP,options);
options = optimset('algorithm','active-set','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
[x,fval,exitflag,output,lambda] = fmincon(objfSP,x,[],[],[],[],lb,ub,constfSP,options);
% options = optimset('algorithm','sqp', 'TolCon',1e-8, 'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
% [x,fval,exitflag,output,lambda] = fmincon(objfSP,x,[],[],[],[],lb,ub,constfSP,options);

save(sprintf('../output/2407_SP_target%d', indic.target))

%%
out_trans=exp(x);
out_trans((find(list.sp=='H')-1)*T+1:find(list.sp=='H')*T)=upbarH./(1+exp(x((find(list.sp=='H')-1)*T+1:find(list.sp=='H')*T)));
% out_trans((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T)=(x((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T)).^2;
% out_trans((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T)=(x((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T)).^2;
% out_trans((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T)=(x((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T)).^2;

if indic.target==1
    out_trans((find(list.sp=='F')-1)*T+1+percon:find(list.sp=='F')*T)=Ftarget'./(1+exp(x((find(list.sp=='F')-1)*T+1+percon:find(list.sp=='F')*T)));
end

[xn,xf,xg,Ag, An, Af,...
Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
F, N, G, E, Y, C, H, A_lag, Emnet, A,muu,...
pn, pg, pf, pee, w, wsn, wsf, wsg, tauf, taul,...
taurese, tauresg, ...
Trans, SWF, PV,PVSWF, objF]...
    = aux_SP(out_trans, list, params, T, initOPT, indic);
gammal = zeros(size(pn));
ws=wsn;
taus =zeros(size(ws)) ; % only for now!
obs =[PV, PVSWF, objF]; % save measures of objective function 

sp_all=eval(symms.allvars);

%% save results
    if indic.target==1
        save(sprintf('../output/SP_target_2407_emnet%d_phii%d_sigma%d_Bop%d_util%d_PV%d.mat',...
            indic.targetWhat, indic.know_spill, indic.sigmaWorker, indic.Bop, indic.util), 'sp_all', 'obs')
    else
        save(sprintf('../output/SP_notarget_2407_phii%d_sigma%d_Bop%d_util%d_PV%d.mat',...
             indic.know_spill, indic.sigmaWorker, indic.Bop, indic.util), 'sp_all', 'obs')
        
    end
end
