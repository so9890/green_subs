function [LF_SIM]=solve_LF_VECT(T, list, params,symms, init201519, helper, indic, Ems, MOM)
%test OPT policy result without target in competitive equilibrium
read_in_params;

varrs=helper.LF_SIM; 
y=log(varrs);
z=sqrt(varrs);

% read in results

H= log((params(list.params=='upbarH')-varrs(list.allvars=='H', :))./(varrs(list.allvars=='H', :)))';
w=y(list.allvars=='w', :)';
Lf=y(list.allvars=='Lf', :)';
Ln=y(list.allvars=='Ln', :)';
Lg=y(list.allvars=='Lg', :)';
C=y(list.allvars=='C', :)';
F=y(list.allvars=='F', :)';
G=y(list.allvars=='G', :)';

Af = y(list.allvars=='Af', :)';
Ag = y(list.allvars=='Ag', :)';
An =y(list.allvars=='An', :)';

sff =z(list.allvars=='sff', :)';
sg =z(list.allvars=='sg', :)';
sn =z(list.allvars=='sn', :)';
ws=z(list.allvars=='ws', :)';
gammal =z(list.allvars=='gammal', :)';
pg=y(list.allvars=='pg', :)';
pn=y(list.allvars=='pn', :)';
pee=y(list.allvars=='pee', :)';
pf=y(list.allvars=='pf', :)';

Trans=varrs(list.allvars=='Trans', :)';

if indic.limit_LF==1
    syms tauf real
    symms.choice=[symms.choice, tauf];
    tauf=varrs(list.allvars=='tauf',:)';    
end
list.choice=string(symms.choice);

x0=eval(symms.choice);
x0=x0(:);

%f=laissez_faire_VECT(x0, params, list, varrs, init201519,T, indic, Ems);
% 
 % if max(abs(f))>1e-9
%    % equations where results are off
     % ind=1:length(f);
     % pos=ind(abs(f)>1e-9);
     % eqset= floor(pos./T);
%     error('optimal policy does not solve laissez faire. For the relevant equations see eqset');
%     % to evaluate stuff
 % end

modFF = @(x)laissez_faire_VECT(x, params, list, varrs, init201519,T, indic, Ems);
options = optimoptions('fsolve', 'TolFun', 10e-10, 'MaxFunEvals',8e3, 'MaxIter', 3e5,  'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[x0, fval, exitf] = fsolve(modFF, x0, options);

% options = optimoptions('fsolve', 'TolFun', 10e-6,'Display','iter', 'MaxFunEvals',8e3, 'MaxIter', 3e5); %,  'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
% [x0, fval, exitf] = fsolve(modFF, x0, options);

lb=[];
ub=[];

objf=@(x)objectiveFAKE(x);
    % if indic.xgrowth==0
%         if indic.sep>=1
            constLF=@(x)laissez_faire_VECT_fmincon(x, params, list, varrs, init201519,T, indic, Ems);
% %         else
% %             constLF=@(x)laissez_faireVECT_fmincon(x, params, list, varrs, init201519,T, indic);
% %         end
%     else
%         constLF=@(x)laissez_faireVECT_xgrowth_fmincon(x, params, list, varrs, init201519, T, indic, MOM);
%     end
% 
%     if (indic.xgrowth==0 && indic.noskill==1) 
        options = optimset('algorithm','sqp','TolCon', 1e-9,'Tolfun',1e-26,'MaxFunEvals',5e10,'MaxIter',6e10,'Display','iter');
    % else
    %     options = optimset('algorithm','active-set','TolCon', 1e-8,'Tolfun',1e-26,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
    % end
[x,fval,exitflag,output,lambda] = fmincon(objf,x0,[],[],[],[],lb,ub,constLF,options);

%  count=0;
%  while exitflag==-2 && count<4
%      count=count+1
% [x,fval,exitflag,output,lambda] = fmincon(objf,x,[],[],[],[],lb,ub,constLF,options);
% end
% test solution to 
% if indic.xgrowth==0
%     if indic.sep==0
%         f=laissez_faireVECT(x, params, list, varrs, init201519,T, indic);
%     else
f=laissez_faire_VECT(x, params, list, varrs, init201519,T, indic, Ems);
%     end
% else
%     f=laissez_faireVECT_xgrowth(x, params, list, varrs, init201519, T, indic);
% end

if max(abs(f))>1e-9
    error('LF function does not solve')
end

% save results
% if indic.xgrowth==0
% if indic.sep==0
%     LF_SIM=aux_solutionLF_VECT(x, list, symms, varrs, params, T, indic);
% else
    LF_SIM=aux_solutionLF_VECT(x, list, symms, varrs, params, T, indic);
% end
% else
%     LF_SIM=aux_solutionLF_VECT_xgrowth(x, list, symms,varrs, params, T , indic, init201519);
% end
end
