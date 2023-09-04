% in contrast to "solve easy" this version with tauf on F only
% to solve easy Sp and op problems
addpath('../tools')
% to save results of loop
keySet={'<1', 'log', 'Bop'};
valueSet= repmat({struct([])},1,length(keySet));
resultsTHETA=containers.Map(keySet, valueSet);
indic.taxsch=1; %
                %==0 then uses linear tax schedule with lump sum transfers
                %    of carbon tax revenues (efficient)
                %==1 linear tax, carbon tax revenues used for green
                %    subsidies (politically feasible)

indic.notaul=0; % relevant for optimal solution
indic.util =0; 
indic.extern =1; % ==0 no externality, ==1 with externality

syms taus real
symms.poleasy =[symms.pol, taus];
list.poleasy= string(symms.poleasy);
%% Sp solution 
clear Sp
x0=log([0.4,0.4]);

modFF = @(x)easy_sp(x, params, list, indic, init201519);
options = optimoptions('fsolve', 'TolFun', 10e-8, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt','DiSplay', 'Iter');%, );%, );%,  );
[sol2, fval, exitf] = fsolve(modFF, x0, options);
options = optimoptions('fsolve', 'TolFun', 10e-10, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'DiSplay', 'Iter');%, );%, );%,  );
[sol3, fval, exitf] = fsolve(modFF, sol2, options);

% save solution
read_in_params;
clear Sp;
Sp.Lf=exp(sol3(2));
Sp.Lg=exp(sol3(1));
Sp.Af=(1+vf)*init201519(list.init=='Af0');
Sp.Ag=(1+vg)*init201519(list.init=='Ag0');
Sp.h=Sp.Lg+Sp.Lf;
Sp.F=Sp.Af*Sp.Lf;
Sp.G=Sp.Ag*Sp.Lg;
Sp.Y=(Sp.F)^(eppsy)*(Sp.G)^(1-eppsy);
    deltaYF=eppsy*(Sp.G./Sp.F) ^(1-eppsy);
    deltaYG=(1-eppsy)*(Sp.F./Sp.G)^(eppsy);
Sp.C=Sp.Y;
Sp.pg= (1-eppsy)*(Sp.F./Sp.G)^eppsy;
if indic.taxsch==1 % with subsidy
    Sp.pf= Sp.Lf./(Sp.Lg+Sp.Lf)*(Sp.G/Sp.F+deltaYF/deltaYG).*Sp.pg;% from free movement of labor, gov budget with taus, and optimality Y
else
    Sp.pf= Sp.pg.*Sp.Ag./Sp.Af;
end
Sp.pigou =deltaYF-Sp.pf;
if indic.taxsch==0
    Sp.taus = 0;
else
    Sp.taus = Sp.pigou.*Sp.F./Sp.Lg;
end
Sp.s = Sp.Lf/Sp.h; % share of labour in fossil sector Lf/h
Sp.w= Sp.pf.*Sp.Af; 
Sp.taul = 1-Sp.h^sigmaa.*chii.*Sp.C^(thetaa)./Sp.w; % from hh optimality
%Sp.pg= eppsy^eppsy*(1-eppsy)^(1-eppsy)*((1-eppsy)/eppsy*Sp.Lf*Sp.Af/Sp.Ag/Sp.Lg)^eppsy;
% %- utility
if indic.util==1
    Sp.Ucon=(Sp.C^(1-thetaa)-1)/(1-thetaa);
else
    Sp.Ucon=log(Sp.C);
end
Sp.Ext = -weightext*(omegaa*Sp.F)^extexpp;
Sp.dEdF = -weightext*extexpp*(omegaa*Sp.F)^extexpp/Sp.F;
Sp.Ulab = -chii*Sp.h^(1+sigmaa)/(1+sigmaa);

Sp.SWF = Sp.Ucon+Sp.Ulab+indic.extern*Sp.Ext;
Sp.Ul = -chii*Sp.h^sigmaa;
% test analytic derivation of Sp.h
% Sp.htest= ((Sp.w^(1-thetaa)+Sp.dEdF*Sp.Af*Sp.s*Sp.h^thetaa)/chii)^(1/(sigmaa+thetaa));
Sp.htest = (Sp.w^(1-thetaa)/chii*(1-eppsy)/(1-Sp.s))^(1/(sigmaa+thetaa));

Sp.Ultest = -Sp.C^(-thetaa)*Sp.Y/Sp.G*(1-eppsy)*Sp.Ag; % marginal disutility of labour 
%% Competitive equilibrium

clear LF
x0=log([Sp.pg,Sp.Lg]);
LF.tauf=Sp.pigou;
LF.taul=Sp.taul; % this tax implements efficient hours worked given tauf=Sp.pigou
LF.taus=Sp.taus; % depends on tax scheme!
tauf=LF.tauf;
taul=LF.taul;
pol=eval(symms.pol);
modFF = @(x)easy_lf_taus(x, params, list, pol,  init201519, indic);
options = optimoptions('fsolve', 'TolFun', 10e-8, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt','DiSplay', 'Iter');%, );%, );%,  );
[sol2, fval, exitf] = fsolve(modFF, x0, options);
options = optimoptions('fsolve', 'TolFun', 10e-10, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'DiSplay', 'Iter');%, );%, );%,  );
[sol3, fval, exitf] = fsolve(modFF, sol2, options);

% save solution
read_in_params;

LF.Af=(1+vf)*init201519(list.init=='Af0');
LF.Ag=(1+vg)*init201519(list.init=='Ag0');
LF.Lf=exp(sol3(1));
LF.Lg=exp(sol3(2));

LF.G = LF.Ag*LF.Lg;
LF.F = LF.Af*LF.Lf;
LF.pg = (1-eppsy)*(LF.F./LF.G)^eppsy; 
LF.w = LF.pg*LF.Ag+LF.taus;
LF.pf = LF.w/(LF.Af);
LF.Y = (LF.F)^(eppsy)*(LF.G)^(1-eppsy);
LF.h = LF.Lf+LF.Lg;
LF.s =LF.Lf/LF.h;

if indic.taxsch==1 % subsidies
    LF.C= LF.w*LF.h;
    LF.T=LF.w*LF.h.*(LF.taul);
    LF.hsup = (LF.w.^(1-thetaa)*(1-LF.taul)./chii).^(1/(thetaa+sigmaa));

else % lump sum tansfers of tf*F
    LF.T=LF.w*LF.h.*(LF.taul)+LF.tauf.*LF.F;
    LF.C= LF.w*LF.h.*(1-LF.taul)+LF.T;
    LF.hsup = ((((1-LF.taul).*LF.w+LF.T./LF.h)^(-thetaa)*LF.w*(1-LF.taul))/(chii))^(1/(sigmaa+thetaa));
end

% %- utility
if indic.util==1
    LF.Ucon=(LF.C^(1-thetaa)-1)/(1-thetaa);
else
    LF.Ucon=log(LF.C);
end
LF.Ext = -weightext*(omegaa*LF.F)^extexpp;
LF.dEdF = LF.Ext*extexpp/(LF.F); % negative!
LF.Ulab = -chii*LF.h^(1+sigmaa)/(1+sigmaa);

LF.SWF = LF.Ucon+LF.Ulab+indic.extern*LF.Ext;

% social cost of carbon as: what would a household be willing to pay
LF.scc = -LF.dEdF*LF.C.^thetaa;

% analytic derivation of transfers so that competitive eqbm = First Best
% only if utility = log
% if indic.util==0
%    LF.Tana=-LF.dEdF*LF.Af*LF.s*LF.h^thetaa/(1+LF.dEdF*LF.Af*LF.s*LF.h^thetaa)*LF.w*LF.h;
%    LF.Tana_gov= LF.tauf/(1-LF.tauf)*LF.w*LF.s*LF.h;
%    LF.taufcheckk=-LF.dEdF*LF.Af*LF.h/(1-LF.dEdF*LF.Af*LF.h*(1-LF.s));
% end

%% optimal policy

clear Opt
if indic.notaul==0
    x0=log([Sp.s,Sp.h]);
else
    x0=log([Sp.s]);
end
% if indic.taxsch<=1
    modFF = @(x)easy_opt_taus(x, params, list,  init201519, indic);
% else
%     modFF = @(x)easy_opt_Gov(x, params, list,  init201519, indic);
% end
options = optimoptions('fsolve', 'TolFun', 10e-8, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt','DiSplay', 'Iter');%, );%, );%,  );
[sol2, fval, exitf] = fsolve(modFF, x0, options);
options = optimoptions('fsolve', 'TolFun', 10e-10, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'DiSplay', 'Iter');%, );%, );%,  );
[sol3, fval, exitf] = fsolve(modFF, sol2, options);

%% save solution 

read_in_params;

Opt.Af=(1+vf)*init201519(list.init=='Af0');
Opt.Ag=(1+vg)*init201519(list.init=='Ag0');

Opt.s = exp(sol3(1));
if indic.notaul==0
    Opt.h = exp(sol3(2));
% auxiliary
Opt.Lg = (1-Opt.s)*Opt.h; 
    Opt.Lf = Opt.s*Opt.h;
    % production
    Opt.G = Opt.Ag*Opt.Lg;
    Opt.F = Opt.Af*Opt.Lf;
    Opt.Y=(Opt.F)^(eppsy)*(Opt.G)^(1-eppsy);
        deltaYF=eppsy*(Opt.G./Opt.F) ^(1-eppsy);
        deltaYG=(1-eppsy)*(Opt.F./Opt.G)^(eppsy);
    Opt.C=Opt.Y;
    Opt.pg= (1-eppsy)*(Opt.F./Opt.G)^eppsy;

    if indic.taxsch==1 % with subsidy
        Opt.pf= Opt.Lf./(Opt.Lg+Opt.Lf)*(Opt.G/Opt.F+deltaYF/deltaYG).*Opt.pg;% from free movement of labor, gov budget with taus, and optimality Y
    else
        Opt.pf= Opt.pg.*Opt.Ag./Opt.Af;
    end
        Opt.tauf =deltaYF-Opt.pf;
    if indic.taxsch==0
        Opt.taus = 0;
    else
        Opt.taus = Opt.tauf.*Opt.F./Opt.Lg;
    end
    Opt.w = Opt.pf.*Opt.Af; 
    Opt.taul = 1-Opt.h^sigmaa.*chii.*Opt.C^(thetaa)./Opt.w; % from hh optimality
    if indic.taxsch ==0
         Opt.T = Opt.tauf.*Opt.F+Opt.taul.*Opt.w.*Opt.h;
    else
         Opt.T = Opt.taul.*Opt.w.*Opt.h;
    end
elseif indic.notaul==1
        FG = Opt.Af*Opt.s/((1-Opt.s)*Opt.Ag);
        Opt.pg = (1-eppsy).*(FG)^eppsy; % optimal demand green
   
        deltaYF=eppsy*(1./FG) ^(1-eppsy);
        deltaYG=(1-eppsy)*(FG)^(eppsy);

    if indic.taxsch==1 % with taus
        Opt.pf = (1./FG+deltaYF./deltaYG).*Opt.pg.*Opt.s;  
    else
        Opt.pf = Opt.Ag.*Opt.pg./Opt.Af;
    end
    Opt.tauf = deltaYF-Opt.pf;
    if indic.taxsch==1
        Opt.taus = Opt.tauf.*Opt.Af.*Opt.s./(1-Opt.s);
        Th=0;
    else
        Opt.taus=0;
        Th=Opt.tauf.*Opt.Af.*Opt.s; % transfers per h
    end
    Opt.w = Opt.pf.*Opt.Af;
    Opt.h = (Opt.w/(chii*(Opt.w+Th)^thetaa))^(1/(sigmaa+thetaa));
    Opt.T = Th*Opt.h;
    Opt.C = Opt.w*Opt.h+Opt.T;
    Opt.Lg = (1-Opt.s).*Opt.h;
    Opt.Lf = Opt.s.*Opt.h;
    Opt.G = Opt.Ag.*Opt.Lg;
    Opt.F = Opt.Af.*Opt.Lf;
    Opt.Y = (Opt.F)^(eppsy)*(Opt.G)^(1-eppsy);
end

% %- utility
if indic.util==1
    Opt.Ucon=(Opt.C^(1-thetaa)-1)/(1-thetaa);
else
    Opt.Ucon=log(Opt.C);
end
Opt.dUcondC = Opt.C^(-thetaa);
Opt.Ext = -weightext*(omegaa*Opt.F)^extexpp;
Opt.dEdF = -weightext*extexpp*(omegaa*Opt.F)^extexpp/Opt.F;
Opt.Ulab = -chii*Opt.h^(1+sigmaa)/(1+sigmaa);

Opt.SWF = Opt.Ucon+Opt.Ulab+indic.extern*Opt.Ext;
Opt.scc = -Opt.dEdF*Opt.C^thetaa;

%% Comparison to analytic results
% check derivative Gov
    
    % derivatives
    dFdh = Opt.Af*Opt.s;
    dFds = Opt.Af*Opt.h;
    dCdh = (Opt.Af*Opt.s)^(eppsy)*(Opt.Ag*(1-Opt.s))^(1-eppsy);
    MPL  = dCdh;
    dCds = dCdh*Opt.h*(eppsy/Opt.s-(1-eppsy)/(1-Opt.s));
    dwdH = 0;
    dwds = (1-eppsy)*Opt.Af^eppsy*Opt.Ag^(1-eppsy)*eppsy*(Opt.s/(1-Opt.s))^(eppsy-1)/(1-Opt.s)^2;
    if indic.taxsch>1
        dGovdh=Opt.tauf*Opt.pf*dFdh;
        dCdh =dCdh-dGovdh;
        dtaufds = -(1-eppsy)/eppsy*(1/(1-Opt.s)^2);
        dpfds =-Opt.pf*(eppsy-1)/(1-Opt.tauf)*dtaufds;
        dGovds = Opt.pf*Opt.F*dtaufds+Opt.tauf*Opt.F*dpfds+Opt.tauf*Opt.pf*dFds;
        dCds = dCds-dGovds;
    % analytic solution to check
    dGovdscheck= Opt.pf*Opt.F*(-(1-eppsy)/eppsy/(1-Opt.s)^2 *(eppsy/Opt.s)+(eppsy-Opt.s)/((1-Opt.s)*eppsy*Opt.s)); %correct!
    taulcheck= (dGovds*Opt.s/Opt.h-dGovdh)/((1-eppsy)*Opt.Y/Opt.G*(-1)*Opt.Ag); 
    taulcheck2= 1+(Opt.h*dCdh-Opt.s*dCds)/(-Opt.h*Opt.w);
    taulcc=dwds*(Opt.s/Opt.w)-dwdH*(Opt.h/Opt.w); 
    % check tauf
    taufcheckk = Opt.scc+dGovds/(Opt.w*Opt.h)*(1-Opt.tauf); 
    taufcheck = Opt.scc+ Opt.tauf-(1-Opt.tauf)/Opt.w*dwds;
    tauff= 1-Opt.scc*Opt.w/dwds; 
    Uh = -chii*Opt.h^(sigmaa);
    d2YdG2=-(1-eppsy)*Opt.Y/Opt.G^2+(1-eppsy)^2*Opt.Y/Opt.G^2;
    taullcheckk= -Opt.h/Opt.w*Opt.Ag^2*d2YdG2;
    end
    
%% save results
%- create structure
st.Sp=Sp;
st.Opt=Opt;
st.LF=LF;
st.thetaa = thetaa;
%- save to map
if indic.util==0
    resultsTHETA('log')= st;
elseif indic.util==1
    if indic.Bop==0
        resultsTHETA('<1')= st;
    else
        resultsTHETA('Bop')= st;

    end
end



% - create table from results 
kk=keys(resultsTHETA);
Table=table(keys(resultsTHETA)',zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1));
Table.Properties.VariableNames={'Thetaa','FB hours', 'FB SWF', 'FB s', 'FB MPL','FB Pigou', 'Only tauf=pigou hours' , 'Only tauf=pigou SWF' , 'Only tauf=pigou s' , 'Only tauf=pigou scc', 'Only tauf=pigou wage', ...
                                       'Optimal hours', 'Optimal SWF','Optimal s', 'Optimal wage', 'Optimal taul', 'Optimal tauf', 'Optimal scc'};

%- only hours, and policy
TableH=table(keys(resultsTHETA)',zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1));
TableH.Properties.VariableNames={'Thetaa','FB hours', 'FB Pigou', 'CE hours',  'CE scc',  ...
                                       'Opt hours','Opt taul', 'Opt tauf', 'Opt scc'};

for i=1:3
    st=resultsTHETA(string(kk(i)));
    Table(i,2:end)={st.Sp.h, st.Sp.SWF,st.Sp.s, st.Sp.w, st.Sp.pigou, st.LF.h, st.LF.SWF, st.LF.s, st.LF.scc, st.LF.w, ...
                    st.Opt.h, st.Opt.SWF, st.Opt.s, st.Opt.w, Opt.taul, st.Opt.tauf, st.Opt.scc};
    TableH(i,2:end)={st.Sp.h, st.Sp.pigou, st.LF.h, st.LF.scc, ...
                    st.Opt.h, Opt.taul, st.Opt.tauf, st.Opt.scc};
end

table2latex(TableH)