function [c, ceq] = constraintsSP(y, T, params, init, list, Ems, indic, percon, MOM)

% pars
read_in_params;
Ftarget=(Ems+deltaa)./omegaa; 
% Ftarg_20s=(MOM.US_Budget20_30+3*deltaa)./omegaa; 
% transform x: all are exponentially transformed
 x=exp(y);
% except for hours

 x((find(list.sp=='H')-1)*T+1:find(list.sp=='H')*T) = upbarH./(1+exp(y((find(list.sp=='H')-1)*T+1:find(list.sp=='H')*T)));

 % x((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T) = (y((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T)).^2;
 % x((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T) = (y((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T)).^2;
 % x((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T) = (y((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T)).^2;

if indic.target==1
 x((find(list.sp=='F')-1)*T+1+percon:find(list.sp=='F')*T)   = Ftarget'./(1+exp(y((find(list.sp=='F')-1)*T+1+percon:find(list.sp=='F')*T)));
end

% variables
[xn,xf,xg,Ag, An, Af,...
Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
F, N, G, E, Y, C, H, A_lag, Emnet, A,muu,...
pn, pg, pf, pee, w, wsn, wsf, wsg, tauf, taul,...
taurese, tauresg, ...
Trans, SWF, PV,PVSWF, objF]...
    = aux_SP(x, list, params, T, init, indic);

% inequality constraints
c=[];

% equality constraints
ceq =[];
ceq(1:1*T)     = C - (Y-xn-xg-xf);
ceq(1*T+1:2*T) = An-An_lag.*(1+gammaa.*(sn./rhon).^etaa.*(A_lag./An_lag).^phii); 
ceq(2*T+1:3*T) = Af-Af_lag.*(1+gammaa.*(sff./rhof).^etaa.*(A_lag./Af_lag).^phii);
ceq(3*T+1:4*T) = Ag-Ag_lag.*(1+gammaa.*(sg./rhog).^etaa.*(A_lag./Ag_lag).^phii);
ceq(4*T+1:5*T) = H - (Ln + Lf+Lg); % labor market clearing
ceq(5*T+1:6*T) = S - (sn + sff+sg); % scientists market clearing

ceq = ceq';
end