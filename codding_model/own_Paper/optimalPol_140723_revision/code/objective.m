function f = objective(y, T, params, list, Ftarget, indic, init201519, percon, MOM, taulFixed)

% pars
read_in_params;

%-- transform variables

x=exp(y);
% hours
x((find(list.opt=='H')-1)*T+1:find(list.opt=='H')*T) = upbarH./(1+exp(y((find(list.opt=='H')-1)*T+1:find(list.opt=='H')*T)));
if indic.target==1
    x((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)   = Ftarget./(1+exp(y((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)));
end

%- auxiliary variables
[xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, H, A_lag, Emnet, A,muu,...
            pn, pg, pf, pee,  wsf, wsn, wsg,  tauf, taul,taus, taurese, tauresg, Trans,...
            w, SWF, PV,PVSWF, objF]= aux_OPT(x, list, params, T, init201519, indic, MOM, taulFixed);

% %Objective function value:
% %!! Dot product!!! so no dot.*
%  f(1) = (-1)*(vec_discount*(Utilcon-Utillab- Utilsci-indic.extern*weightext*(omegaa.*F).^extexpp)+indic.PV*PV);
f=-1*objF;
end





