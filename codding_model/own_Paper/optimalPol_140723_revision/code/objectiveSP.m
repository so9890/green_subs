function f = objectiveSP(y, T, params, list, Ftarget, indic, init, percon)

% read in stuff

% pars
read_in_params;

% transform x: all are exponentially transformed
 x=exp(y);
% % except for hours
x((find(list.sp=='H')-1)*T+1:find(list.sp=='H')*T) = upbarH./(1+exp(y((find(list.sp=='H')-1)*T+1:find(list.sp=='H')*T)));

% x((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T) = (y((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T)).^2;
% x((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T) = (y((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T)).^2;
% x((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T) = (y((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T)).^2;
 
if indic.target==1
    x((find(list.sp=='F')-1)*T+1+percon:find(list.sp=='F')*T) = Ftarget'./(1+exp(y((find(list.sp=='F')-1)*T+1+percon:find(list.sp=='F')*T)));
end


%- auxiliary

[xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, H, A_lag, Emnet, A,muu,...
            pn, pg, pf, pee, w, wsn, wsf, wsg, tauf, taul,...
            taurese, tauresg, ...
            Trans, SWF, PV,PVSWF, objF]= aux_SP(x, list, params, T, init, indic);


f = (-1)*objF;


end





