function [c, ceq] = constraints(y, T, params, init, list, Ems, indic, MOM, percon, taulFixed)
% function to read in constraints on government problem

% pars
read_in_params;

% Ftarg_20s=(MOM.US_Budget20_30+3*deltaa)./omegaa; 

% transform x: all are exponentially transformed
 x=exp(y);

% except for hours
x((find(list.opt=='H')-1)*T+1:find(list.opt=='H')*T) = upbarH./(1+exp(y((find(list.opt=='H')-1)*T+1:find(list.opt=='H')*T)));

if indic.target==1
    Ftarget=(Ems'+deltaa)./omegaa; 
    x((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)   = Ftarget./(1+exp(y((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)));
end


%- auxiliary variables

[xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, H, A_lag, Emnet, A,muu,...
            pn, pg, pf, pee,  wsf, wsn, wsg,  tauf, taul,taus, taurese, tauresg, Trans,...
            w, SWF, PV,PVSWF, objF]= aux_OPT(x, list, params, T, init, indic, MOM, taulFixed);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%    Inequality Constraints    %%%
 % only for direct periods
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % time period specific constraints! 
c = []; %  periods and 2 additional ones: 
c(1:T) = -Trans; % lump- sum transfers have to be positive: IMPORTANT NOTE 24/07: gov would subsidize fossil sector but minimally!
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%    Equality Constraints    %%%
 % include missing equations here %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ceq = [];

     if indic.xgrowth==0
            ceq(T*0+1:T*1) = Ln -pn.*(1-alphan).*N./w; % labour demand by sector 
            ceq(T*1+1:T*2) = Lg - pg.*(1-alphag).*G./w;
            ceq(T*2+1:T*3) = Lf- H./(1+Ln./Lf+Lg./Lf); % labour market clearing 
            ceq(T*3+1:T*4) = -C+(1-taul).*w.*H+Trans;
            ceq(T*4+1:T*5) = wsg-wsf; 
            ceq(T*5+1:T*6) = wsg-wsn; % free movement of scientists
            ceq(T*6+1:T*7) = S-sn-sg-sff; % scientists market clearing
            if indic.notaul==1 || indic.notaul==3 % then taul is fixed
                ceq(T*7+1:T*8) = chii*H.^sigmaa-(muu.*(1-taul).*w);
            end

     elseif indic.xgrowth==1
            ceq(T*0+1:T*1) = Ln -pn.*(1-alphan).*N./w; % labour demand by sector 
            ceq(T*1+1:T*2) = Lg -pg.*(1-alphag).*G./w;
            ceq(T*2+1:T*3) = Lf- H./(1+Ln./Lf+Lg./Lf); % labour market clearing 
            ceq(T*3+1:T*4) = -C+(1-taul).*w.*H+Trans;

            if indic.notaul==1 || indic.notaul==3 
                ceq(T*4+1:T*5)=  chii*H.^sigmaa-(muu.*(1-taul).*w);
            end

    end
 

 %
ceq = ceq';
end