function [xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, H, A_lag, Emnet, A,muu,...
            pn, pg, pf, pee, w, wsn, wsf, wsg, tauf, taul,...
            taurese, tauresg, ...
            Trans, SWF, PV,PVSWF, objF]= aux_SP(x, list, params, T, init, indic)

read_in_params;


Lg    = x((find(list.sp=='Lg')-1)*T+1:find(list.sp=='Lg')*T);
Ln    = x((find(list.sp=='Ln')-1)*T+1:find(list.sp=='Ln')*T);
%Lf     = x((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T);
H     = x((find(list.sp=='H')-1)*T+1:find(list.sp=='H')*T);

xn     = x((find(list.sp=='xn')-1)*T+1:find(list.sp=='xn')*T);
xf     = x((find(list.sp=='xf')-1)*T+1:find(list.sp=='xf')*T);
xg     = x((find(list.sp=='xg')-1)*T+1:find(list.sp=='xg')*T);

Af     = x((find(list.sp=='Af')-1)*T+1:find(list.sp=='Af')*T);
Ag     = x((find(list.sp=='Ag')-1)*T+1:find(list.sp=='Ag')*T);
An     = x((find(list.sp=='An')-1)*T+1:find(list.sp=='An')*T);
sff     = x((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T);
sg      = x((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T);
sn      = x((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T);

C      = x((find(list.sp=='C')-1)*T+1:find(list.sp=='C')*T);
F      = x((find(list.sp=='F')-1)*T+1:find(list.sp=='F')*T);

if indic.xgrowth==0
    % initial values
    An0=init(list.init=='An0');
    Ag0=init(list.init=='Ag0');
    Af0=init(list.init=='Af0');

    % aux variables

    %A       = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);
    Af_lag  = [Af0;Af(1:T-1)]; % shift Af backwards
    Ag_lag  = [Ag0;Ag(1:T-1)];
    An_lag  = [An0;An(1:T-1)];
else
  if indic.zero==0
    Af=zeros(T,1);
    Af_lag=[init(list.init=='Af0'); Af(1:T)]; % drop last value later
    Ag=zeros(T,1);
    Ag_lag=[init(list.init=='Ag0'); Ag(1:T)];
    An=zeros(T,1);
    An_lag=[init(list.init=='An0'); An(1:T)];
    A_lag=zeros(T,1);

        for i=1:T
        A_lag(i)   = (rhof*Af_lag(i)+rhon*An_lag(i)+rhog*Ag_lag(i))./(rhof+rhon+rhog);

        Af(i)=Af_lag(i).*(1+gammaa*(sff(i)/rhof).^etaa.*(A_lag(i)/Af_lag(i))^phii);
        Ag(i)=Ag_lag(i).*(1+gammaa*(sg(i)/rhog).^etaa.*(A_lag(i)/Ag_lag(i))^phii);
        An(i)=An_lag(i).*(1+gammaa*(sn(i)/rhon).^etaa.*(A_lag(i)/An_lag(i))^phii);

        %-update lags

        Af_lag(i+1)=Af(i);
        Ag_lag(i+1)=Ag(i);
        An_lag(i+1)=An(i);

        end

        Af_lag=Af_lag(1:end-1);
        An_lag=An_lag(1:end-1);
        Ag_lag=Ag_lag(1:end-1);
    else % version with zero growth
        An_lag=init(list.init=='An0');
        Ag_lag=init(list.init=='Ag0');
        Af_lag=init(list.init=='Af0');

        Af=zeros(size(F));
        Ag=zeros(size(F));
        An=zeros(size(F));

        for i=1:T
            An(i)=(1+vn)*An_lag;
            Ag(i)=(1+vg)*Ag_lag;
            Af(i)=(1+vf)*Af_lag;
            %- update laggs
            An_lag=An(i);
            Af_lag=Af(i);
            Ag_lag=Ag(i);
        end

        An_lag=An;
        Af_lag=Af;
        Ag_lag=Ag;
    end
end

G       = xg.^alphag.*(Ag.*Lg).^(1-alphag); 
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));
Lf      = (F./xf.^alphaf).^(1/(1-alphaf))./Af;
pg      = (G./(Ag.*Lg)).^((1-alphag)/alphag)./alphag;
pf      = (F./(Af.*Lf)).^((1-alphaf)/alphaf)./(alphaf);
tauf    = (G./F).^(1/eppse).*pg-pf; 
pee     = (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse));


N       = xn.^alphan.*(An.*Ln).^(1-alphan); 
pn      = (N./(An.*Ln)).^((1-alphan)/alphan)./alphan;
A_lag   = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)./(rhof+rhon+rhog);
A       = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);
Y       = (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % final output production
 
% prices compatible with sp solution 
% scientists wages paid by mashine producers
wsn       = (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N*(1-alphan).*An_lag)./(An.*rhon^etaa); 
wsf_tilde = (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*F*(1-alphaf).*Af_lag)./(Af.*rhof^etaa); 
wsg_tilde = (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag).*Ag_lag)./(Ag.*rhog^etaa);  
% wages scienties should be different! as sp acknowledges dynamic
% efficiency)=> can be used to derive research subsidies
taurese   = 1-wsf_tilde./wsn;
tauresg   = 1-(wsg_tilde./wsf_tilde);
% wage paid to scientists
wsf       = wsf_tilde./(1-taurese); 
wsg       = wsg_tilde./((1-taurese).*(1-tauresg));

% wages labor
w      = pg.*(1-alphag).*G./Lg; % wages labor should be equalized in efficient allocation 
muu = C.^(-thetaa);


taul = 1-(chii.*H.^(sigmaa))./(muu.*w);
Trans    = ((1-taul).*w.*H)+tauf.*F;
Emnet     = omegaa*F-deltaa; % net emissions


% growth rate consumption

%gammac = gammaa.*(sff(T)./rhof).^etaa.*(A(T)./Af(T)).^phii;

% utility
if thetaa~=1
    Utilcon = (C.^(1-thetaa))./(1-thetaa);
elseif thetaa==1
    Utilcon = log(C);
end

Utillab = chii*(H.^(1+sigmaa))./(1+sigmaa);
SWF = Utilcon-Utillab;
 %- create discount vector
     disc=repmat(betaa, 1,T);
     expp=0:T-1;
     vec_discount= disc.^expp;
     PVSWF = vec_discount*SWF;
    
    % continuation value
    %- last period growth rate as proxy for future growth rates
    gammay = Y(T)/Y(T-1)-1;
    PVconsump= 1/(1-betaa*(1+gammay)^(1-thetaa))*Utilcon(T);
    PVwork = 1/(1-betaa)*(Utillab(T)); % this decreases last period work and science 
    PV= betaa^T*(PVconsump-PVwork);

    %Objective function value:
    %!! Dot product!!! so no dot.*
    % f = (-1)*(vec_discount*(Utilcon-Utillab- Utilsci)+PVcontUtil);+
    objF=(vec_discount*(SWF)+indic.PV*PV);


end