function [xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, H, A_lag, Emnet, A,muu,...
            pn, pg, pf, pee,  wsf, wsn, wsg,  tauf, taul, taus, taurese, tauresg, Trans,...
            w, SWF, PV,PVSWF, objF]= aux_OPT(x, list, params, T, init201519, indic, MOM, taul)

read_in_params;


% lambdaa = x((find(list.opt=='lambdaa')-1)*T+1:(find(list.opt=='lambdaa'))*T);
Lf     = x((find(list.opt=='Lf')-1)*T+1:find(list.opt=='Lf')*T);
Lg     = x((find(list.opt=='Lg')-1)*T+1:find(list.opt=='Lg')*T);
C      = x((find(list.opt=='C')-1)*T+1:find(list.opt=='C')*T);
F      = x((find(list.opt=='F')-1)*T+1:find(list.opt=='F')*T);
G      = x((find(list.opt=='G')-1)*T+1:find(list.opt=='G')*T);
H      = x((find(list.opt=='H')-1)*T+1:find(list.opt=='H')*T);
if indic.xgrowth==0
    sn      = x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T);
    sg      = x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T);
    sff      = x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T);
    
else
   sff=sff0*ones(size(F));
    sn=sn0*ones(size(F));
    sg=sg0*ones(size(F));
end

%% auxiliary variables
% loop over technology
    Af=zeros(T,1);
    Af_lag=[init201519(list.init=='Af0'); Af(1:T)]; % drop last value later
    Ag=zeros(T,1);
    Ag_lag=[init201519(list.init=='Ag0'); Ag(1:T)];
    An=zeros(T,1);
    An_lag=[init201519(list.init=='An0'); An(1:T)];
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

muu = C.^(-thetaa); % same equation in case thetaa == 1
% prices and taus
% machines green => comes first to get to taus

pf      = (F./(Af.*Lf)).^((1-alphaf)/alphaf)./(alphaf); % production fossil

if indic.notaul<2 % no earmarking
    taus = zeros(size(F));
    pg   = (G./(Ag.*Lg)).^((1-alphag)/alphag).*(1-taus)./alphag; % from production function green
    tauf = (G./F).^(1/eppse).*pg-pf; % optimality energy producers
    xg    = (alphag*pg).^(1/(1-alphag)).*Lg.*Ag;
else
    xg   = G.*(G./(Lg.*Ag)).^((1-alphag)/alphag); % from production green and demand green machines
    tauf = ((G./F).^(1/eppse).*xg./(G.*alphag)-pf)./(1+(G./F).^(1/eppse).*F./(G.*alphag)); % from demand final good producers and production green and gov. budget
    taus = tauf.*F./xg;
    pg   = (G./(Ag.*Lg)).^((1-alphag)/alphag).*(1-taus)./alphag; % from production function green
end


pee     = ((pf+tauf).^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse));
pn      = ((1-deltay.*pee.^(1-eppsy))./(1-deltay)).^(1/(1-eppsy)); % definition prices and numeraire

% output
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));
N       = (1-deltay)/deltay.*(pee./pn).^(eppsy).*E; % demand N final good producers 
Y       = (deltay^(1/eppsy)*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % final output production

Ln      = N./(An.*(pn.*alphan).^(alphan./(1-alphan))); % production neutral

% wages and policy elements
w       = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*(pf).^(1/(1-alphaf)).*Af; % labour demand fossil
    
%- wages scientists  
wsn         = (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N*(1-alphan).*An_lag)./(An.*rhon^etaa); 
%ws          =  wsn; % received by scientists => equal due to free movement of scientists
wsf_tilde   = (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*F*(1-alphaf).*Af_lag)./(Af.*rhof^etaa); 
wsg_tilde   = (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag).*Ag_lag)./(Ag.*rhog^etaa)./(1-taus);  % to include taus

if indic.subsres == 1 % research subsidies financed lump-sum are available 
    taurese   = 1-wsf_tilde./wsn;
    tauresg   = 1-(wsg_tilde./wsf_tilde);
else 
    taurese = zeros(size(wsn));
    tauresg = zeros(size(wsn));
end

% wage paid to scientists
wsf       = wsf_tilde./(1-taurese); 
wsg       = wsg_tilde./((1-taurese).*(1-tauresg));

xn      = (alphan*pn).^(1/(1-alphan)).*Ln.*An;
xf      = (alphaf*pf).^(1/(1-alphaf)).*Lf.*Af;

% in case taul is chosen optimally
if indic.notaul~=1 && indic.notaul~=3 % i.e. taul is not fixed
    taul = 1-(chii.*H.^(sigmaa))./(muu.*w);
end

if indic.notaul<2
    Trans    = ((taul).*w.*H)+tauf.*F;      % gov budget  
else
    Trans    = ((taul).*w.*H);      % gov budget; carbon tax revenues distributed as subsidies
end
Emnet     = omegaa*F-deltaa; % net emissions
A  = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);

% gammac =(1+gammaa.*(sff(T)./rhof).^etaa.*(A(T)./Af(T)).^phii)-1;

% utility
if thetaa~=1
    Utilcon = (C.^(1-thetaa))./(1-thetaa);
elseif thetaa==1
    Utilcon = log(C);
end
  

Utillab = chii.*(H.^(1+sigmaa))./(1+sigmaa);
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
    PVwork =1/(1-betaa)*(Utillab(T)); % this decreases last period work and science 
    PV= betaa^T*(PVconsump-PVwork);

    %Objective function value:
    %!! Dot product!!! so no dot.*
    % f = (-1)*(vec_discount*(Utilcon-Utillab- Utilsci)+PVcontUtil);+
    objF=(vec_discount*(SWF)+indic.PV*PV);

end