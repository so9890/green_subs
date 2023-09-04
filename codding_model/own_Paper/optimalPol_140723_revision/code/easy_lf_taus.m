% simple model to compare social planner and 
% competitive equilibrium



%- model for numeric solution as function of pg and Lg
function f = easy_lf_taus(x, params, list, pol,  init201519, indic)

    Lf = exp(x(1));
    Lg = exp(x(2));
    
    % auxiliary

    read_in_params;
    tauf    = pol(:,list.pol=='tauf');
    taul    = pol(:,list.pol=='taul');

    % exogenous variables: tauf, taul
    % endogenous variables: pg, pf, w, h, Lf, Lg, C, lambdaa, G, F, Y

    Af=(1+vf)*init201519(list.init=='Af0');
    Ag=(1+vg)*init201519(list.init=='Ag0');

    G = Ag*Lg;
    F = Af*Lf;  % production
    if indic.taxsch==1 % with subsidies
        taus = tauf*F./Lg; % gov environmental budget
    elseif indic.taxsch==0
        taus=0;
    end
    pg = (1-eppsy)*(F./G)^eppsy; % optimal demand green input
    w = pg*Ag+taus; % labor demand green
    pf = w/(Af); % wage clearing 
    %Y = (F)^(eppsy)*(G)^(1-eppsy);
    hdem = Lf+Lg;
   if indic.taxsch==1
       T=taul.*w.*hdem;
   else
       T=taul.*w.*hdem+ tauf.*F;
   end
     hsup = ((((1-taul).*w+T./hdem)^(-thetaa)*w*(1-taul))/(chii))^(1/(sigmaa+thetaa));
    
        
    % 
    f(1)= pf+tauf - eppsy.*(G/F)^(1-eppsy); %1-(pf/eppsy)^eppsy*(pg/(1-eppsy))^(1-eppsy); % determines pg; could have alternatively 
    f(2)= hdem-hsup; % labour market clearing determines Lg; (goods market clears by walras law)
end