function f = easy_opt_taus(x, params, list,  init201519, indic)
% solves for optimal allocation by choosing s and h, where s=Lf/h 

  read_in_params;
    % exogenous variables: tauf, taul
    % endogenous variables: pg, pf, w, h, Lf, Lg, C, lambdaa, G, F, Y

    Af=(1+vf)*init201519(list.init=='Af0');
    Ag=(1+vg)*init201519(list.init=='Ag0');

    
if indic.notaul==0
    s = exp(x(1));
    h = exp(x(2));
    
    % auxiliary

  
    % labour market clearing: 
    Lg = (1-s)*h; 
    Lf = s*h;
    % production
    G = Ag*Lg;
    F = Af*Lf;
    Y=(F)^(eppsy)*(G)^(1-eppsy);
        deltaYF=eppsy*(G./F) ^(1-eppsy);
        deltaYG=(1-eppsy)*(F./G)^(eppsy);
    C=Y;
    pg= (1-eppsy)*(F./G)^eppsy;

    if indic.taxsch==1 % with subsidy
        pf= Lf./(Lg+Lf)*(G/F+deltaYF/deltaYG).*pg;% from free movement of labor, gov budget with taus, and optimality Y
    else
        pf= pg.*Ag./Af;
    end
        tauf =deltaYF-pf;
    if indic.taxsch==0
        taus = 0;
    else
        taus = tauf.*F./Lg;
    end
    w= pf.*Af; 
    taul = 1-h^sigmaa.*chii.*C^(thetaa)./w; % from hh optimality

    
    % derivatives
    dFdh = Af*s;
    dFds = Af*h;
    dCdh = (Af*s)^(eppsy)*(Ag*(1-s))^(1-eppsy); % from aggregate output
    dCds = dCdh*h*(eppsy/s-(1-eppsy)/(1-s));
     
    Uc = C^(-thetaa);
    Uh = -chii*h^sigmaa;
    Uf = -weightext*extexpp*(omegaa)^extexpp*F.^(extexpp-1);
    
    % Optimality conditions planner 
    f(1) = Uc*dCdh + Uh +indic.extern* Uf*dFdh; % optimality wrt h => taul
    f(2) = Uc*dCds + indic.extern*Uf*dFds; % optimality government wrt s => tauf

elseif indic.notaul==1
    s = exp(x(1));

    FG = Af*s/((1-s)*Ag);
    pg = (1-eppsy).*(FG)^eppsy; % optimal demand green
   
        deltaYF=eppsy*(1./FG) ^(1-eppsy);
        deltaYG=(1-eppsy)*(FG)^(eppsy);

    if indic.taxsch==1 % with taus
        pf = (1./FG+deltaYF./deltaYG).*pg.*s;  
    else
        pf = Ag.*pg./Af;
    end
    tauf = deltaYF-pf;
    if indic.taxsch==1
        taus = tauf.*Af.*s./(1-s);
        Th=0;
    else
        taus=0;
        Th=tauf.*Af.*s; % transfers per h
    end
    w = pf.*Af;
    h = (w/(chii*(w+Th)^thetaa))^(1/(sigmaa+thetaa));
    T = Th*h;
    C = w*h+T;
    Y = C;
    Lg = (1-s).*h;
    Lf = s.*h;
    G = Ag.*Lg;
    F = Af.*Lf;
    %%%% old

    % derivatives
    dFdh = Af*s;
    dFds = Af*h;
    dCdh = (Af*s)^(eppsy)*(Ag*(1-s))^(1-eppsy);
    dCds = dCdh*h*(eppsy*(1-s)-s*(1-eppsy))/(s*(1-s));
        
    Uc = C^(-thetaa);
%     Uh = -chii*h^sigmaa; Satisfied in auxiliary variables
    Uf = -weightext*extexpp*(omegaa)^extexpp*F.^(extexpp-1);
    
    % Optimality conditions planner 
    %f(1) = Uc*dCdh + Uh +indic.extern* Uf*dFdh; % optimality wrt h => taul
    f(1) = Uc*dCds + indic.extern*Uf*dFds; % optimality government wrt s => tauf

end
end