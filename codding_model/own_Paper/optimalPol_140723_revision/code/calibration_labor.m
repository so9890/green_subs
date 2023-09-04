function f = calibLabour(x,  MOM, C, Lnw, Lgw, Lfw, pf, F, paramss, list, poll)

% parameters
upbarH=paramss(list.paramsdir=='upbarH');
thetaa=paramss(list.paramsdir=='thetaa');
sigmaa=paramss(list.paramsdir=='sigmaa');

taul=poll(list.poldir=='taul');
tauf=poll(list.poldir=='tauf'); 

% variables

gammal = x(list.calib=='gammal')^2;
w = exp(x(list.calib=='w'));
Trans = exp(x(list.calib=='Trans'));
chii = exp(x(list.calib=='chii'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


% auxiliary equations
muu = C^(-thetaa);
H = MOM.targethour;
% equations
q=0;
        
% Government := > T (Transfers) 

q=q+1;
f(q) = -Trans- MOM.Debt + w.*H.*taul +tauf.*F;
         
% % budget => C
% q=q+1;
% f(q) = -C +(1-taul).* w.*H+T;

% skill market clearing
%9)
q=q+1;
f(q) = w.*H  - (Lnw+Lgw+Lfw); % 

%- supply: 
q=q+1;
f(q)= chii.*H.^sigmaa- (muu*(1-taul)*w-gammal); %=> determines chii

%12
q=q+1;
f(q)= gammal*(upbarH-H);

%fprintf('number equations: %d; number variables %d', q, length(list.calib));

end
