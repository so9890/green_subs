% if indic.GOV==1
%     GovRev=params(list.params=='GovRev');
% else
GovRev=0;
% end
S=params(list.params=='S');
if indic.sigmaWorker==0
    sigmaa = params(list.params=='sigmaa');
elseif indic.sigmaWorker==1 % more elastic labor supply
    sigmaa = 1/1.5;
elseif indic.sigmaWorker==2
    sigmaa = 1/0.5;
elseif indic.sigmaWorker==3 % inelasitc labor supply
    sigmaa = 1/0.0001;
end

sigmaas = params(list.params=='sigmaas');

if indic.util==0
    thetaa = params(list.params=='thetaa');
else
    if indic.Bop==1 % income effect dominates! => labor supply less responsive
        thetaa=(0.2+1/sigmaa)./((1-0.2)/sigmaa);
    else
        thetaa=0.4;
    end
end

upbarH  = params(list.params=='upbarH');
chii  = params(list.params=='chii');
if sum(list.params=='betaa')==1
    betaa = params(list.params=='betaa');
end

%if indic.labshareequ==0
    alphag = params(list.params=='alphag');
    alphaf = params(list.params=='alphaf');
    alphan = params(list.params=='alphan');
% else
%     alphag = 1/3*(params(list.params=='alphag')+params(list.params=='alphaf')+params(list.params=='alphan'));
%     alphaf = alphag;
%     alphan = alphag; 
% end
%if indic.subs==0
    eppsy = params(list.params=='eppsy');
% else
%     eppsy=1.3;
% end

if indic.elasE==0
    eppse = params(list.params=='eppse');

else 
    eppse =10;
end

deltay = params(list.params=='deltay');
gammaa = params(list.params=='gammaa');
etaa = params(list.params=='etaa');
rhon = params(list.params=='rhon');

%if indic.sizeequ ==0
    rhof = params(list.params=='rhof');
    rhog = params(list.params=='rhog');
% else
%     rhog=rhon;
%     rhof=rhon;
% end
if indic.know_spill==1
    phii=0;
elseif indic.know_spill ==0
    phii = params(list.params=='phii');
elseif indic.know_spill==2
    phii =0.25;
elseif indic.know_spill==3
    phii =0.75; % 
end

omegaa = params(list.params=='omegaa'); % carbon content of fossil energy
deltaa = params(list.params=='deltaa'); % natural sink

% utility externality from emissions
extexpp=1.02; 
weightext=0.01; % high weight: 

% growth rates in case of exogenous growth
vn=0;
vg=0;
vf=0;

% for exogenous growth fix research inputs at initial ss value
% if indic.xgrowth==1
%     error('update level of researchers for exogenous growth version')
% % hhelper= load(sprintf('params_0209_sep%d', indic.sep));
% if indic.sep==0 || indic.sep==2
%     sg0 =1.0305e-06; %hhelper.x0LF(hhelper.list.choice=='sg');% 
%     sff0=5.5660e-08; %hhelper.x0LF(hhelper.list.choice=='sff');
%     sn0= 0.3364; %hhelper.x0LF(hhelper.list.choice=='sn'); %
% elseif indic.sep==1
%     sg0 =0.1220;
%     sn0= 0.7922;
%     sff0= 0.0949; 
% else
%     error('initial values for scientists not given')
% end
% end