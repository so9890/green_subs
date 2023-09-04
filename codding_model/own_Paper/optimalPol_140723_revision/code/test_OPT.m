function [f]=test_OPT(T, list,  params,symms, init201519, helper, indic, Ems)

read_in_params;

varrs=helper.LF_SIM; 
y=log(varrs);
z=sqrt(varrs);

% read in results

H= log((params(list.params=='upbarH')-varrs(list.allvars=='H', :))./(varrs(list.allvars=='H', :)))';
w=y(list.allvars=='w', :)';
Lf=y(list.allvars=='Lf', :)';
Ln=y(list.allvars=='Ln', :)';
Lg=y(list.allvars=='Lg', :)';
C=y(list.allvars=='C', :)';
F=y(list.allvars=='F', :)';
G=y(list.allvars=='G', :)';

Af = y(list.allvars=='Af', :)';
Ag = y(list.allvars=='Ag', :)';
An =y(list.allvars=='An', :)';

sff =z(list.allvars=='sff', :)';
sg =z(list.allvars=='sg', :)';
sn =z(list.allvars=='sn', :)';
ws=z(list.allvars=='ws', :)';
gammal =z(list.allvars=='gammal', :)';
pg=y(list.allvars=='pg', :)';
pn=y(list.allvars=='pn', :)';
pee=y(list.allvars=='pee', :)';
pf=y(list.allvars=='pf', :)';

Trans=y(list.allvars=='Trans', :)';

if indic.limit_LF==1
    syms tauf real
    symms.choice=[symms.choice, tauf];
    tauf=varrs(list.allvars=='tauf',:)';    
end
list.choice=string(symms.choice);

x0=eval(symms.choice);
x0=x0(:);

f=laissez_faire_VECT(x0, params, list, varrs, init201519,T, indic, Ems);
% 
%  if max(abs(f))>1e-9
% %    % equations where results are off
    % ind=1:length(f);
    %   pos=ind(abs(f)>1e-9);
    %   eqset= floor(pos./T);
% %     error('optimal policy does not solve laissez faire. For the relevant equations see eqset');
% %     % to evaluate stuff
%  end

 if max(abs(f))>1e-8
     error('LF function does not solve at 1e-8')
 else
     fprintf('Solution solves LF problem')
 end

end
