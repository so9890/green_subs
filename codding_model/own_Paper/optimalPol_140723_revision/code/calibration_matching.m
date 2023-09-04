function [x0LF, SL, SP, SR, Sall, Sinit201014, init201014 , Sinit201519, init201519, Sparams, Spol, params, pol, symms, MOM, indexx, list]...
    = calibration_matching(MOM, symms, list, parsHelp, polhelp, indic)
% function to match moments to model equations
% i.e. solving model plus additional equations for paramters
% to this end, the model is solved 
% 1) from an aggregate pproduction side
% 2) the household probelm and labour producing firms
% 3) research side. 

% To be chosen: 
% An0, Af0, Ag0 (2019)
% deltay: weight on energy in production function 

% output
% indexx: map of indices indicating variable transformation for code.
%         contains laissez-faire index and calibration index. 

if polhelp(list.poldir=='tauf')~=0
error('tauf not =0; have to change calibration and pay attention how env revs are recycled!')
end
%% 1) initaliase stuff
% a) define symbolic model variables and parameters which have to be calibrated
% b) define indexx for model variables for transformation

%- calibration research side
syms Af_lag Ag_lag An_lag sff sg sn  ws gammaa real
symms.calibRes = [Af_lag Ag_lag An_lag sff sg sn ws gammaa];
 
list.calibRes =string(symms.calibRes);

%- all variables: to save base year variables!
syms muu H C F G N Y E Lf Ln Lg Af Ag An sff sg sn ...
     w pg pn pee pf gammal xn xg xf Emnet A ...
     tauf taul taus taurese tauresg Trans Ln Lg Lf SWF real
symms.allvars= [muu, H, C, F, G, N, Y, E, Af, Ag, An, sff, sg, sn, ...
                w,ws, pg, pn, pee, pf, gammal, xn, xg, xf, Emnet, A, ...
                tauf, taul, taus, taurese, tauresg, Trans, Ln, Lg, Lf, SWF];
list.allvars  = string(symms.allvars);

%- variables and index for laissez faire
symms.choice = [H, C, F, G, Lf, Ln, Lg, Af, Ag, An, sff, sg, sn,...
                w, pg, pn, pee, pf, gammal, Trans, ws];
list.choice  = string(symms.choice);

%- calibration productivity
syms deltay omegaa real
symms.prod= [pn, pg, deltay omegaa];
list.prod = string(symms.prod);

%- calibration labour firms and HH side
syms chii real
symms.calib = [Trans gammal w chii];
list.calib = string(symms.calib);

%-b) generate index to transform model variables for quantitative solution
%-- index LF solution
indexxLF.lab = boolean(zeros(size(list.choice)));
indexxLF.exp = boolean(zeros(size(list.choice)));
indexxLF.sqr = boolean(zeros(size(list.choice)));

indexxLF.lab(list.choice=='H')=1;
indexxLF.exp(list.choice~='gammal'&list.choice~='H' )=1;
indexxLF.sqr(list.choice=='gammal')=1;

%-- index for calibration
indexxcalib.lab = boolean(zeros(size(list.calib)));
indexxcalib.exp = boolean(zeros(size(list.calib)));
indexxcalib.sqr = boolean(zeros(size(list.calib)));

indexxcalib.exp( list.calib~='gammal' )=1;
indexxcalib.sqr( list.calib=='gammal')=1;

%- save indices to map 
indexx = containers.Map({'LF', 'calib'}, {indexxLF, indexxcalib});


%% 2) Calibration
% a) production side
% b) household and labor supply side
% c) research side
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %-a) Production side: deltay

%- initial values
pn=log(1);
pg=log(2);
deltay =  log((1-0.4)/0.4);
x0=eval(symms.prod(list.prod~='omegaa'));

%- solve
f = calibration_production(x0, MOM, list, parsHelp, polhelp);
prodf = @(x)calibration_production(x,  MOM, list, parsHelp, polhelp);
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[solProd, fval] = fsolve(prodf, x0, options);

%- transform variables
trProd=exp(solProd);
trProd(list.prod=='deltay')=1/(1+trProd(list.prod=='deltay'));

% required for next functions
[C, Lnw, Lgw, Lfw, pf, F, pn, pg, pee, E, Y, N, G, xn, xg, xf, ...
 AfLf, AgLg, AnLn, omegaa, deltay]=results_production(list, trProd, MOM, parsHelp, polhelp, 'calib');
trProd =eval(symms.prod);
%% Labour side

%- initial guess
H = MOM.targethour;
w = 8;
gammal = 0;
chii =10;
Trans  = w*H*polhelp(list.poldir=='taul')+polhelp(list.poldir=='tauf')*F+MOM.Debt;

x0=eval(symms.calib);

%- transform
guess_trans=trans_guess(indexxcalib, x0, parsHelp, list.paramsdir);

%test:
f =calibration_labor(guess_trans,  MOM, C, Lnw, Lgw, Lfw, pf, F, parsHelp, list, polhelp);

%- solve
Labf = @(x)calibration_labor(x,  MOM, C, Lnw, Lgw, Lfw, pf, F, parsHelp, list, polhelp);
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[solLab, fval] = fsolve(Labf, guess_trans, options);
trLab=trans_allo_out(indexxcalib, solLab, parsHelp, list.paramsdir, indic);

%% Research side
%- from previous
% cell_par=arrayfun(@char, symms.calib, 'uniform', 0);
% Sparams=cell2struct(num2cell(trLab), cell_par, 2);
w=trLab(list.calib=='w');

%- initial guess
 sg  = log(0.001);
 ws  = log(3);
 gammaa =log(3);
 sff = log(0.001);
 sn  = log(parsHelp(list.paramsdir=='S')-exp(sg)-exp(sff)); 
    Af = AfLf/Lfw*w;
    An = AnLn/Lnw*w;
    Ag = AgLg/Lgw*w;
 Af_lag = log(Af/(1.02));
 An_lag = log(An/(1.02));
 Ag_lag = log(Ag/(1.02));

x0 = eval(symms.calibRes);

%- test
% f= calibration_research(x0, MOM, list, trProd, parsHelp, polhelp, Af, An, Ag);
% 
% % solving model
modF3 = @(x)calibration_research(x, MOM, list, trProd, parsHelp, polhelp, Af, An, Ag); 
options = optimoptions('fsolve', 'TolFun', 10e-10, 'MaxFunEvals',8e4, 'MaxIter', 3e5); %, 'Display', 'Iter','Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[x, fval, exitf] = fsolve(modF3, x0, options);
% save(sprintf('calib_0501_sep'), 'xSc');
trRes=exp(x);

%% save results
 
[Sparams, Spol, params, pol,...
    x0LF, SL, SP, SR, Sall, ...
    Sinit201014, init201014 , Sinit201519, init201519]...
    =calibration_results(symms,trProd, trLab, trRes, parsHelp, list, polhelp, MOM, indic);

%% Test if is calibration and baseline model solve LF in baseyear
guess_transLF=trans_guess(indexxLF, x0LF, params, list.params);
f=laissez_faire(guess_transLF, params, list, pol, init201014, indic);

if max(abs(f))>1e-9
    error('calibration is not a solution to LF')
else
    fprintf('Hurray!!! LF solves at baseline calibration!!!');
end
end