function [c,ceq]=laissez_faire_VECT_fmincon(x, params, list, varrs, laggs,T, indic, Ems)
c=[];

ceq=laissez_faire_VECT(x, params, list, varrs, laggs,T, indic, Ems);
end