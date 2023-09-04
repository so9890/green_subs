    %  save(sprintf('2309_results_opt_main_notaul0_target%d', indic.target), 'x')
    % gg=  load('0308_results_opt_noskill_notaul0_notarget', 'x');
    % x=gg.x;

    %% add further direct optimization periods => goal: so that continuation value does not impact allocation anymore
    hhh= load(sprintf('2309_results_opt_main_notaul0_target%d', indic.target), 'x');
    x=hhh.x;
    T=12;
    count=0; % count number of iterations
    mm=10;
    Emsnew=Ems;
    Tinit=T;

    while mm>1e-5
        % save starting values
        Told=T;
        xold=x;
        %- transform and save for comparison
        out_trans=exp(x);
        if indic.noskill==0
            out_trans((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)=upbarH./(1+exp(x((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)));
            out_trans((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)=upbarH./(1+exp(x((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)));
        else
            out_trans((find(list.opt=='H')-1)*T+1:find(list.opt=='H')*T)=upbarH./(1+exp(x((find(list.opt=='H')-1)*T+1:find(list.opt=='H')*T)));
        end

        if indic.sep==1
            out_trans((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = upbarH./(1+exp(x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)));
            out_trans((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = upbarH./(1+exp(x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)));
            out_trans((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = upbarH./(1+exp(x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)));
        elseif indic.sep==0
            out_trans((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = (x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)).^2;
            out_trans((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = (x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)).^2;
            out_trans((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = (x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)).^2;
            out_trans((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T) = upbarS./(1+exp(x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T)));
        end    
        if indic.target==1
            out_trans((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)=Ftarget./(1+exp(x((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)));
        end
        out_transold=out_trans; 

        Emsnew=[Emsnew,0]; % add another net zero period to emissions limit
        Ftarget =  (Emsnew'+deltaa)/omegaa;

        % sequentially increase number of explicit periods
        count=count+1; 
        T=T+1;

        %- update initial values: add last period value as guess for new direct
        % optimization period
        x0 = zeros(nn*T,1);
            for ll=list.opt
                x0(T*(find(list.opt==ll)-1)+1:T*(find(list.opt==ll)))  = [xold(Told*(find(list.opt==ll)-1)+1:Told*(find(list.opt==ll)));xold(Told*(find(list.opt==ll)))]; 
            end
        %- optimization for new horizon
         constf=@(x)constraints_flexetaa(x, T, params, init201519, list, Emsnew, indic, MOM, percon, taulFixed);
         objf=@(x)objective(x, T, params, list, Ftarget, indic, init201519, percon, MOM, taulFixed);

        if ~isfile(sprintf('2309_results_opt_main_notaul0_target%d_Tplus%d.mat',indic.target, count))

             options = optimset('algorithm','sqp','TolCon',1e-9,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
             [x,fval,exitflag,output,lambda] = fmincon(objf,x0,[],[],[],[],lb,ub,constf,options);
             save(sprintf('2309_results_opt_main_notaul0_target%d_Tplus%d',indic.target, count), 'x')
        else
           helper=load(sprintf('2309_results_opt_main_notaul0_target%d_Tplus%d', indic.target, count));
          x=helper.x;
        end
          abbs=zeros(length(list.opt),1);
            
          %- transform
            out_trans=exp(x);
            if indic.noskill==0
                out_trans((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)=upbarH./(1+exp(x((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)));
                out_trans((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)=upbarH./(1+exp(x((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)));
            else
                out_trans((find(list.opt=='H')-1)*T+1:find(list.opt=='H')*T)=upbarH./(1+exp(x((find(list.opt=='H')-1)*T+1:find(list.opt=='H')*T)));
            end

            if indic.sep==1
                out_trans((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = upbarH./(1+exp(x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)));
                out_trans((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = upbarH./(1+exp(x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)));
                out_trans((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = upbarH./(1+exp(x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)));
            elseif indic.sep==0
                out_trans((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = (x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)).^2;
                out_trans((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = (x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)).^2;
                out_trans((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = (x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)).^2;
                out_trans((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T) = upbarS./(1+exp(x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T)));
            end    
        if indic.target==1
            out_trans((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)=Ftarget./(1+exp(x((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)));
        end


          for ll=1:length(list.opt)
              %- slice by variable
              old = out_transold(Told*(ll-1)+1:Told*(ll));
              new = out_trans(T*(ll-1)+1:T*(ll));
              %- compares only first 12 periods
                abbs(ll)=max(abs((new(1:Tinit)-old(1:Tinit))./old(1:Tinit))); % x= new, xold=old
          end
          mm=max(abbs);
          fprintf('max deviation in percent %d, number of finished iterations %d', mm*100, count)

    end