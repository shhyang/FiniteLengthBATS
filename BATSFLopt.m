classdef BATSFLopt < BATSFiniteLength & handle
   
    properties
        DegreeDist_inact
        DegreeDist_OH
    end
    
    properties (Hidden)
        
        filename
        
    end
    
    
    methods
        function obj = BATSFLopt(M,K,q,h)
            obj = obj@BATSFiniteLength(M,K,q,h);
        end
        
        function [dd_out, dd_diff] = ddOpt(obj,n,ddin,numiter,flag)
           % optimization 
            if nargin == 3
                numiter = 100;
                flag = 'BP';
                warning('Run 100 iterations for BP decoding!');
            elseif nargin == 4
                flag = 'BP';
                warning('Use BP decoding!');
            end
            
            if ~strcmp(flag,'BP') && ~strcmp(flag,'inac')
                error('flag = ''BP'' or ''inac''!');
            end
            
            degreeDist = ddin;
            if length(ddin) > obj.InputNumber
                degreeDist=ddin(1:obj.InputNumber)/sum(ddin(1:obj.InputNumber));
            elseif length(ddin) < obj.InputNumber
                degreeDist(obj.InputNumber) = 0;
            end
            
            obj.filename = [date '_K' num2str(obj.InputNumber) '_M' num2str(obj.BatchSize) '_n' num2str(n)];
            
            switch flag
                case 'BP'
                    obj.filename = ['ddBP' obj.filename];
                    objfunc = @() obj.errorProb(obj.PoissonRec(n,'BP'));
                    disp('Start to optimize the degree distribution for BP decoding.');
                case 'inac'
                    obj.filename = ['ddInac' obj.filename];
                    objfunc = @() obj.expInac(obj.PoissonRec(n,'inac'));
                    disp('Start to optimize the degree distribution for inactivation decoding.');
            end
            
            tic;
                                    
            [dd_out, dd_diff] = obj.runDDOpt(objfunc,degreeDist,numiter);
            
            timeused = toc;
            
            disp(['Optimization runs for ' num2str(timeused) ' seconds:']);                        
        end
        
        function [ddout, ddiff] = runDDOpt(obj,objfunc,ddin,numiter)
            iter = 0;
            delta = 0.01;
            
            ddiff = zeros(size(ddin)); 
            
            obj.setDegreeDist(ddin);
            best_ei = objfunc();
                        
            while iter <= numiter
                disp(['Iteration ' num2str(iter)]);
                new_diff = ddiff;
                better_diff = 0;
                
                for d = 1:obj.InputNumber
                    % try to increase the prob mass of d
                    ii = 1;
                    df2 = ddiff;
                    while 1
                        df2(d) = ddiff(d) + ii*delta;
                        ddout = ddin + df2;
                        ddout = ddout/sum(ddout);
                        
                        obj.setDegreeDist(ddout);
                        new_ei = objfunc();
                        
                        if new_ei >= best_ei
                            break;
                        end
                        
                        best_ei = new_ei;
                        new_diff = df2;
                        ii = ii+1;
                    end
                    
                    if ii > 1
                        better_diff = 1;
                        disp(['Better degree dist found (' num2str(best_ei) ') + d = ' num2str(d)]);
                    end
                    
                    % try to decrease the prob mass of d
                    if ddin(d) + ddiff(d) == 0
                        continue;
                    end
                    ii = 1;
                    last_run = 0;
                    while last_run == 0
                        df2(d) = ddiff(d) - ii*delta;
                        if df2(d)+ddin(d)<0
                            last_run = 1;
                            df2(d) = - ddin(d);
                        end
                        ddout = ddin + df2;
                        ddout = ddout/sum(ddout);
                        
                        obj.setDegreeDist(ddout);
                        new_ei = objfunc();
                        
                        if new_ei >= best_ei
                            break;
                        end
                        
                        best_ei = new_ei;
                        new_diff = df2;
                        ii = ii+1;
                    end
                    
                    if ii > 1
                        better_diff = 1;
                        disp(['Better degree dist found (' num2str(best_ei) ')- d = ' num2str(d)]);
                    end
                end
                
                if better_diff == 1
                    ddiff = new_diff;
                    save(obj.filename,'ddiff','-ASCII','-append');
                else
                    disp(['No better degree dist find at iteration: ' num2str(iter)]);
                    break;
                end
                iter = iter + 1;              
            end
            
            ddout = ddin + ddiff;
            ddout = ddout/sum(ddout);
        end

        
        function degout = minExpInact(obj,n,degin,numiter,flag)
            % optimization using poission reception
            if nargin == 3
                numiter = 100;
                warning('Run 100 iterations!');
                flag = 'poisson';
            elseif nargin == 4
                flag = 'poisson';
            end
            
            degreeDist = degin;
            if length(degin) > obj.InputNumber
                degreeDist=degin(1:obj.InputNumber)/sum(degin(1:obj.InputNumber));
            elseif length(degin) < obj.InputNumber
                degreeDist(obj.InputNumber) = 0;
            end
            
            obj.filename = ['ddInac' date '_K' num2str(obj.InputNumber) '_M' num2str(obj.BatchSize) '_n' num2str(n)];
            
            switch flag
                case 'poisson'
                    obj.filename = ['P' obj.filename];
                    objfunc = @() sum(obj.PoissonRec(n,'inac'));
                case 'exact'
                    obj.filename = ['E' obj.filename];
                    objfunc = @() sum(obj.BPDecoding(n,'inac'));
                otherwise
                    error('flag = ''poisson'' or ''exact'' only');
            end
            
            tic;
                                    
            [obj.DegreeDist_inact, numdeg, iter] = obj.runOpt_Inac(objfunc,degreeDist,numiter);
            
            timeused = toc;
            
            disp(['Optimization runs for ' num2str(timeused) ' seconds:']);
            disp([num2str(iter) ' iterations']);
            disp([num2str(numdeg) ' degree values are modified']);
            
            degout = obj.DegreeDist_inact;
        end  
        
        function [ddout, numdeg, iter] = runOpt_Inac(obj,objfunc,ddin,numiter)
            
            delta = 0.01;
            
            iter = 0;
            
            % try to find a better degree distribution by disturbing degin
            % a little from psi_1 to psi_K
            
            ddout = ddin;         
            newDD = 0;        
            numdeg = 0;
            d = 1;           
            sign = 1;
            
            obj.setDegreeDist(ddin);
            best_ei = objfunc();
            best_d = ddout;
            newD = 0;
            
            while iter <= numiter
                                
                dd2 = obj.disturb(best_d,d,sign*delta);
                
                obj.setDegreeDist(dd2);
                ei = objfunc();
                
                if ei < best_ei 
                    % a better degree dist is found                
                    best_d = dd2;
                    best_ei = ei;
                    if newD == 0
                        newD = 1;
                        disp(['Better degree dist found at d = ' num2str(d)]);
                    end
                else % ei <= expInac
                    if newD == 1
                        best_dd = best_d;
                        best_d = ddout;
                        newDD = 1;
                        newD = 0;
                    end
                    if sign > 0 && best_d(d)>delta
                        sign = -1;
                    elseif d < obj.InputNumber
                        d = d + 1;
                        sign = 1;
                    else % d == K
                        if newDD == 1
                            numdeg = numdeg + 1;
                            d = 1;
                            ddout = best_dd;
                            best_d = ddout;
                            newDD = 0;
                            save(obj.filename,'ddout','-ASCII','-append');
                            disp('New degree dist saved');
                        else
                            disp('No better degree dist find');
                            break;
                        end
                    end
                end
                iter = iter + 1;
            end     
            
            if iter > numiter
                iter = numiter;
                if newDD == 1
                    ddout = best_dd;
                    numdeg = numdeg + 1;
                    save(obj.filename,'ddout','-ASCII','-append');
                    disp('New degree dist saved');
                end
            end
        end
        
       
        function [dd_out, dd_diff] = minErrProb(obj,n,degin,numiter,flag)
           % optimization using poission reception
            if nargin == 3
                numiter = 100;
                warning('Run 100 iterations with Poisson evaluation!');
                flag = 'poisson';
            elseif nargin == 4
                flag = 'poisson';
            end
            
            degreeDist = degin;
            if length(degin) > obj.InputNumber
                degreeDist=degin(1:obj.InputNumber)/sum(degin(1:obj.InputNumber));
            elseif length(degin) < obj.InputNumber
                degreeDist(obj.InputNumber) = 0;
            end
            
            obj.filename = ['ddEP' date '_K' num2str(obj.InputNumber) '_M' num2str(obj.BatchSize) '_n' num2str(n)];
            
            switch flag
                case 'poisson'
                    obj.filename = ['P' obj.filename];
                    objfunc = @() obj.errorProb(obj.PoissonRec(n,'BP'));
                case 'exact'
                    obj.filename = ['E' obj.filename];
                    objfunc = @() obj.errorProb(obj.BPDecoding(n,'BP'));
                otherwise
                    error('flag = ''poisson'' or ''exact'' only');
            end
            
            tic;
                                    
            [dd_out, dd_diff] = obj.runOpt_ErrProb2(objfunc,degreeDist,numiter);
            
            timeused = toc;
            
            disp(['Optimization runs for ' num2str(timeused) ' seconds:']);                        
        end
        
        function [ddout, ddiff] = runOpt_ErrProb(obj,objfunc,ddin,numiter)
            iter = 0;
            delta = 0.01;
            
            ddiff = zeros(size(ddin)); 
            
            obj.setDegreeDist(ddin);
            best_ei = objfunc();
                        
            while iter <= numiter
                disp(['Iteration ' num2str(iter)]);
                better_diff = 0;
                
                for d = 1:obj.InputNumber
                    % try to increase the prob mass of d
                    ii = 1;
                    df2 = ddiff;
                    while 1
                        df2(d) = ddiff(d) + delta;
                        ddout = ddin + df2;
                        ddout = ddout/sum(ddout);
                        
                        obj.setDegreeDist(ddout);
                        new_ei = objfunc();
                        
                        if new_ei >= best_ei
                            break;
                        end
                        
                        best_ei = new_ei;
                        ddiff(d) = df2(d);
                        ii = ii+1;
                    end
                    
                    if ii > 1
                        better_diff = 1;
                        disp(['Better degree dist found + d = ' num2str(d)]);
                    end
                    
                    % try to decrease the prob mass of d
                    if ddin(d) + ddiff(d) == 0
                        continue;
                    end
                    ii = 1;
                    last_run = 0;
                    while last_run == 0
                        df2(d) = ddiff(d) - delta;
                        if df2(d)+ddin(d)<0
                            last_run = 1;
                            df2(d) = - ddin(d);
                        end
                        ddout = ddin + df2;
                        ddout = ddout/sum(ddout);
                        
                        obj.setDegreeDist(ddout);
                        new_ei = objfunc();
                        
                        if new_ei >= best_ei
                            break;
                        end
                        
                        best_ei = new_ei;
                        ddiff(d) = df2(d);
                        ii = ii+1;
                    end
                    
                    if ii > 1
                        better_diff = 1;
                        disp(['Better degree dist found - d = ' num2str(d)]);
                    end
                end
                
                if better_diff == 0
                    disp(['No better degree dist find at iteration: ' num2str(iter)]);
                    break;
                end
                iter = iter + 1;              
            end
            
            ddout = ddin + ddiff;
            ddout = ddout/sum(ddout);
        end

        
        function [ddout, ddiff] = runOpt_ErrProb2(obj,objfunc,ddin,numiter)
            iter = 0;
            delta = 0.01;
            
            ddiff = zeros(size(ddin)); 
            
            obj.setDegreeDist(ddin);
            best_ei = objfunc();
                        
            while iter <= numiter
                disp(['Iteration ' num2str(iter)]);
                new_diff = ddiff;
                better_diff = 0;
                
                for d = 1:obj.InputNumber
                    % try to increase the prob mass of d
                    ii = 1;
                    df2 = ddiff;
                    while 1
                        df2(d) = ddiff(d) + ii*delta;
                        ddout = ddin + df2;
                        ddout = ddout/sum(ddout);
                        
                        obj.setDegreeDist(ddout);
                        new_ei = objfunc();
                        
                        if new_ei >= best_ei
                            break;
                        end
                        
                        best_ei = new_ei;
                        new_diff = df2;
                        ii = ii+1;
                    end
                    
                    if ii > 1
                        better_diff = 1;
                        disp(['Better degree dist found + d = ' num2str(d)]);
                    end
                    
                    % try to decrease the prob mass of d
                    if ddin(d) + ddiff(d) == 0
                        continue;
                    end
                    ii = 1;
                    last_run = 0;
                    while last_run == 0
                        df2(d) = ddiff(d) - ii*delta;
                        if df2(d)+ddin(d)<0
                            last_run = 1;
                            df2(d) = - ddin(d);
                        end
                        ddout = ddin + df2;
                        ddout = ddout/sum(ddout);
                        
                        obj.setDegreeDist(ddout);
                        new_ei = objfunc();
                        
                        if new_ei >= best_ei
                            break;
                        end
                        
                        best_ei = new_ei;
                        new_diff = df2;
                        ii = ii+1;
                    end
                    
                    if ii > 1
                        better_diff = 1;
                        disp(['Better degree dist found - d = ' num2str(d)]);
                    end
                end
                
                if better_diff == 1
                    ddiff = new_diff;
                else
                    disp(['No better degree dist find at iteration: ' num2str(iter)]);
                    break;
                end
                iter = iter + 1;              
            end
            
            ddout = ddin + ddiff;
            ddout = ddout/sum(ddout);
        end

        
        function [ddout,Delta]= maxErrExponent(obj,me)
            D = obj.InputNumber;
            K = obj.InputNumber;
            
            if nargin == 1
                A = [];
                b = [];
            elseif nargin == 2
                A = [0 1:D];
                b = me;
            end
                        
            
            f = [1; zeros(D,1)];
            
            A = [A; -[ones(K,1) obj.getCoeff()]];
            
            b = [b; -ones(K,1)];
            
            Aeq = [0 ones(1,D)];
            
            beq = 1;
            
            lb = zeros(D+1,1);
            
            %ub = ones(D+1,1);
            
            options = optimoptions(@linprog,'Algorithm','simplex');
            
            [x,Delta] = linprog(f,A,b,Aeq,beq,lb,[],[],options);
            
            %[x,Delta] = linprog(f,A,b,Aeq,beq,lb);
            
            ddout=x(2:end)';
        end
        
        function [best_dd,best_c,best_delta] = robustSolitonOpt(obj,n,numiter,flag)
            
            if nargin == 2
                numiter = 100;
                flag = 'BP';
                warning('Run 100 iterations for BP decoding!');
            elseif nargin == 3
                flag = 'BP';
                warning('Use BP decoding!');
            end
            
            switch flag
                case 'BP'
                    objfunc = @() obj.errorProb(obj.PoissonRec(n,'BP'));
                case 'inac'
                    objfunc = @() obj.expInac(obj.PoissonRec(n,'inac'));
            end
            
            iter = 0;
            best_delta = 0.5;
            K = obj.InputNumber;
            best_ei = 1;
            
            disp('Start to optimize the robust soliton distribution.');
            
            tic;
            
            while iter <= numiter
                newD = 0;
                max_c = sqrt(K)/log(K/best_delta)/2;
                min_c = 1/sqrt(K)/log(K/best_delta);
                for c = min_c:0.01:max_c
                    dd2 = solitonDist(c,best_delta);
                    obj.setDegreeDist(dd2);
                    ei = objfunc();
                    
                    if ei < best_ei
                        best_c = c;
                        best_dd = dd2;
                        best_ei = ei;
                        newD = 1;
                    end
                end
                
                min_delta = K/exp(sqrt(K)/best_c);
                max_delta = K/exp(1/best_c/sqrt(K));
                for delta = min_delta:0.05:max_delta
                    dd2 = solitonDist(best_c,delta);
                    obj.setDegreeDist(dd2);
                    ei = objfunc();
                    
                    if ei < best_ei
                        best_delta = delta;
                        best_dd = dd2;
                        best_ei = ei;
                        newD = 1;
                    end
                end
                
                if newD == 0
                    disp('No better degree dist found!');
                    break;
                else
                    disp(['New degree dist find (' num2str(best_c) ', ' num2str(best_delta), ', obj = ' num2str(best_ei),')!']);
                end
                
                iter = iter + 1;
            end
            
            timeused = toc;
            
            if iter > numiter
                iter = numiter;
            end
            
            disp(['Optimization runs for ' num2str(timeused) ' seconds:']);
            disp([num2str(iter) ' iter']);
            
        end
        
        function A = getCoeff(obj)
            D = obj.InputNumber;
            K = obj.InputNumber;
            A = zeros(K,D);
            
            TO = tril(ones(K,K));
                                     
            tmtau = kron(ones(1,K),(0:K-1)') - kron(ones(K,1),(0:K-1));
            as = tmtau*diag(1./(K:-1:1));
            
            for s = 1:obj.BatchSize
                % p_{[0:K-1],s} = Cs * Psi 
                Cs = zeros(K,D);
                Cs(1,s) = obj.hbar_sp(s+1);
                a = obj.hype_std(s+1,s+1:D,1:K-1);
                a = reshape(a,D-s,K-1);
                Cs(2:end,s+1:D) = obj.hbar_s(s+1)*a'*diag((s+1:D)/K);
                
                A = A + (TO - tril(as,-s))*Cs;
                
                if s < obj.BatchSize
                    tmtau = tmtau - 1;
                    as = as .* tmtau * diag([1./(K-s:-1:1) zeros(1,s)]);
                end
            end
        end
        
    end
    
    methods(Static)
        function degout = disturb(degin,d,delta)
            degout = degin;
            degout(d) = max(degout(d)+delta,0);
            degout = degout/sum(degout);
        end
    end
end