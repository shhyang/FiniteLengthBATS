classdef BATSAsymp < BATSBasic & handle
   
    properties
        DegreeDist
        Rate
        DegreeDist_minrate
        Rate_minrate
        DegreeDist_maxratio
        Rate_maxratio
    end
    
    methods
        function obj = BATSAsymp(M,K,q,h)
            obj = obj@BATSBasic(M,K,q,h);
        end
        
        function [ddist,rate,exitflag] = opt(obj,eta)
            % This function generates the degree distribution for a single input rank
            % distribution
            %
            % eta : fraction of unrecovered packets, e.g., 0.01
            %
            
            z = obj.hbar_s(:,2:end);
            
            [ddist,rate,exitflag] = obj.degreeDistForEffectiveRank(eta, z(1,:));
                    
            obj.DegreeDist = ddist;
            obj.Rate = rate;
            
        end
           
        function [ddist,rate,exitflag] = optMulti_minrate(obj,eta)
            % This function generates the degree distribution for multiple rank
            % distributions
            %
            % eta : fraction of unrecovered packets, e.g., 0.01
            %
            
            z = obj.hbar_s(:,2:end);
            
            [ddist,rate,exitflag] = obj.degreeDistForMinimalRate(eta, z);
            
            obj.DegreeDist_minrate = ddist;
            obj.Rate_minrate = rate;
        end
        
        function [ddist,rate,exitflag] = optMulti_maxratio(obj,eta)
            % This function generates the degree distribution for multiple rank
            % distributions
            %
            % eta : fraction of unrecovered packets, e.g., 0.01
            %
            
            z = obj.hbar_s(:,2:end);
            ehb = z * (1:M)';
            u = min(ehb)/2;
            z = diag(u./ehb)*z;
            [ddist,rate,exitflag] = obj.degreeDistForMinimalRate(eta, z);
            rate = rate/u;
            
            obj.DegreeDist_maxratio = ddist;
            obj.Rate_maxratio = rate;              
        end
                
        function [x,fval,exitflag] = degreeDistForEffectiveRank(obj, eta, z)
            
            M = obj.BatchSize;
            
            Ndd = 100; % number of sample
            
            v = linspace(0,1-eta,Ndd);
            
            v = v';
            
            D = min(obj.InputNumber, ceil(M/eta)-1);
            
            % This function uses linear programming to find the degree distribution
            %
            
            % z is the effective rank distribution of H
            % z = z(1:M)
            % D is the maximum degree
            
            %
            % min_x f^T x such that
            % A x <= b,
            % Aeq x = beq,
            % lb <= x <= ub.
            tic
            
            f = [-1 zeros(1,D)]';
            
            Aeq = [0 ones(1,D)];
            beq = 1;
            
            % generate A, b
            
            b = zeros(length(v),1);
            
            
            A = getAforLinearConstraints();
            
            % for testing only
            %A(:,2) = 0;
            % %%%%%%%%%%%%%%%
            
            lb = zeros(D+1,1);
            % ub = 2*[M;ones(D,1)];
            
            
            options = optimoptions('linprog','Algorithm','simplex');
            %options = optimset('LargeScale','off','Simplex','on');
            %options = optimset('LargeScale','off');
            
            [x,fval,exitflag] = linprog(f,A,b,Aeq,beq,lb,[],[],options);
            
            %[x,fval] = linprog(f,A,b,Aeq,beq,lb,[],[],optimset('Display','iter','MaxIter',200));
            toc
            
            fval = -fval;
            x = x(2:end);
            x = x'/sum(x);
            
            
            function A = getAforLinearConstraints()
                
                A = zeros(length(v),D+1);
                
                sz = zeros(1,M);
                
                for r = 1:M
                    sz(r) = sum(z(r:M),2);
                end
                
                % first column of A
                
                theta = 1;
                
                A(:,1) = theta * log(1-v);
                
                % other columns
                
                for d = 1:M
                    for r = 1:d-1
                        A(:,d+1) = A(:,d+1) + z(r)*betainc(v,d-r,r);
                    end
                    A(:,d+1) = A(:,d+1) + sz(:,d)*ones(size(v));
                    A(:,d+1) = A(:,d+1)*d;
                end
                
                for d = (M+1) : D
                    for r = 1:M
                        A(:,d+1) = A(:,d+1) + z(r)*betainc(v,d-r,r);
                    end
                    A(:,d+1)=A(:,d+1)*d;
                end
                A = -A;
            end
            
        end
        
        function [x,fval,exitflag] = degreeDistForMinimalRate(obj, eta, z)
            
            M = obj.BatchSize;
            
            Ndd = 100; % number of sample
            
            
            v = linspace(0,1-eta,Ndd);
            
            v = v';
            
            D = min(obj.InputNumber,ceil(M/eta)-1);
            
            % This function uses linear programming to find the degree distribution
            %
            
            % z is the effective rank distribution of H
            % z = z(1:M)
            % D is the maximum degree
            
            %
            % min_x f^T x such that
            % A x <= b,
            % Aeq x = beq,
            % lb <= x <= ub.
            tic
            
            f = [-1 zeros(1,D)]';
            
            
            
            Aeq = [0 ones(1,D)];
            beq = 1;
            
            % generate A, b
            N = size(z,1);
            
            b = zeros(N*length(v),1);
            
            
            A = getAforLC();
            
            % for testing only
            %A(:,2) = 0;
            % %%%%%%%%%%%%%%%
            
            lb = zeros(D+1,1);
            % ub = 2*[M;ones(D,1)];
            
            options = optimoptions('linprog','Algorithm','simplex');
            %options = optimset('LargeScale','off','Simplex','on');
            %options = optimset('LargeScale','off');
            
            [x,fval,exitflag] = linprog(f,A,b,Aeq,beq,lb,[],[],options);
            
            %[x,fval] = linprog(f,A,b,Aeq,beq,lb,[],[],optimset('Display','iter','MaxIter',200));
            toc
            
            fval = -fval;
            x = x(2:end);
            x = x'/sum(x);
            
            
            
            function A = getAforLC()
                A = zeros(length(v)*N,D+1);
                
                sz = zeros(N,M);
                
                for r = 1:M
                    sz(:,r) = sum(z(:,r:M),2);
                end
                
                % first column of A
                
                theta = ones(N,1);
                
                A(:,1)=kron(theta,log(1-v));
                
                % other columns
                
                for d = 1:M
                    for r = 1:d-1
                        A(:,d+1) = A(:,d+1) + kron(z(:,r),betainc(v,d-r,r));
                    end
                    A(:,d+1) = A(:,d+1) + kron(sz(:,d),ones(size(v)));
                    A(:,d+1) = A(:,d+1)*d;
                end
                
                for d = (M+1) : D
                    for r = 1:M
                        A(:,d+1) = A(:,d+1) + kron(z(:,r),betainc(v,d-r,r));
                    end
                    A(:,d+1)=A(:,d+1)*d;
                end
                A = -A;
            end
            
        end
        
    end
end