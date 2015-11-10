classdef BATSFiniteLength < BATSBasic & handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        BatchNumber
        DegreeDist
        ExpBatchNumber
    end
    
    properties (Hidden)
        p % p(t+1,s+1) = p_s^t
        psum % 
        rho
%        Q % Q(t+1,i,j) = Q_t(i,j)
        r_BP
    end
    
    properties (SetAccess = private)

    end
    
    methods
        function obj = BATSFiniteLength(M,K,q,h)
            obj = obj@BATSBasic(M,K,q,h);
        end
        
        function y = Q_t(obj,t)
            if obj.psum(t+1) == 0
                y = getQt(obj.InputNumber - t,obj.BatchSize, 0);
            else
                y = getQt(obj.InputNumber - t,obj.BatchSize,obj.p(t+1,:)/obj.psum(t+1));
            end
        end
        
        function y = Qp_t(obj,t)
            % this function return p_t Q_t
            y = getQt(obj.InputNumber - t,obj.BatchSize,obj.p(t+1,:));
        end
        
        function [Delta,y] = errorExponent(obj)
            y = 1 - cumsum(obj.psum);
            for tau = 0:obj.InputNumber
                lam = diag(obj.Qp_t(tau));
                y(tau+1:end) = y(tau+1:end) + lam(1:obj.InputNumber+1-tau);
            end
            %y = -log(y);
            y(end)=0;
            y(2:obj.r_BP) = 0;
            Delta = -log(max(y));
        end     
                
        function stop = FixedRec(obj,n,flag)
            % flag = 'BP' : compute BP
            % flag = 'inac' : with inactivation
            
            if nargin == 2
                flag = 'BP';
            end
            
            if isempty(obj.DegreeDist)
                error('Set digree distribution first!');
            end
            
            if ~strcmp(flag,'BP') && ~strcmp(flag,'inac')
                error('flag = ''BP'' or ''inac''!');
            end
            
            K = obj.InputNumber;
            
            stop = zeros(1,K+1);
            
            % t = 0
            
            lambda = zeros(n+1,K+1);            
            rho_t = obj.rho(1);
            lambda(n+1,:) = [(1-rho_t)^n, zeros(1,K)];
            Qt = obj.Q_t(0);
            for c = n:-1:1
                lambda(c,:) = lambda(c+1,:) * Qt * c/(n-c+1) * rho_t/(1-rho_t);
            end
            
            stop(1) = sum(lambda(:,1));
            
            % t > 0
            for t=1:K
                rho_t = obj.rho(t+1);
                lambdatmp = zeros(n+1,K+1-t);
                bi = (1-rho_t).^(0:n);
                Qt = obj.Q_t(t);
                QtP = eye(K-t+1);
                
                for dc = 0:n
                    if strcmp(flag, 'BP')
                        lambdatmp(1:(n+1-dc),:) = lambdatmp(1:(n+1-dc),:) + diag(bi) * lambda((1+dc):(n+1),2:end) * QtP;
                    else 
                        lambdatmp(1:(n+1-dc),:) = lambdatmp(1:(n+1-dc),:) + diag(bi) * [lambda((1+dc):(n+1),1)+lambda((1+dc):(n+1),2) lambda((1+dc):(n+1),3:end)] * QtP;
                    end
                    bi(end)=[];
                    bi = bi .* (dc+1+(0:n-dc-1))/(dc+1)*rho_t;
                    QtP = QtP*Qt;
                end
                
                lambda = lambdatmp;
                stop(t+1) = sum(lambda(:,1));
            end
            
        end

        function stop = FixedRec_acc(obj,N,flag)
            % flag = 'BP' : compute BP
            % flag = 'inac' : with inactivation
            % This function output the stopping time distribution for n =
            % 0,...,N
            
            if nargin == 2
                flag = 'BP';
            end
            
            if isempty(obj.DegreeDist)
                error('Set digree distribution first!');
            end
            
            if ~strcmp(flag,'BP') && ~strcmp(flag,'inac')
                error('flag = ''BP'' or ''inac''!');
            end
            
            K = obj.InputNumber;            
            stop = zeros(N+1,K+1);            
            
            % initialize the first rows
            Lambda = zeros(N+1,K+1);
            Qtp = obj.Qp_t(0);
            Lambda(1,1) = 1;
            
            for n = 1:N
                Lambda(n+1,:) =  Lambda(n,:)*Qtp;
                %Lambda(n+1,:) =  mulQt(Lambda(n,:),Qtp,obj.BatchSize);
            end
            
            % start the loop for 0 <= t < K
            prodrhot = 1;
            rhov = ones(1,N+1);
            
            for t=0:K
                prodrhot = prodrhot*(1 - obj.rho(t+1));                             
                for i = 0:N
                     stop(i+1,t+1) = rhov(N+1-i:end) * Lambda(1:i+1,t+1);
                     if i < N
                        rhov(N-i) = rhov(N+1-i) * prodrhot;
                        rhov(N+1-i:end) = rhov(N+1-i:end) * (i+1)./(1:i+1);
                     end
                end
 
                if t < K                   
                    Qtp = obj.Qp_t(t+1);                    
                    nchoosec = ones(N+1,1);                    
                    Lambda_t = Lambda(:,t+2:end);
                    if strcmp(flag,'inac')
                        Lambda_t(:,1) = Lambda_t(:,1) + Lambda(:,t+1);
                    end
                    Lambda(:,t+2:end) = 0;                    
                    for c = 0:N
                        Lambda(c+1:end,t+2:end) = Lambda(c+1:end,t+2:end) + diag(nchoosec) * Lambda_t;
                        if c < N
                            nchoosec(1) = [];
                            nchoosec = nchoosec .* (1:N-c)' / (c+1);
                            Lambda_t = Lambda_t(1:end-1,:) * Qtp;
                            %Lambda_t = mulQt(Lambda_t(1:end-1,:), Qtp, obj.BatchSize);
                        end
                    end                   
                end
            end                     
        end
        
        function stop = PoissonRec(obj,n,flag)
            % flag = 'BP' : compute BP
            % flag = 'inac' : with inactivation
            if nargin == 2
                flag = 'BP';
            end
            
            if isempty(obj.DegreeDist)
                error('Set digree distribution first!');
            end
            
            if ~strcmp(flag,'BP') && ~strcmp(flag,'inac')
                error('flag = ''BP'' or ''inac''!');
            end
            
            K = obj.InputNumber;
                        
            stop = zeros(1,K+1);
            
            Gamma = [1 zeros(1,K)];
            
            for t = 0:K
                QtmI = obj.Qp_t(t) - obj.psum(t+1)*eye(K-t+1);
                Gamma = (expmv(n,QtmI',Gamma'))';
                %Gamma = Gamma * expm(n*QtmI);
                
                stop(t+1) = Gamma(1);
                if strcmp(flag,'inac') && t < K
                     Gamma(2) = Gamma(1) + Gamma(2);
                end
                Gamma(1) = [];
            end           
        end
        
        function stop = PoissonRec_acc(obj,N,n0,flag)
            % flag = 'BP' : compute BP
            % flag = 'inac' : with inactivation
            % this function output the stopping time distribution for 
            % n = 1,...,N
            
            if nargin == 3
                flag = 'BP';
            elseif nargin == 2
                n0 = 1;
                flag = 'BP';
            end
            
            if isempty(obj.DegreeDist)
                error('Set digree distribution first!');
            end
            
            if ~strcmp(flag,'BP') && ~strcmp(flag,'inac')
                error('flag = ''BP'' or ''inac''!');
            end
            
            K = obj.InputNumber;
                                    
            stop = zeros(N,K+1);
            
            Gamma = [ones(N,1) zeros(N,K)];
            
            for t = 0:K
                QtmI = obj.Qp_t(t) - obj.psum(t+1)*eye(K-t+1);
                Q0 = expm(n0*QtmI);
                for n = 1:N
                    Gamma(n:N,:) = Gamma(n:N,:) * Q0;
                end
                
                stop(:,t+1) = Gamma(:,1);
                if strcmp(flag,'inac') && t < K
                     Gamma(:,2) = Gamma(:,1) + Gamma(:,2);
                end
                Gamma(:,1) = [];
            end
        end
                 
        function setDegreeDist(obj,Psi)
            obj.DegreeDist = Psi;
            if length(Psi) > obj.InputNumber
                obj.DegreeDist=Psi(1:obj.InputNumber)/sum(Psi(1:obj.InputNumber));
            elseif length(Psi) < obj.InputNumber
                obj.DegreeDist(obj.InputNumber) = 0;
            end

            % generate p_{t,s}, rho and p_t
            obj.p = zeros(obj.InputNumber+1,obj.BatchSize+1);
            
            % p_s^0 = hbar_s'Psi_s
            obj.p(1,2:end) = obj.hbar_sp(2:end).*obj.DegreeDist(1:obj.BatchSize);
            
            D = length(obj.DegreeDist);
            % 
            for s = 0:obj.BatchSize
                % p_s^t = hbar_s * 
                a = obj.hype_std(s+1,s+1:D,1:obj.InputNumber);
                a = reshape(a,D-s,obj.InputNumber);
                obj.p(2:end,s+1) = obj.hbar_s(s+1)*(obj.DegreeDist(s+1:D).*(s+1:D)/obj.InputNumber * a)';
            end
            
            % get rho
            obj.psum = sum(obj.p, 2);
            tmp = cumsum(obj.psum);
            obj.rho = obj.psum./[1; 1-tmp(1:end-1)];
            obj.rho(isnan(obj.rho)) = 1;
            
%             % set r_BP
%             obj.r_BP = find(Psi(1:obj.BatchSize) .* obj.hbar_sp(2:end)>0);
%             if ~isempty(obj.r_BP)
%                 obj.r_BP = obj.r_BP(1);
%             else
%                 obj.r_BP = 0;
%             end
        end
        
    end
    
    methods(Static)
       
        function [err, eco] = errorProb(stop_dist)
            err = sum(stop_dist(:,1:end-1),2);
            eco = cumsum(err);
        end
        function y = expInac(inact_p)
            y = sum(inact_p(:,1:end-1),2);
        end
    end
end

