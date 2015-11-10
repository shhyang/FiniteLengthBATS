classdef BATSBasic < handle
    properties (SetAccess = private)
        BatchSize
        InputNumber
        FieldSize
        RankDist
    end
    properties (Hidden, SetAccess = private)
        hbar_s % Effective Rank Distribution
        hbar_sp % accumulative sum of effective rank distribution
        hype_std % hype_std(s+1,d,t) = hype(K-1,d-1,t-1,d-s-1)
%        h_Q % h_Q(t+1,i+1,s+1,j-i+1) = hype(K-t,i,s,i+s-j)
    end
    
    methods
        function BB = BATSBasic(M,K,q,h)
            % M : batch size
            % K : number of input
            % q : field size
            % h : rank distributions in a matrix form. Each row of h is a rank
            % distribution. The number of columns of h is M+1. The first column is for
            % the probability of rank zero and the last column is for the probability
            % of rank M, full rank.
            %
            BB.BatchSize = M;
            BB.InputNumber = K;
            BB.FieldSize = q;
            
            BB.setRankDist(h);
            BB.hyper(K,M,K);
        end
        
        function setRankDist(BB,h)
            Mp = BB.BatchSize+1;
            if length(h)>Mp
                error('Length of rank distribution is larger than batch size + 1.');
            elseif length(h)<Mp
                error('Length of rank distribution is smaller than batch size + 1.');
            else 
                a = h;
            end
            
            BB.RankDist = a;
            
            if BB.BatchSize == 1
                BB.hbar_s = a;
                BB.hbar_sp = [1,a(2)];
            else
                BB.hbar_s = zeros(size(a));
                BB.hbar_sp = BB.hbar_s;
                for s = 0:BB.BatchSize
                    for i = s:BB.BatchSize
                        tmp = BATSBasic.zeta_mrq(i,s,BB.FieldSize)*a(i+1);
                        BB.hbar_sp(s+1) = BB.hbar_sp(s+1) + tmp;
                        BB.hbar_s(s+1) = BB.hbar_s(s+1) + tmp/(BB.FieldSize)^(i-s);
                    end
                end
            end
            
        end
        
        function hyper(obj,K,M,D)
            % hype(K-1,d-1,t-1,d-s-1), here max{1,d+t-K} <= d-s <= t, s = 0,...,M
            obj.hype_std = zeros(M+1,D,K);
            for s = 0:M
                tmp = 1; % hype(K-1,s,0,0) = 1;
                for d = s+1:D
                    tmp2 = tmp; % hype(K-1,d-1,d-s-1,d-s-1) = (d-1 d-s-1)/(K-1 d-s-1) = prod_{i=0}^{d-s-2} (d-1-i)/(K-1-i)
                    for t = d-s:K-s
                        obj.hype_std(s+1,d,t) = tmp2; % hype(K-1,d-1,t-1,d-s-1)
                        tmp2 = tmp2 * t/(K-t) * (K-t-s)/(t-d+s+1); % hype(K-1,d-1,t,d-s-1) = t/(K-t) * (K-t-s)/(t-d+s+1) * hype(K-1,d-1,t-1,d-s-1)
                    end
                    tmp = tmp * d/(K-d+s); % hype(K-1,d,d-s,d-s) = (d d-s)/(K-1 d-s) = prod_{i=0}^{d-s-1} (d-i)/(K-1-i) = d/(K-d+s) * (d-1 d-s-1)/(K-1 d-s-1)
                end
            end
        end
        
        
        function A = getpmatrix(obj)
            K = obj.InputNumber;
            M = obj.BatchSize;
            
            A = zeros(K+1,K);
            
            % p_0
            A(1,1:M) = obj.hbar_sp(2:end); 
            
            % p_t, t>0
            
            for s = 1:obj.BatchSize
                % p_s^t = hbar_s *
                a = obj.hype_std(s+1,s+1:K,1:K);
                a = reshape(a,K-s,K);
                a = a';
                A(2:end, s+1:K) = A(2:end, s+1:K) + obj.hbar_s(s+1) * a ;
            end
            A(2:end,:) = A(2:end,:) * diag((1:K)/K);
        end
        
    end
    
    
    
    methods(Static)
        function y = zeta_mrq(m,r,q)
            if m < r
                tmp = m;
                m = r;
                r = tmp;
            end
            if (r == 0)
                y = 1;
            else
                y = 1;
                for i = 0:r-1
                    y = y * (1 - q^(-m+i));
                end
            end
        end
    end
    
end