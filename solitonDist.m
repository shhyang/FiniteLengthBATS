function y = solitonDist(K,c,delta)
    y = [1/K 1./((2:K).*(1:(K-1)))];
    
    if nargin == 1
        return;
    elseif nargin ~= 3
        warning('One or three inputs!');
        return;
    end
    
    S = c*log(K/delta)*sqrt(K);
    
    D = ceil(K/S);
    
    tau = [S/K./(1:D-1) S/K*log(S/delta) zeros(1,K-D)];
    
    y = y + tau;
    
    y = y/sum(y);
end