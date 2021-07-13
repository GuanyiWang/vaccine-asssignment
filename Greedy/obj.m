function f = obj(V,W,C,g)

%the default choice of weight is equal-weight
if nargin < 4, g = ones(length(V),1); end
    
%Define the objective function
N = length(V);
f = (V'*W*V + C'*V - ones(1,N)* W * V - V'*W*ones(N,1))/sum(g);
end