function f = obj(V,W,C)
%Define the objective function
N = length(V);
f = (V'*W*V + C'*V - ones(1,N)* W * V - V'*W*ones(N,1))/N;
end