function X_path = geometric_brownian(N,r,sigma,T,startvalue)

t = (0:1:N)'/N;                   % t is the column vector [0 1/N 2/N ... 1]
W = [0; cumsum(randn(N,1))]/sqrt(N); % S is running sum of N(0,1/N) variables
t = t*T;
W = W*sqrt(T);
Y = (r-(sigma^2)/2)*t + sigma * W;
X = startvalue*exp(Y);


%%% Med valder:
% t = (0:1:N)'/N;                   % t is the column vector [0 1/N 2/N ... 1]
% W = [0; cumsum(randn(N,1))];%/sqrt(N); % S is running sum of N(0,1/N) variables
% t = mtimes(t,T);
% h = mrdivide(T,N);
% W = mtimes(W,sqrt(h));
% Y = plus(mtimes(minus(r,mrdivide(mpower(sigma,2),2)),t) , mtimes(sigma,W));
% X = startvalue*exp(Y);

X_path = X;

end