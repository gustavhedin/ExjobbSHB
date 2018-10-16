clear
% inputs:
% volatility sigma
% interest rate r
% maturity time tau
sigma_der = 0;
r_der = 0;
tau_der = 0;

startvalue = 10; % Strike price

K = 10; 

% constants:
N = 100; % # steps
nbrMC = 10000; % # MC simulations
X = zeros(N+1,1);

for i = 1:nbrMC
    sigma      = ADRev(0.2);
    r       =  ADRev(0.05);
    tau      =     ADRev(1);
    
    eRev_forward = geometric_brownian_adr(N,r,sigma,tau,startvalue,K);
    reverse_sweep = chainRule(eRev_forward); 
    X = X + reverse_sweep.value;
    sigma_der = sigma_der + sigma.derivative;
    r_der = r_der + r.derivative;
    tau_der = tau_der + tau.derivative;
end


X = X/nbrMC;
sigma_der = sigma_der/nbrMC;
r_der = r_der/nbrMC;
tau_der = tau_der/nbrMC;

% disp('payoff')
% X(end,:)

t = ((0:1:N)'/N)*tau.value;
t_ = flipud(t);
optionprice = X(end)*exp(-r.value*t_);

figure
plot(X,'o--r')
hold on
plot(optionprice,'m')
plot(sigma_der,'g')
plot(r_der,'b')
plot(tau_der,'c')
legend('payoff','optionprice','Vega','Rho','Theta')

formatSpec = 'Optionvalue = %4.8f, Vega= %4.8f, Rho= %4.8f, Theta= %4.8f';
fprintf(formatSpec,optionprice(1),sigma_der(end),r_der(end),tau_der(end))


%-------------------------------------------------------------------------%
function P = payoff(X,K)
    P = X;
    for i = 1:size(P.val,1)
        if X.val(i)-K > 0
            P.val(i) = X.val(i)-K;
        else
            P.val(i) = 0;
            P.der(i,1) = 0;
            P.der(i,2) = 0;
            P.der(i,3) = 0;
            P.der(i,4) = 0;
        end
    end
end

function P = payoff_barrier(X,K)
    P = X;
    
    if X.val(end) > K
        P = X;
    else
       for i = 1:size(P.val,1)
                P.val(i) = 0;
                P.der(i,1) = 0;
                P.der(i,2) = 0;
                P.der(i,3) = 0;
       end
    end
end