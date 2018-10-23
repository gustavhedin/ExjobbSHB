clear all

UNDERLYING = [];
OPTIONPRICE = [];
DELTA = [];
VEGA = [];
RHO = [];

for s = 5:15
    %inputs:
    % volatility sigma
    % interest rate r
    % maturity time tau

    sigma      = valder(0.2,[1 0 0 0]);
    r       =  valder(0.05, [0 1 0 0]);
    tau      =     valder(1,[0 0 1 0]);
    startvalue = s; % Strike price

    K = 10; 

    % constants:
    N = 100; % # steps
    nbrMC = 10000; % # MC simulations
    X = valder(0);

    for i = 1:nbrMC
        gbm = geometric_brownian(N,r,sigma,tau,startvalue);
        path = valder(gbm.val,[gbm.der(:,1) gbm.der(:,2) gbm.der(:,3) ones(N+1,1)]);
        P = payoff(path,K);
        X = X + P;
    end


    X = mrdivide(X,nbrMC);

    output = double(X);
    %disp(['payoff Vega Rho -Theta Delta']);
    output(end,:);

    t = ((0:1:N)'/N)*tau.val;
    t_ = flipud(t);
    optionprice = X.val(end)*exp(-r.val*t_);

%     figure
%     plot(X.val,'o--r')
%     hold on
%     plot(optionprice,'m')
%     plot(X.der(:,1),'g')
%     plot(X.der(:,2),'b')
%     plot(X.der(:,3),'c')
%     plot(X.der(:,4),'k')
%     legend('payoff','optionprice','Vega','Rho','Theta')
% 
%     formatSpec = 'Optionvalue = %4.8f, Vega= %4.8f, Rho= %4.8f, Theta= %4.8f';
%     fprintf(formatSpec,optionprice(1),X.der(end,1),X.der(end,2),X.der(end,3))

    UNDERLYING = [UNDERLYING startvalue];
    OPTIONPRICE = [OPTIONPRICE optionprice(1)];
    DELTA = [DELTA X.der(end,4)];
    VEGA = [VEGA X.der(end,1)];
    RHO = [RHO X.der(end,2)];

end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    plot(UNDERLYING,OPTIONPRICE,'o--c')
    hold on
    plot(UNDERLYING,DELTA,'o--r')
    plot(UNDERLYING,VEGA,'o--m')
    plot(UNDERLYING,RHO,'o--b')
    %legend('OPTIONPRICE','DELTA','VEGA','RHO')

    d1 = (log(UNDERLYING./K) + (r.val+0.5*sigma.val^2)*(tau.val))./(sigma.val*sqrt(tau.val));
    d2 = d1 - sigma.val*sqrt(tau.val);

    delta_verification = cdf('Normal',d1,0,1);
    vega_verification = UNDERLYING.*pdf('Normal',d1,0,1).*sqrt(tau.val);
    rho_verification = K*tau.val*exp(-r.val*tau.val).*cdf('Normal',d2,0,1);
    
    plot(UNDERLYING,delta_verification,'k')
    plot(UNDERLYING,vega_verification,'y')
    plot(UNDERLYING,rho_verification,'g')
    legend('OPTIONPRICE','DELTA','VEGA','RHO','delta_ver','vega_ver','rho_ver')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
function P = payoff(X,K)
    P = X;
    for i = size(P.val,1)%1:size(P.val,1)
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
