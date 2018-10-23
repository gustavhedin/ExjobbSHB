clear all

UNDERLYING = [];
OPTIONPRICE = [];
DELTA = [];
VEGA = [];
RHO = [];


for s=5:15
    % inputs:
    % volatility sigma
    % interest rate r
    % maturity time tau
    sigma_der = 0;
    r_der = 0;
    tau_der = 0;

    startvalue = s; 

    K = 10; % Strike price

    % constants:
    N = 100; % # steps
    nbrMC = 10000; % # MC simulations
    X = zeros(N+1,1);

    for i = 1:nbrMC
        sigma      = ADRev(0.2);
        r       =   ADRev(0.05);
        tau      =     ADRev(1);

        eRev_forward = geometric_brownian_adr(N,r,sigma,tau,startvalue,K);
        reverse_sweep = chainRule(eRev_forward); 
        X = X + reverse_sweep.value;
        sigma_der = sigma_der + sigma.derivative(end);
        r_der = r_der + r.derivative(end);
        tau_der = tau_der + tau.derivative(end);
    end


    X = X/nbrMC;
    sigma_der = sigma_der/nbrMC;
    r_der = r_der/nbrMC;
    tau_der = tau_der/nbrMC;
    
    optionprice = X(end)*exp(-r.value*tau.value);

    
%     figure
%     plot(X,'o--r')
%     hold on
%     plot(optionprice,'m')
%     plot(sigma_der,'g')
%     plot(r_der,'b')
%     plot(tau_der,'c')
%     legend('payoff','optionprice','Vega','Rho','Theta')
% 
%     formatSpec = 'Optionvalue = %4.8f, Vega= %4.8f, Rho= %4.8f, Theta= %4.8f\n';
%     fprintf(formatSpec,optionprice(1),sigma_der(end),r_der(end),tau_der(end))

    UNDERLYING = [UNDERLYING startvalue];
    OPTIONPRICE = [OPTIONPRICE optionprice(1)];
    DELTA = [DELTA 0];
    VEGA = [VEGA sigma_der];
    RHO = [RHO r_der];

end

    RHO = (RHO-tau.value*OPTIONPRICE)*exp(-r.value*tau.value);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    plot(UNDERLYING,OPTIONPRICE,'o--c')
    hold on
    plot(UNDERLYING,DELTA,'o--r')
    plot(UNDERLYING,VEGA,'o--m')
    plot(UNDERLYING,RHO,'o--b')
    %legend('OPTIONPRICE','DELTA','VEGA','RHO')

    d1 = (log(UNDERLYING./K) + (r.value+0.5*sigma.value^2)*(tau.value))./(sigma.value*sqrt(tau.value));
    d2 = d1 - sigma.value*sqrt(tau.value);

    delta_verification = cdf('Normal',d1,0,1);
    vega_verification = UNDERLYING.*pdf('Normal',d1,0,1).*sqrt(tau.value);
    rho_verification = K*tau.value*exp(-r.value*tau.value).*cdf('Normal',d2,0,1);
    
    plot(UNDERLYING,delta_verification,'k')
    plot(UNDERLYING,vega_verification,'y')
    plot(UNDERLYING,rho_verification,'g')
    legend('OPTIONPRICE','DELTA','VEGA','RHO','delta_ver','vega_ver','rho_ver')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

