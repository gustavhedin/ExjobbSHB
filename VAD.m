
% User defined parameters:
r = 0.05;           % risk free interest rate
sigma = 0.2;        % volatility
Startvalue = 10;    % Starting value for the underlying asset at time 0.
T = 1;              % Time horizon
N = 1000;           % # simulation points on [0,T];
K = 10;             % Strike price
nbr_MC = 10000;         % # of Monte Carlo simulations
nbrMC_z = 10;       % # of samples over the barrier

% Declaring nedded variables:
h = T/N;

X = zeros(N,nbr_MC);
X(1,:) = Startvalue*ones(1,nbr_MC);

Y_delta = zeros(N,nbr_MC);
Y_delta(1,:) = ones(1,nbr_MC);

Y_vega = zeros(N,nbr_MC);

 for n=2:N
         Z = randn(1,nbr_MC);
        % Asset
        X(n,:) = X(n-1,:) + r*h*X(n-1,:) + sigma*sqrt(h)*X(n-1,:).*Z;

        % Tangent Processes: 
        % Delta: theta = X0:
        Y_delta(n,:) = Y_delta(n-1,:) + r*h*Y_delta(n-1,:) + 0*X(n-1,:) ...
            + (Y_delta(n-1,:)*sigma*sqrt(h)+ 0*X(n-1,:)).*Z; 
        
        % Vega: theta = sigma 
        Y_vega(n,:) = Y_vega(n-1,:) + r*h*Y_vega(n-1,:) + 0*X(n-1,:) ...
            + (Y_vega(n-1,:)*sigma*sqrt(h)+ sqrt(h)*X(n-1,:)).*Z; 
end
 
 % Find price 
 Z = randn(1,nbr_MC);
 Xend = X(end,:)+r*h*X(end,:) + sigma*sqrt(h)*X(end,:).*Z;
 Price = mean(max(Xend-K,0))*exp(-r*T);   
 
 % Generate random nubers for the step over the barrier.
 Z = randn(nbrMC_z,nbr_MC);
    
 % Obtain valued for the process and derivatives on the other side of
 % the barrier:
   
 firstpart = repmat(X(end,:)+r*h*X(end,:),nbrMC_z,1);
 lastpart = repmat(sigma*X(end,:)*sqrt(h),nbrMC_z,1);
 % First derivatives:
    X_Tplus = firstpart + Z.*lastpart;
    X_Tminus = firstpart - Z.*lastpart;
    X_Tdot = firstpart;

    V_Tplus = payoff(X_Tplus,K);
    V_Tminus = payoff(X_Tminus,K);
    V_Tdot = payoff(X_Tdot,K);

    % Delta
     dmu_dtheta = repmat(Y_delta(end,:)*(1+r*h) + X(end,:)*0,nbrMC_z,1);
     dsig_dtheta = repmat(Y_delta(end,:)*sigma*sqrt(h) + X(end,:)*0,nbrMC_z,1);
     divfactor =repmat(1./(X(end,:)*sigma*sqrt(h)),nbrMC_z,1);
     Delta = mean(mean(dmu_dtheta.*(1/2).*(V_Tplus-V_Tminus).*(Z.*divfactor)...
            + dsig_dtheta.*(V_Tplus-2*V_Tdot+V_Tminus).*((Z.^2-1).*divfactor)))*exp(-r*T);

    % Vega
      dmu_dtheta =  repmat(Y_vega(end,:)*(1+r*h) + X(end,:)*0,nbrMC_z,1);
      dsig_dtheta =  repmat(Y_vega(end,:)*sigma*sqrt(h) + X(end,:)*sqrt(h),nbrMC_z,1);
      Vega = mean(mean(dmu_dtheta.*(1/2).*(V_Tplus-V_Tminus).*(Z.*divfactor)...
            + dsig_dtheta.*(V_Tplus-2*V_Tdot+V_Tminus).*((Z.^2-1).*divfactor)))*exp(-r*T);
 
    % Here, apply an AD-tequnique to get 2nd order derivatives. 

    
% Here, redo simulation until 1 step before next barrier:
X = [X;Xend];

% payoff = V(end)*prodsum(I1, I2 ...) , Ii = 1 if X > K else 0, at a specific barrier checkingpoint.  

%% Verification:

tt = 0;
d1 = (log(Startvalue/K) + (r+0.5*sigma^2)*(T-tt))./(sigma*sqrt(T-tt));
d2 = d1 - sigma*sqrt(T-tt);

delta_verification = cdf('Normal',d1,0,1);
vega_verification = Startvalue.*pdf('Normal',d1,0,1).*sqrt(T-tt);
price_verification = Startvalue.*cdf('Normal',d1,0,1) - exp(-r*(T-tt)).*K.*cdf('Normal',d2,0,1);


formatSpec = 'Optionvalue = %4.8f, Vega= %4.8f, Rho= %4.8f, Theta= %4.8f';
fprintf(formatSpec,price_verification(1),vega_verification(1),0,0)

%delta_verification
%Delta
%vega_verification 
%Vega
%price_verification
%Price
%%
%-------------------------------------------------------------------------%
function P = payoff(x,K)
   P=max(x-K,0);
end

