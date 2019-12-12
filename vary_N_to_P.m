% res1 is the number of different values for C_to_N that we try
res1 = 200;

% xres is the number of iterations used in finding best solution for 
% a given environment and a given rate of sythesis
xres = 200;

% sres is the number of iterations used in finding best solution for 
% a given environment and a given rate of sythesis
sres = 80;

N_to_P_min = 1;
N_to_P_max = 100;

C_to_N = 20;
kappa = 2;

% tau is the time in hours for hydrolases to digest their own mass
tau = 2;

% dry weight in grams per ml of the substrate
density = 0.5;

% Ci is the mass of carbon needed per unit volume of growth, in g per ml
Ci = 0.33;
Ni = 0.032;
Pi = 0.005;

% delta is the minimal fraction of C that must be digested in order to 
% digest the available N
delta = 0;

% epsilon is the efficiency of recycling for autolytic cells
epsilon = 0.5;

% alpha is the mass of machinery needed for cell mobility, relative
% to the mass of essential machinery
alpha = 0.02;

% beta is the mass of material in vesicles, relative to the mass of
% the rest of the fungus, including hydrolases
beta = 0.1;

% lambda is the maximum rate of resource use per unit volume,
% in g per ml per hour
lambda = 0.3;

% L is the maximum rate of resource use per unit volume,
% relative to Ci + Ni + Pi
L = lambda/(Ci + Ni + Pi);

% phi is the correction term
phi = (Ci + Ni)/(Ci + Ni + Pi);

N_to_P_vector = zeros(res1, 1);

% apparent growth rates for each category of organism
Mu_immobile = zeros(res1, 1);
Mu_motile = zeros(res1, 1);
Mu_autolytic = zeros(res1, 1);
Mu_fungal = zeros(res1, 1);

% rate of digestion per unit volume for each category of organism
D_immobile = zeros(res1, 1);
D_motile = zeros(res1, 1);
D_autolytic = zeros(res1, 1);
D_fungal = zeros(res1, 1);
DC_fungal = zeros(res1, 1);
DP_fungal = zeros(res1, 1);
DN_fungal = zeros(res1, 1);

for i = 1:res1
    
    N_to_P = exp(log(N_to_P_min) + ...
        (log(N_to_P_max) - log(N_to_P_min))*(i-1)/(res1-1));
    
    N_to_P_vector(i) = N_to_P;
    
    M_tot = 12*C_to_N*N_to_P + 14*N_to_P + 31;
    
    Ce = density*(kappa^2)*12*C_to_N*N_to_P/M_tot;
    Ne = density*(kappa^2)*14*N_to_P/M_tot;
    Pe = density*(kappa^2)*31/M_tot;
    
    [M, x] = find_best_immobile...
    (Ci, Ni, Pi, Ce, Ne, Pe, kappa, tau, delta, L, sres, xres);
     
    Mu_immobile(i) = M;
    D_immobile(i) = phi*x/tau;
    
    [M, x] = find_best_motile(Ci, Ni, Pi, Ce, Ne, Pe, ...
        kappa, tau, delta, L, alpha, sres, xres);
    
    Mu_motile(i) = M;
    D_motile(i) = phi*x/tau;
    
    [M, ~, x] = find_best_autolytic(Ci, Ni, Pi, Ce, Ne, Pe, ...
        kappa, tau, delta, L, epsilon, sres, xres);

    Mu_autolytic(i) = M;
    D_autolytic(i) = phi*x/tau;
    
    [M, x, xC, xN, xP] = find_best_fungi...
        (Ci, Ni, Pi, Ce, Ne, Pe, kappa, tau, delta,L,beta,sres,10*xres);
    
    Mu_fungal(i) = M;
    DC_fungal(i) = phi*xC/tau;
    DN_fungal(i) = phi*xN/tau;
    DP_fungal(i) = phi*xP/tau;
    D_fungal(i) = phi*x/tau;
    
    clc
    percent_finished = 100*i/res1
end

figure(3)
semilogx(N_to_P_vector, Mu_motile*24, 'c', 'LineWidth', 2)
hold on
semilogx(N_to_P_vector, Mu_autolytic*24, 'g', 'LineWidth', 2)
semilogx(N_to_P_vector, Mu_fungal*24, 'm', 'LineWidth', 2)
semilogx(N_to_P_vector, Mu_immobile*24, 'k--', 'LineWidth', 2)
xlabel('N:P ratio')
ylabel('Apparent growth rate, day^{-1}')
legend('Motile cells', 'Autolytic cells', 'Fungi', 'Immobile cells')
 
figure(4)
semilogx(N_to_P_vector, D_motile*(Ci+Ni+Pi)*24, 'c', 'LineWidth', 2)
hold on
semilogx(N_to_P_vector, D_autolytic*(Ci+Ni+Pi)*24, 'g', 'LineWidth', 2)
semilogx(N_to_P_vector, D_fungal*(Ci+Ni+Pi)*24, 'm', 'LineWidth', 2)
semilogx(N_to_P_vector, D_immobile*(Ci+Ni+Pi)*24, 'k--', 'LineWidth', 2)
xlabel('N:P ratio')
ylabel('Digestion Rate, g ml^{-1} day^{-1}')
legend('Motile cells', 'Autolytic cells', 'Fungi', 'Immobile cells')
