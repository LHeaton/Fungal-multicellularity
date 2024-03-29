% res1 is the number of different values for the C:N ratio and 
% recalcitrance that are tried
res1 = 200;

% xres determines the number of different values for x that are tried in
% finding the optimal solution. Final value found is accurate to
% xres^2
xres = 200;

% res is the number of iterations used in finding best solution for 
% a given environment and a given rate of sythesis
sres = 60;

N_to_P = 50;
    
kappa = 2;
    
C_to_N_min = 5;
C_to_N_max = 300;
 
tau_min = 0.05;
tau_max = 50;

% Ci is the mass of carbon needed per unit volume of growth, in g per ml
Ci = 0.33;
Ni = 0.032;
Pi = 0.005;

% dry weight in grams per ml of the substrate
density = 0.5;

% epsilon is the efficiency of recycling for senscent cells
epsilon = 0.9;

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

delta = 0;
    
C_to_N_vector = zeros(res1, 1);
Recalcitrance = zeros(res1, 1);

Mu_immobile = zeros(res1);
Mu_motile = zeros(res1);
Mu_autolytic = zeros(res1);
Mu_fungal = zeros(res1);
Mu_cell = zeros(res1);

x_immobile = zeros(res1);
x_motile = zeros(res1);
x_autolytic = zeros(res1);
x_fungal = zeros(res1);
xC_fungal = zeros(res1);
xN_fungal = zeros(res1);
xP_fungal = zeros(res1);
x_cell = zeros(res1);

for i = 1:res1
    
    C_to_N = exp(log(C_to_N_min) + ...
        (log(C_to_N_max) - log(C_to_N_min))*(i-1)/(res1-1));
    
    C_to_N_vector(i) = C_to_N;
    
    M_tot = 12*C_to_N*N_to_P + 14*N_to_P + 31;
    
    Ce = density*(kappa^2)*12*C_to_N*N_to_P/M_tot;
    Ne = density*(kappa^2)*14*N_to_P/M_tot;
    Pe = density*(kappa^2)*31/M_tot;
        
    for j = 1:res1
        
        tau = exp(log(tau_min) + ...
            (log(tau_max) - log(tau_min))*(j-1)/(res1-1));
        
        Recalcitrance(j) = tau;
        
        [Mu_cell(i,j), x_cell(i,j)] = ...
            find_best_cell(Ci, Ni, Pi, Ce, Ne, tau, delta, L);
        
        [Mu_immobile(i,j), x_immobile(i,j)] = find_best_immobile...
            (Ci, Ni, Pi, Ce, Ne, Pe, kappa, tau, delta, L, sres, xres);
        
        [Mu_motile(i,j), x_motile(i,j)] = find_best_motile...
            (Ci, Ni, Pi, Ce, Ne, Pe, kappa, tau, delta, ...
            L, alpha, sres, xres);
        
        [Mu_autolytic(i,j), ~, x_autolytic(i,j)] = find_best_autolytic...
            (Ci, Ni, Pi, Ce, Ne, Pe, kappa, tau, delta, ...
            L, epsilon, sres, xres);
        
        [Mu_fungal(i,j), x_fungal(i,j), xC_fungal(i,j), ...
            xN_fungal(i,j), xP_fungal(i,j)] = find_best_fungi(Ci, Ni, ...
            Pi, Ce, Ne, Pe, kappa, tau, delta, L, beta, sres, 2*xres);
    end
    
    percent_finished = 100*i/res1
end

clear tau
clear M_tot
clear percent_finished

save('vary_CNP_and_tau_NP50_k2_beta05.mat')
