eta = 1;
res = 100000;
t_end = 10/eta;

t = linspace(0, t_end, res);
V = exp(eta*t);
mu = eta*ones(res, 3);

dt = t(2);

mu_lim = zeros(2, 3);

for dum = 1:3
    
    if dum == 1
        
        omega = 1.1;
        
        % tt is the number of time steps taken to exhaust 
        % the local supply
        tt = find(t >= omega/eta, 1, 'first');
        
    elseif dum == 2
        
        omega = 1.5;
        tt = find(t >= omega/eta, 1, 'first');
        
    elseif dum == 3
        
        omega = 3;
        tt = find(t >= omega/eta, 1, 'first');
    end
    
for i = 2:res
    
    if i > tt
    
        V(i) = V(i-1)*(1 + mu(i-1, dum)*dt);
        f = V(i-tt)/V(i);
        
        mu(i, dum) = eta*(1-f);
    end
end

mu_lim(:, dum) = eta*(1 + lambertw(-omega*exp(-omega))/omega); 
end

figure
plot(t, mu(:,1), 'r')
hold on
plot(t, mu(:,2), 'b')
plot(t, mu(:,3), 'g')
plot([t(1), t(end)], mu_lim(:,1), 'r--')
plot([t(1), t(end)], mu_lim(:,2), 'b--')
plot([t(1), t(end)], mu_lim(:,3), 'g--')
xlabel('Time, hours')
ylabel('Apparent growth rate')
legend('\Omega = 1.1','\Omega = 1.5','\Omega = 2')
axis([t(1), t(end), 0, 1.1*eta])