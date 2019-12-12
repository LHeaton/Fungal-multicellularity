
t_end = 20*24;
Omega = 2;

res = 100*t_end;
eta = [1, 0.5, 0.2, 0.1];
v = 1;
T = Omega./eta;

V0 = 10^(-7);

f = zeros(res, 4);
mu = zeros(res, 4);

time = linspace(0, t_end, res);
dt = time(2) - time(1);

for j = 1:4
    
    tt = find(time >= T(j), 1, 'first');
    T(j) = time(tt);

    V = zeros(res, 1);
    r = zeros(res, 1);
    V(1) = V0;
    r(1) = (3*V0/(4*pi))^(1/3);
    
    mu(1,j) = eta(j);
    
    for i = 2:res
        
        V(i) = V(i-1)*(1 + eta(j)*(1 - f(i-1,j))*dt);
            
        r(i) = (3*V(i)/(4*pi))^(1/3);
        
        if r(i) > r(i-1) + v*dt
            
            r(i) = r(i-1) + v*dt;
            V(i) = (4/3)*pi*r(i)^3;
        end
        
        if i > tt
            
            f(i,j) = V(i - tt)/V(i);
            
        end
        
        mu(i,j) = (V(i) - V(i-1))/(dt*V(i));
    end
end

figure(3)
plot(time/24, f(:,1), 'r')
hold on
plot(time/24, f(:,2), 'b')
plot(time/24, f(:,3), 'g')
plot(time/24, f(:,4), 'm')
xlabel('Time, days')
ylabel('Fraction with exhausted resource')
axis([0 t_end/24 0 1])

figure(4)
plot(time/24, mu(:,1), 'r')
hold on
plot(time/24, mu(:,2), 'b')
plot(time/24, mu(:,3), 'g')
plot(time/24, mu(:,4), 'm')
xlabel('Time, days')
ylabel('Apparent growth rate')
