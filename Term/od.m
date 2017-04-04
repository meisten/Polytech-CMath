function od
options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5]);
[t,y] = ode45(@rigid,[0 80],[0 0 0],options);
x=t;
s1 = dsolve ('4*D3y+67.2*D2y+4.2*Dy+y=10', 'D2y(0) = 0', 'Dy(0) = 0', 'y(0) = 0');
s1 = subs (s1, x);
subplot(2,1,1);
plot(x,y(:,1),'-');
grid on;
subplot(2,1,2);
plot(x,y(:,1)-s1,'-');
grid on;
