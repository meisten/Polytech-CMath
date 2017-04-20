function AdamBashforth_Ode
options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5]);
[t,y] = ode45(@AdamBashforth_Koshi,[0 80],[0 0 0],options);
s1 = subs(vpa(dsolve('4*D3y+27*D2y+3*Dy+y=10','y(0)=0','Dy(0)=0','D2y(0)=0')), t);

subplot(2,1,1);
plot(t,y(:,1),'-');
grid on;

subplot(2,1,2);
plot(t,y(:,1)-s1,'-');
grid on;
