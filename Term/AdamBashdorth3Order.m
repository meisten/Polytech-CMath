function AdamBashdorth3Order
h = 0.01;
DT = 80;
n = 3; 
a = [
    0 1 0; 
    0 0 1; 
    -1/4 -3/4 -27/4
];
y = [0; 0; 0]; 
m = floor (DT/h);
t = 0;
f = [0; 0; 10/4]; 
eig (a)
 
out = zeros (m + 1, n + 1);
s=3;
f1=zeros(n,s);
out (1, :) = [t, y'];
f1(:,1)=mf(t,y,a,f);
 
for i = 2: s
   k1=mf(t,y,a,f);
   k2=mf(t+h/2, y+h*k1/2, a,f);
   k3=mf(t+h, y-k1*h+2*k2*h, a,f);
   y=y+h*(k1+4*k2+2*k3)/6;
   t = t + h;
   out(i,:)=[t,y'];
   f1(:,i)=mf(t,y,a,f);
end
 for i=s+1:m+1
     y=y+(h/12)*(23*f1(:,s)-16*f1(:,s-1)+5*f1(:,s-2));
     t=t+h;
     out (i, :) = [t, y'];
     for j=1:s-1
         f1(:,j)=f1(:,j+1);
     end
     f1(:,s)=mf(t,y,a,f);
 end
 
 x=out(:,1); 
 z=out(:,2);
 s1 = dsolve ('4*D3y+27*D2y+3*Dy+y=10', 'D2y(0) = 0', 'Dy(0) = 0', 'y(0) = 0');
 s1 = subs (s1, x);
 
 %построение графиков
 subplot(2,1,1);
 plot (x,z);
 axis tight;
 grid on;
 title ('Решение')
 
 subplot(2,1,2);
 plot(x, z-s1);
 axis tight;
 grid;
 title ('Погрешность')
