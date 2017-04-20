function AdamBashdorth_Solve
% block (initial conditions)
    h = 0.01;
    Dt = 80;
    n = 3; 
    A = [
        0 1 0; 
        0 0 1; 
        -1/4 -3/4 -27/4
    ];
    Y = [0; 0; 0]; 
    m = floor (Dt/h);
    t = 0;
    B = [0; 0; 10/4]; 
    eig (A)
% endblock (initial conditions) 

% block (AdamsBashforth 3rd)
    % block (AdamsBashforth Solution)
        AB.solve = zeros (m + 1, n + 1);
        F = zeros(n,3);
        F(:,1) = RightPart(t,Y,A,B);
    	AB.solve(1, :) = [t, Y'];
        
        for i = 2:3
            k1 = RightPart(t,Y,A,B);
            k2 = RightPart(t+h/2, Y+h*k1/2, A,B);
            k3 = RightPart(t+h, Y-k1*h+2*k2*h, A,B);
            Y = Y + h*(k1+4*k2+k3)/6;
            t = t + h;  
            AB.solve(i,:) = [t,Y'];
            F(:,i) = RightPart(t,Y,A,B);
        end  
        
        for i = 4:m+1
            F(:,i) = RightPart(t,Y,A,B);
            Y = Y + (h/12)*(23*F(:,i) - 16*F(:,i-1) + 5*F(:,i-2));
            t = t + h;
            AB.solve(i,:) = [t,Y'];
        end       
    % endblock (AdamsBashforth Solution)
    
    % block (AdamsBashforth 3rd dsolve)
        AB.dsolve = subs(vpa(dsolve('4*D3y+27*D2y+3*Dy+y=10','y(0)=0','Dy(0)=0','D2y(0)=0')),AB.solve(:,1));
    % endblock (AdamsBashforth 3rd dsolve)
% endblock (AdamsBashforth 3rd)

% block (output)
     subplot(2,1,1);
     plot (AB.solve(:,1), AB.solve(:,2));
     axis tight;
     xlim([0 Dt])
     grid on;
     title ('Adams Bashdorth 3rd: Solution')  
     
     subplot(2,1,2);
     plot(AB.solve(:,1), AB.solve(:,2) - AB.dsolve);
     axis tight;
     xlim([0 Dt])
     grid;
     title ('Adams Bashdorth 3rd: Estimate')
 % endblock (output)
