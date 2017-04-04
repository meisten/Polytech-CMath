function RungeKutta
% block (initial conditions)
    h = 0.1; 
    Dt = 80; 
    n = 3;
    A = [
        0 1 0; 
        0 0 1; 
        -1/4 -3/4 -27/4;
    ];
    eig(A)
    B = [0; 0; 10/4]; 
    Y = [0; 0; 0]; 
    t = 0; 
    m = floor(Dt/h);    
% endblock (initial conditions)

% block (Runge Kutta method solution 4th degree)
    RK.solution4th = zeros(m+1,n+1);    
    RK.solution4thHalf = zeros(m + 1,n + 1);
    
    RK.solution4th(1,:) = [t,Y'];
    for i =2:m + 1
        [t,Y] = RungeQt(t, h, Y, A, 1, B);
        RK.solution4th(i,:)=[t,Y']; 
    end
% endblock (Runge-Kutta method solution)

% block (Runge-Kutta-Half method solution)    
    Y=[0;0;0]; 
    t=0;
    RK.solution4thHalf(1,:) = [t,Y']; 
    for i =2:m + 1
        [t,Y] = RungeQt(t, h, Y, A, 2, B);
        RK.solution4thHalf(i,:)=[t,Y']; 
    end 
% endblock (Runge-Kutta-Merson solution with a half step) 

% block (Runge Kutta method solution dsolve)    
    RK.dsolve = subs(vpa(dsolve('4*D3y+27*D2y+3*Dy+y=10','y(0)=0','Dy(0)=0','D2y(0)=0')),RK.solution4th(:,1)); 
% endblock (Runge Kutta method solution dsolve)   

% block (output)
    subplot(3, 1, 1); 
    plot(RK.solution4th(:,1),RK.solution4th(:,2));
    grid; 
    xlim([0 Dt])
    title 'RK Method: Solution';

    subplot(3, 1, 2);
    plot(RK.solution4th(:,1), (2^4)*(RK.solution4th(:,2) - RK.solution4thHalf(:,2))/(2^4-1));  
    grid;
    xlim([0 Dt])
    title 'RK Method: Error estimate Runge-Kutta method';

    subplot(3, 1, 3)
    plot(RK.solution4th(:,1), RK.solution4th(:,2) - RK.dsolve);
    grid;
    xlim([0 Dt])
    title 'RK Method: Reference error';
% endblock (output)

% block (RungeQt function)
function [t,Y] = RungeQt(t, H, Y, A, k, B)
    h = H/k;
    for i=1:k
        K1 = RightPart(t,Y,A,B);
        K2 = RightPart(t + (h/2),Y + (h/2)*K1, A, B);
        K3 = RightPart(t + (h/2),Y + (h/2)*K2, A, B);
        K4 = RightPart(t + (h),Y + (h)*K3, A, B);

        Y = Y + (h/6)*(K1+2*K2+2*K3+K4);      
        t = t + h;
    end
% endblock (RungeQt function)

