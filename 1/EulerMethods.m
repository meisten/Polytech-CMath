function EulerMethods
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

% block (Euler Forward Method)
    % block (Forward solution)
        EF.Solution = zeros(m + 1,n + 1);    
        EF.SolutionHalfStep = zeros(m + 1,n + 1);

        EF.Solution(1,:) = [t,Y'];

        for i = 2:m + 1
            [t, Y] = EulerForward(t, h, Y, A, 1, B);
            EF.Solution(i,:) = [t,Y']; 
        end
    % endblock (Forward solution)
 
    % block (Forward solution with a half step)
        Y=[0;0;0]; 
        t = 0;
        EF.SolutionHalfStep(1,:) = [t,Y'];
        for i = 2:m + 1
            [t, Y] = EulerForward(t, h, Y, A, 2, B);
            EF.SolutionHalfStep(i,:)=[t,Y'];
        end
    % endblock (Forward solution with a half step)

    % block (Forward solution dsolve)
        EF.dsolve = subs(vpa(dsolve('4*D3y+27*D2y+3*Dy+y=10','y(0)=0','Dy(0)=0','D2y(0)=0')),EF.Solution(:,1));
    % endblock (Forward solution dsolve)
% endblock (Euler Forward Method)

% block (Euler Backward Method)     
    % block (Backward solution)   
        Y=[0;0;0]; 
        t = 0;
        EB.Solution = zeros(m + 1,n + 1);    
        EB.SolutionHalfStep = zeros(m + 1,n + 1);

        EB.Solution(1,:) = [t,Y'];

        for i = 2:m + 1
            [t, Y] = EulerBackward(t, h, Y, A, 1, B);
            EB.Solution(i,:) = [t,Y']; 
        end      
    % endblock (Backward solution)
    
    % block (Backward solution with a half step)
        Y=[0;0;0]; 
        t = 0;
        EB.SolutionHalfStep(1,:) = [t,Y'];
        for i = 2:m + 1
            [t, Y] = EulerBackward(t, h, Y, A, 2, B);
            EB.SolutionHalfStep(i,:)=[t,Y'];
        end
    % endblock (Backward solution with a half step)
    
    % block (Forward solution dsolve)
        EB.dsolve = subs(vpa(dsolve('4*D3y+27*D2y+3*Dy+y=10','y(0)=0','Dy(0)=0','D2y(0)=0')),EB.Solution(:,1));
    % endblock (Forward solution dsolve)
% endblock (Euler Backward Method)

% block (output)
    figure('Position',[20 50 1100 600])
    % block (output Forward Method)
        subplot(3, 2, 1); 
        plot(EF.Solution(:,1),EF.Solution(:,2)); 
        grid; 
        title 'Forward Method: Solution';

        subplot(3, 2, 3)
        plot(EF.Solution(:,1),EF.Solution(:,2)-EF.SolutionHalfStep(:,2)); 
        grid; 
        title 'Forward Method: Error estimate';

        subplot(3, 2, 5)
        plot(EF.Solution(:,1),EF.SolutionHalfStep(:,2)-EF.dsolve);
        grid; 
        title 'Forward Method: Reference error';
    % endblock (output Forward Method)
    % block (output Backward Method)
        subplot(3, 2, 2); 
        plot(EB.Solution(:,1),EB.Solution(:,2)); 
        grid;
        xlim([0 Dt])
        title 'Backward Method: Solution';

        subplot(3, 2, 4)
        plot(EB.Solution(:,1),EB.Solution(:,2) - EB.SolutionHalfStep(:,2));
        xlim([0 Dt])
        grid; 
        title 'Backward Method: Error estimate';

        subplot(3, 2, 6)
        plot(EB.Solution(:,1),EB.SolutionHalfStep(:,2) - EB.dsolve);
        xlim([0 Dt])
        grid; 
        title 'Backward Method: Reference error';
    % endblock (output Backward Method)
%endblock (output)

% block (Euler Forward & Backward functions)
    % block (Euler Forward function) 
    function [t,Y] = EulerBackward(t, H, Y, A, k, B)
        h = H/k;
        for u = 1:k
            Y = (inv(eye(3) - h*A)) * Y + (inv(eye(3) - h*A)) * h * B;
            t = t + h;
        end
    % endblock (Euler Forward function)
    
    % block (Euler Backward function) 
    function [t,Y] = EulerForward(t, H, Y, A, k, B)
        h = H/k;
        for u = 1:k
            Y = Y + h*RightPart(t,Y,A,B); 
            t = t + h;
        end
    % endblock (Euler Backward function)
% endblock (Euler Forward & Backward functions)