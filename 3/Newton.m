function Newton

% {block initial conditions}
    disp('------- START ---------');
    x0 = [-656; -656];
    abs = 1e-5; 
    dx = x0;
    iterations = 0;
% {endblock initial conditions}

% {block main}
    while (norm(dx) > abs && (iterations < 1000))
        [F,Fx] = Function(x0);
        display(F);
        dx=-(Fx\F);
        iterations = iterations + 1;
        x0 = x0 + dx;
    end
% {endblock main}

% {block output}
    [F,~] = Function(x0);
    display(F);
    display(x0);
    display(iterations);
    if norm(dx) < abs
        disp('Solution had been found');
    else
        disp('Solution had not been found');
    end
% {endblock output}

% {block fsolve}
    fun = @Function;
    x0 = fsolve(fun,x0);
    display(x0);
% {enblock fsolve}

% {block solve}
    syms x y; 
    r = (solve('x^2*y^2 - 3*x^3 -6*y^3 + 8 = 0', 'x^4 - 9*y + 2 = 0'));
    display(vpa(r.x))
    % eqn2 = x^4 - 9*y + 2 == 0;
% {endblock solve}
    disp('------- END ---------');
    
function [Inition,Jacoby] = Function(x)
% f1 = x^2*y^2 - 3x^3 - 6y^3 + 8 = 0
% f2 = x^4 - 9y + 2 = 0
% 13 ???????

% (df1/dx, df1/dy) = (2xy^2 - 9x^2, 2x^2y - 18y^2)
% (df2/dx, df2/dy) = (4*x^3,        -9)

Inition=[
    (x(1)^2)*(x(2)^2) - 3*(x(1)^3) - 6*(x(2)^3) + 8; 
    (x(1)^4) - 9*x(2) + 2
    ];

Jacoby=[
    2*x(1)*(x(2)^2) - 9*(x(1))^2, 2*(x(1)^2)*x(2)-18*x(2)^2;
    4*(x(1)^3), -9
];


