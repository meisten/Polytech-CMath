function QR

% {block initial conditions}
    n = 3;
    A = rand(n);
    A1 = A;
    Q1 = eye (n);
% {endblock initial conditions}

% {block main}
    for k = 1:n-1
        for i = k+1:n
            if A(i,k)~=0
                d = sqrt (A(k,k)*A(k,k) + A(i,k)*A(i,k));
                C = A(k,k)/d;
                S = A(i,k)/d;
            end
            T = eye (n);
            T(k,k) = C;
            T(k,i) = S;
            T(i,i) = C;
            T(i,k) = -S;
            A = T*A;
            Q1 = T*Q1;
        end
    end
    R=A;
    disp(R)
    Q = Q1';
% {endblock main}

% {block verification}
    disp(Q*R)
    disp(A1)
    b = A*ones(n, 1);
    x = zeros(n,1);
    y=Q'*b;
% {endblock verification}

for i=n:-1:1
    S=0;
    for k=i+1:n
        S=S+R(i,k)*x(k);
    end
    x(i)=(y(i)-S)/R(i,i);
end