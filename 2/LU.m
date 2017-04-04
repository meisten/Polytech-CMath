function LU

% {block initial conditions}
    n=6; % ?????? ??????? ?
    a=randn(n); % ??????? ?
    m=eye(n); % ????????? ??????? ???????? n
    aCopy=a;
    p1=eye(n); % ????????? ??? "P" 
% {endblock initial conditions}

% {block main}
    for k=1:n
        [~, w]=max(abs(a(k:n, k))); 
        w=w+k-1;
        v(k)=w; 
        e=eye(n);   
        p=e;
        p(k,:)=e(w,:);
        p(w,:)=e(k,:);
        a=p*a;
        L=eye(n);
        
        for j=k+1:n
            g=a(j,k)/a(k,k);
            L(j,k)=-g;
            m(j,k)=g;
        end
        
        p1=p*p1;
        a=L*a;
    end
    for j=1:n-2 
        for k=j+1:n-1
            if v(k)~= k
                d=m(k,j);
                m(k,j)=m(v(k),j);
                m(v(k),j)=d;
            end
        end
    end
    u=a;
    L=m;
    p=p1;
    PA=p*aCopy;
    LU=L*u;
% {endblock main}

% {block calculation x,y b}
    b=aCopy*ones(n, 1); % ?????? b
    b1=p*b;
    y=zeros(n,1);

    for i=1:n
        s=0;
        for k=1:i-1
            s=s+L(i,k)*y(k);
        end
        y(i)=b1(i)-s;
    end

    x=zeros(n,1);

    for i=n:-1:1
        s=0;
        for k=i+1:n
            s=s+u(i,k)*x(k);
        end
        x(i)=(y(i)-s)/u(i,i);
    end
% {endblock calculation x,y b}

% {block verification}
    display(PA)
    display(LU)
    display(x)
    display(y)
% {endblock verification}
display(b)