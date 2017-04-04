function dy = rigid(t,y)
dy = zeros(3,1);  
dy(1) = y(2);
dy(2) = y(3);
dy(3) =(-1/4)*y(1)-(3/4)*y(2)-(27/4)*y(3)+10/4;
