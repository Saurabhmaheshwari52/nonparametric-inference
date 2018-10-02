function [c, ceq] = confun(x)
% x
ceq = [];
c = [1e6*(ones(1,(size(x,1)-8))*x(9)*x(9)-x(5)*x(6)); 
    1e6*(ones(1,(size(x,1)-8))*x(9)*x(9)-x(7)*x(8));
    abs(ones(1,(size(x,1)-8))*x(9)) - x(5);
    abs(ones(1,(size(x,1)-8))*x(9)) - x(6);
    abs(ones(1,(size(x,1)-8))*x(9)) - x(7);
    abs(ones(1,(size(x,1)-8))*x(9)) - x(8)];
% c
end