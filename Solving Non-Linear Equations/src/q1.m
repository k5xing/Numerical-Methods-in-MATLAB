%define symbolic variables
syms x;
syms y;
syms var;
var=[x;y];
%define symbolic function
syms f(x,y);
%system of functions to be evaluted
%f=[2*x-y-exp(-x); -x+2*y-exp(-y)];
f=[4*x^2+y^2-13; x^2+y^2-10];
%initial value
init =[0.1;0.1];
%invoke the Newton Raphson function
newton(f,var,init);
function newton (f, var,init)
    %set up prev
    prev=init;
    %set up intial values for error 
    error=ones(length(init),1);
    %looping until satisfy error tolerance
    while max(error)>0.5
     
        %get jacobian
        jac = jacobian(f, var);
        %evaluate jacobian matrix
        j_value=double(subs(jac,var,prev));
        %evaluate f
        f_value=double(subs(f,var,prev)); 
        %update result
        result = prev - j_value\f_value;
        %update error
        error=100*abs((result-prev)./result);
        %update prev
        prev=result;
    end
     disp(result);
  
end
