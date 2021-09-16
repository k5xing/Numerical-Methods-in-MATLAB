syms T;
syms f1(T);
syms f2(T);
%take resistance
R = input('Enter the resistance value of R:');

f1= 5.775*(10^(-7))*T^2-3.9083*(10^(-3))*T+(R/100-1);
f2=5.775*(10^(-7))*T^2-3.9083*(10^(-3))*T+4.183*(10^(-12))*(T-100)*T^3+R/100-1;

%check which function to use depends on R
if R<100
    newrap(f2, -100,1);
    bisection(f2,-200,0,1);
else 
     newrap(f1, 300,1);
     bisection(f1,0,850,1);
end

%bisection function
function bisection(f,init_1,init_2,counter)
    value_1=subs(f,init_1);
    value_2=subs(f,init_2);
    
    %update prev
    prev=init_1;
    
    %new root
    result=(init_1+init_2)/2;
    
    %get f(T) using new root
    value_new=subs(f,result);
    
    %recursively call biesection function until relative error < 0.1% 
    if(abs((result-prev)/result)<=0.001)
        disp(['Bisection: ',num2str(double(result))]);
         disp(['Iterations using Bisection: ', num2str(counter)]);
    %update the range of the root
    elseif value_new*value_1<0
        %keep tracking of iterations
        counter = counter + 1;
        bisection(f,init_1,result,counter);
    elseif value_new*value_2<0
        counter = counter + 1;
        bisection(f,result,init_2,counter);
    end
end

%Newton Raphson function
function newrap(f, init, counter)
    %update prev
    prev = init;
    %differentiate the function
    df = diff(f);
    %new root
    new_T = double(prev - subs(f, prev)/subs(df, prev));
    
    %recursively call Newton Raphson function until relative error < 0.1% 
    if(abs(new_T - prev) <= 0.001)
        disp(['Newton Raphson: ' ,num2str(double(new_T))]);
        disp(['Iterations using Newton Raphson: ', num2str(counter)]);
    else
        %keep tracking of iterations
        counter = counter + 1;
        newrap(f, new_T,counter);
    end
    
    
end
    
    



