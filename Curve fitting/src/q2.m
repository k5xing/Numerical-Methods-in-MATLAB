%Loading the test file to X
X = load('test1.txt');
%Get R^2 of each method
[la0, la1,Ylinear,linearrsqr]=linearmethod(X);
[polya, Ypoly,polyrsqr,degree]=polynomial(X);
[expa, expb, Xexp,Yexp,exprsqr]=exponential(X);
[powa, powb, Xpow,Ypow,powrsqr]=powerMethod(X);
rsqr = min([abs(1-linearrsqr),abs(1-polyrsqr),abs(1-exprsqr),abs(1-powrsqr)]);

%Determine which function is the best fit
if(rsqr == 1-linearrsqr)
    X1 = X(:,1);
    X2= X(:,2);
    plot(X1, X2, X1, Ylinear);
    grid on;
    xlabel('x'); ylabel('y');
    xt = max(X1)*0.1;
    yt = max(X2)*0.8;
    caption = sprintf('y = %f * x + %f \nR^2 = %f', la1, la0,linearrsqr);
    text(xt,yt, caption, 'FontSize', 12, 'Color', 'black', 'FontWeight', 'bold');

   disp('Linear fit is the best');

elseif(rsqr == 1-polyrsqr)
    display(['Degree: ', num2str(degree)]);
    k = 1;
    while k <= degree+1
        display(['a',num2str(k-1),': ', num2str(polya(k,1))]);
        k = k+1;
    end
    plot(X(:,1), X(:,2),  X(:,1), Ypoly);
    grid on;
    xlabel('x'); ylabel('y');
    xt = max(X(:,1))*0.1;
    yt = max(X(:,2))*0.8;
    caption = sprintf('R^2 = %f', polyrsqr);
    text(xt,yt, caption, 'FontSize', 12, 'Color', 'black', 'FontWeight', 'bold');
    disp('Polynomial fit is the best');


elseif(rsqr == 1-exprsqr)
    
    plot(X(:,1), Xexp(:,2), X(:,1), Yexp);
    grid on;
    xlabel('x'); ylabel('y');
    xt = max(X(:,1))*0.1;
    yt = max(Xexp(:,2))*0.8;
    caption = sprintf('y = %f * e^{%f*x} \nR^2 = %f', expa, expb,exprsqr);
    text(xt,yt, caption, 'FontSize', 12, 'Color', 'black', 'FontWeight', 'bold');
    disp('Exponential fit is the best');

elseif(rsqr == 1- powrsqr)
    
    plot(Xpow(:,1), Xpow(:,2), Xpow(:,1), Ypow);
    grid on;
    xlabel('x'); ylabel('y');
    xt = max(Xpow(:,1))*0.1;
    yt = max(Xpow(:,2))*0.8;
    caption = sprintf('y = %f * x^{%f} \nR^2 = %f', powa, powb,powrsqr);
    text(xt,yt, caption, 'FontSize', 12, 'Color', 'black', 'FontWeight', 'bold');
    disp('Power fit is the best');

end

%Linear_helper
function [a0,a1,X1,X2,X1bar,X2bar]=linear(X)
    X1 = X(:,1);
    X2= X(:,2);
    Sx= sum(X1.*X2);
    Sx1= sum(X1);
    Sx2=sum(X2);
    Sx1sqr=sum(X1.*X1);
    a1=(length(X)*Sx-Sx1*Sx2)/(length(X)*Sx1sqr-Sx1^2);
    X1bar=Sx1/length(X);
    X2bar=Sx2/length(X);
    a0=X2bar-a1*X1bar;
   
   
end

%Linear fit
function [a0, a1, Ylinear,linearrsqr]=linearmethod(X)
   [a0,a1,X1,X2,X1bar,X2bar]=linear(X);
    Ylinear=a0+a1*X1;
    Stprep=X2 - X2bar;
    St=sum(Stprep.*Stprep);
    Srprep=X2-Ylinear;
    Sr=sum(Srprep.*Srprep);
    linearrsqr=(St-Sr)/St;
end

%Polynomial fit
function [xa,Ypoly,polyrsqr, degree] =polynomial(X)
    polyrsqr=2;
    prevrsqr=1;
    degree=0;
    while(abs((polyrsqr-prevrsqr)/polyrsqr)>0.01)
        degree=degree+1;
        prevrsqr=polyrsqr;
        xpower = zeros(length(X), degree+1);
        j = 1;

        %Matrix which contains the powers of x
        while j <= degree*2
            i = 1;
            while i <= length(X)
                xpower(i,j) = power(X(i,1), j);
                i = i+1;
            end
            j = j+1;
        end

        xleft = zeros(degree+1, degree+1);
        offset = -1;
        column = 1;

        %Matrix on the left side of coefficients
        while column <= degree+1
            row = 1;
            %set up the first column
            if row == 1 && column == 1
                xleft(row,column) = length(X);
                row = row + 1;

                while row <= degree + 1
                      xleft(row, column) = sum(xpower(:,row-1));
                      row = row + 1;
                end
            %set up the rest of the columns
            else
                while row <= degree + 1
                    xleft(row,column) = sum(xpower(:,row+offset));
                    row = row + 1; 
                end
            end
            column = column + 1; 
            offset = offset + 1;
        end

        %Matrix on the right side of coefficients
        xright = zeros(degree+1, 1);
        row = 1;
        while row <= degree+1
            %xright(row,1) = sum(xpower(:,row).*X(:, 2)); 
            xright(row,1) = sum(power(X(:,1),row-1).*X(:, 2)); 
            row = row+1;
        end

        %Coefficient matrix
        xa =xleft\xright;
        Ypoly = 0;
        d = 1;
        while d <= degree+1
            Ypoly = Ypoly + xa(d,1)*power(X(:,1),d-1);
            d = d +1;
        end
    
    
    %Calculate R^2
    Ybar=sum(X(:,2))/length(X);
    Stprep= Ypoly - Ybar;
    St=sum(Stprep.*Stprep);
    Srprep=X(:,2)-Ypoly;
    Sr=sum(Srprep.*Srprep);
    polyrsqr=(St-Sr)/St;
    
    end
end

%Exponential fit
function [a,b, Xnew,Yexp,exprsqr]=exponential(X)
    Xnew =X(X(:,2)>0,:);  
        
    Xln=[Xnew(:,1) log(Xnew(:,2))];
    [a0,a1,X1,X2,X1bar,X2bar]=linear(Xln);
    a=exp(a0);
    b=a1;
    Yexp=a*exp(b*X1);
    Xnewbar=sum(Xnew(:,2))/length(Xnew);
    Stprep=Xnew(:,2) - Xnewbar;
    St=sum(Stprep.*Stprep);
    Srprep=Xnew(:,2)-Yexp;
    Sr=sum(Srprep.*Srprep);
    exprsqr=(St-Sr)/St;
end

%Power fit
function [a, b, Xnew,Ypow,powrsqr]= powerMethod(X)
    Xnew =X(X(:,1)>0 & X(:,2)>0, :); 
    Xlog=[log10(Xnew(:,1)) log10(Xnew(:,2))];
    [a0,a1,X1,X2,X1bar,X2bar]=linear(Xlog);
    a=power(10,a0);
    b=a1;
    X1new=Xnew(:,1);
    Ypow=a*(power(X1new,b));
    Xnewbar=sum(Xnew(:,2))/length(Xnew);
    Stprep=Xnew(:,2) - Xnewbar;
    St=sum(Stprep.*Stprep);
    Srprep=Xnew(:,2)-Ypow;
    Sr=sum(Srprep.*Srprep);
    powrsqr=(St-Sr)/St;
end