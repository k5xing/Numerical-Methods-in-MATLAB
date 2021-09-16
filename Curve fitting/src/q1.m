%Loading the test file to X
X = load('test3.txt');
method=input('Select the function to fit the data: ');

%Linear fit
if(method==1)
    %Call linear function
    [a0,a1,X1,X2,X1bar,X2bar]=linear(X);
    %Calculate Y using coeffs found
    Y=a0+a1*X1;
    %Calculate R^2
    Stprep=X2 - X2bar;
    St=sum(Stprep.*Stprep);
    Srprep=X2-Y;
    Sr=sum(Srprep.*Srprep);
    rsqr=(St-Sr)/St;
    %Plot the graph and show equation on the graph
    plot(X1, X2, X1, Y);
    grid on;
    xlabel('x'); ylabel('y');
    xt = max(X1)*0.1;
    yt = max(X2)*0.8;
    caption = sprintf('y = %f * x + %f \nR^2 = %f', a1, a0,rsqr);
    text(xt,yt, caption, 'FontSize', 12, 'Color', 'black', 'FontWeight', 'bold');

%Polynomial fit
elseif(method==2)
    degree=input('Determine the degree of the polynomial: ');
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
    Y = 0;
    d = 1;
    %Calculate Y using coeffs found
    while d <= degree+1
        Y = Y + xa(d,1)*power(X(:,1),d-1);
        d = d +1;
    end
    
    %Calculate R^2
    Ybar=sum(X(:,2))/length(X);
    Stprep= Y - Ybar;
    St=sum(Stprep.*Stprep);
    Srprep=X(:,2)-Y;
    Sr=sum(Srprep.*Srprep);
    rsqr=(St-Sr)/St;
    
    %Display degree and coeffs
    display(['Degree: ', num2str(degree)]);
    k = 1;
    while k <= degree+1
        display(['a',num2str(k-1),': ', num2str(xa(k,1))]);
        k = k+1;
    end
    
    %Plot the graph and show equation on the graph
    plot(X(:,1), X(:,2),  X(:,1), Y);
    grid on;
    xlabel('x'); ylabel('y');
    xt = max(X(:,1))*0.1;
    yt = max(X(:,2))*0.8;
    caption = sprintf('R^2 = %f', rsqr);
    text(xt,yt, caption, 'FontSize', 12, 'Color', 'black', 'FontWeight', 'bold');
    
elseif(method==3)
    %Filter out negative values
    Xnew =X(X(:,2)>0,:);
    
    %If zero or negative values are found
    if(length(Xnew)~=length(X))
        disp('Warning: expect positive values, ignoring non-positive values');
    end
    
    %Take ln of y
    Xln=[Xnew(:,1) log(Xnew(:,2))];
    [a0,a1,X1,X2,X1bar,X2bar]=linear(Xln);
    a=exp(a0);
    b=a1;
    
    %Calculate Y using coeffs found
    Y=a*exp(b*X1);
    
    %Calculate R^2
    Xnewbar=sum(Xnew(:,2))/length(Xnew);
    Stprep=Xnew(:,2) - Xnewbar;
    St=sum(Stprep.*Stprep);
    Srprep=Xnew(:,2)-Y;
    Sr=sum(Srprep.*Srprep);
    rsqr=(St-Sr)/St;
    
    %Plot the graph and show equation on the graph
    plot(X1, Xnew(:,2), X1, Y);
    grid on;
    xlabel('x'); ylabel('y');
    xt = max(X1)*0.1;
    yt = max(Xnew(:,2))*0.8;
    caption = sprintf('y = %f * e^{%f*x} \nR^2 = %f', a, b,rsqr);
    text(xt,yt, caption, 'FontSize', 12, 'Color', 'black', 'FontWeight', 'bold');

elseif(method==4)
    %Filter out negative values
    Xnew =X(X(:,1)>0 & X(:,2)>0, :); 
    
    %If zero or negative values are found
    if(length(Xnew)~=length(X))
        disp('Warning: expect positive values, ignoring non-positive values');
    end
    
    %Take log of x,y
    Xlog=[log10(Xnew(:,1)) log10(Xnew(:,2))];
    
    %Call linear function
    [a0,a1,X1,X2,X1bar,X2bar]=linear(Xlog);
    a=power(10,a0);
    b=a1;
    X1new=Xnew(:,1);
    
    %Calculate Y using coeffs found
    Y=a*(power(X1new,b));
    
    %Calculate R^2
    Xnewbar=sum(Xnew(:,2))/length(Xnew);
    Stprep=Xnew(:,2) - Xnewbar;
    St=sum(Stprep.*Stprep);
    Srprep=Xnew(:,2)-Y;
    Sr=sum(Srprep.*Srprep);
    rsqr=(St-Sr)/St;
    
    %Plot the graph and show equation on the graph
    plot(X1new, Xnew(:,2), X1new, Y);
    grid on;
    xlabel('x'); ylabel('y');
    xt = max(X1)*0.1;
    yt = max(Xnew(:,2))*0.8;
    caption = sprintf('y = %f * x^{%f} \nR^2 = %f', a, b,rsqr);
    text(xt,yt, caption, 'FontSize', 12, 'Color', 'black', 'FontWeight', 'bold');

end

%Linear fit method
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