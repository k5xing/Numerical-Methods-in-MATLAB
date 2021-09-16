%Loading 'A.txt' and 'B.txt'
A4 = load('A.txt');
B4 = load('B.txt');

n = length(B4);

%initial guess values
X_curr = ones(n,1);
X_prev = ones(n,1);
Error = ones(n,1);

%set output conditions
flag = 0;
i = 0;

while max(Error) >= 0.0001
    
    row = 1;
    %solve from X1 to Xn
    while row <= n
        col = 1;
        X_curr(row, 1) = B4(row, 1);
        
        %solve specific Xi
        while col <= n
          %skip diagonal element
          if col ~= row
             X_curr(row, 1)=  X_curr(row, 1)- A4(row, col) * X_curr(col, 1);
          end

          col=col+1;
        end
        X_curr(row,1) = X_curr(row, 1)/A4(row, row);
        
        %update absolute approximate relative error for each Xi
        Error(row,1) = abs((X_curr(row,1)- X_prev(row,1))/X_curr(row,1));
        %update previous Xi values, in order to calculate errors
        X_prev(row, 1) = X_curr(row, 1);

        row = row + 1;
    end
    i = i + 1;
    %change flag to avoid displaying the output of the same condition multiple times
    if max(Error)<0.01 && flag == 0
        display(['Required iterations for error less than 1%: ', num2str(i)]);
        I = round(X_curr, 3, 'significant');
        %Call function to display the calculated currents
        print(I,n);
        flag = 1;
    elseif max(Error)<0.005 && flag == 1
        display(['Required iterations for error less than 0.5%: ', num2str(i)]);
        I = round(X_curr, 3, 'significant');
        %Call function to display the calculated currents
        print(I,n);
        flag = 2;
    elseif max(Error)<0.001 && flag == 2
        display(['Required iterations for error less than 0.1%: ', num2str(i)]);
        I = round(X_curr, 3, 'significant');
        %Call function to display the calculated currents
        print(I,n);
        flag = 3;
    elseif max(Error) < 0.0001 && flag == 3
        display(['Required iterations for error less than 0.01%: ', num2str(i)]);
        I = round(X_curr, 3, 'significant');
        %Call function to display the calculated currents
        print(I,n);
        flag = 4;
    end
   
end

function print(I, n)
    numCurrent = 1;
        while  numCurrent <= n
            disp(['I', num2str(numCurrent), ' = ', num2str(I(numCurrent, 1))]);
            numCurrent = numCurrent + 1;
        end
end