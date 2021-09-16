%Loading 'A.txt' and 'B.txt'
A2 = load('A.txt');
B2 = load('B.txt');

%Append vector B2 to the last column of A2
A2 = [A2 B2];
%Perform row reduction on A2
X2 = rref(A2);
%Extract last column of the row-reduced matrix
Y2 = X2(:,end);
%Keep 3 sig digits
Y2 = round(Y2, 3, 'significant');
disp('Gaussian elimination: ');
print(Y2, length(B2));

%Append vector B2 to the last column of A3 with increased resistance
A3 = [1.05 * load('A.txt') B2];
%Perform row reduction on A3
X3 = rref(A3);
%Extract last column of the row-reduced matrix
Y3 = X3(:,end);
disp('After the values of all resistors increased by 5%: ');
%Keep 3 sig digits
Y3 = round(Y3, 3, 'significant');
print(Y3, length(B2));

%Initiate the Error vector
Diff= zeros(length(B2), 1);

i = 1;
%Calculate the difference between the original current and increased-resistance current
while i <= length(B2) 
    
    Diff(i,1) = 100*abs((Y3(i,1)-Y2(i,1))/Y3(i,1));
    i = i + 1;
    
end

disp('The percent change of currents after increasing the values of resistotrs by 5%: ')

%Keep 3 sig digits
Diff = round(Diff, 3, 'significant');

%display %change of current
numCurrent = 1;
    while  numCurrent <= n
        disp(['%change of I', num2str(numCurrent), ' = ', num2str(Diff(numCurrent, 1)), '%']);
        numCurrent = numCurrent + 1;
    end
    
disp('Conclusion: after calculating the percent change of currents, the difference is within 5%, which means the system of equations is not in ill-condition.')

%define print function
function print(I, n)
    numCurrent = 1;
    while  numCurrent <= n
        disp(['I', num2str(numCurrent), ' = ', num2str(I(numCurrent, 1))]);
        numCurrent = numCurrent + 1;
    end
end