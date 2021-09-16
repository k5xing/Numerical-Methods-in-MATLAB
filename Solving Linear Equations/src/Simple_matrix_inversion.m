%Loading 'A.txt' and 'B.txt'
A1 = load('A.txt');
B1 = load('B.txt');
%Compute inverse of A multiply by B
X1 = A1\B1; 
%Keep 3 sig digits
X1 = round(X1, 3, 'significant');
%Call function to display the calculated currents
disp('Simple matrix inversion: ');
print(X1, length(B1));


function print(I, n)
    numCurrent = 1;
    while  numCurrent <= n
        disp(['I', num2str(numCurrent), ' = ', num2str(I(numCurrent, 1))]);
        numCurrent = numCurrent + 1;
    end
end