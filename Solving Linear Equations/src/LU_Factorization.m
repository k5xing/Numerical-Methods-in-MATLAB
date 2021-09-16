%Loading 'A.txt' and 'B.txt'
A3 = load('A.txt');
B3 = load('B.txt');
%Find the lower & higher triangular matrix
[L,U] = lu(A3);
C = L\B3;
D = U\C;
%Keep 3 sig digits
D =  round(D, 3, 'significant');
%Call function to display the calculated currents
print(D, length(B3));

function print(I, n)
    numCurrent = 1;
    disp('LU Factorization: ');
    while  numCurrent <= n
        disp(['I', num2str(numCurrent), ' = ', num2str(I(numCurrent, 1))]);
        numCurrent = numCurrent + 1;
    end
end