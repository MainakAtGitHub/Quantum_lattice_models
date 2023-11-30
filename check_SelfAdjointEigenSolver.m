% check_SelfAdjointEigenSolver
% set up a random matrix to be diagonalized
A=rand(5,5)+1i*rand(5,5);
% make sure that it is self-adjoined
A=A+A';
% calculate Eigenvectors and eigenvalues
[V,E]=SelfAdjointEigenSolver(A);
% calculate the same with the built-in function
[Vi,Ei]=eig(A);
Ei=diag(Ei);
disp('difference in the eigenvalues: ')
delta=sum(abs(Ei-E))
% eigenvectors might differ by a phase
disp('difference in the eigenvectors (phase): ')
deltaev=V./Vi
disp('difference in the eigenvalues calculated by V''*A*V: ')
deltamatrix=sum(sum(abs(V'*A*V-Vi'*A*Vi)))

% same check for symmetric matrix
B=rand(5,5);
B=B+B';
[V,E]=SelfAdjointEigenSolver(B);
% calculate the same with the built-in function
[Vi,Ei]=eig(B);
Ei=diag(Ei);
disp('difference in the eigenvalues: ')
delta=sum(abs(Ei-E))
% eigenvectors might differ by a phase
disp('difference in the eigenvectors (phase): ')
deltaev=V./Vi
disp('difference in the eigenvalues calculated by V''*B*V: ')
deltamatrix=sum(sum(abs(V'*B*V-Vi'*B*Vi)))

% performance checks
n=input('dimension for performance check n= ');
if n> 10000
    n=10000;
end;
C=rand(n,n)+1i*rand(n,n);
C=C+C';
CV=C;
CE=C;
clear('CV','CE');
disp('complex matrix: diagonalization with SelfAdjointEigenSolver...')
tic;
[CV,CE]=SelfAdjointEigenSolver(C);
toc;
% clear result to avoid bias from memory usage
clear('CV','CE');
disp('Diagonalization with eig...')
tic;
[CV1,CE1]=eig(C);
toc;

% performance checks
C=rand(n,n);
C=C+C';
CV=C;
CE=C;
clear('CV','CE');
disp('real matrix: diagonalization with SelfAdjointEigenSolver...')
tic;
[CV,CE]=SelfAdjointEigenSolver(C);
toc;
% clear result to avoid bias from memory usage
clear('CV','CE');
disp('Diagonalization with eig...')
tic;
[CV,CE]=eig(C);
toc;