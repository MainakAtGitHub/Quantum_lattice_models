libfound=false;
eigeninclude_default='/usr/include/eigen3/';
while ~libfound
    eigeninclude=input('Directory to sources of Eigen [default: /usr/include/eigen3/] :','s');
    if isempty(eigeninclude)
        eigeninclude=eigeninclude_default;
    end;
    d=dir(eigeninclude);
    for i=1:numel(d)
        if isequal(d(i).name,'Eigen')
            libfound=true;
        end;
    end;
    if ~libfound
        disp(['Libraries not foun in', eigeninclude])
    end;
end
mex('SelfAdjointEigenSolver.cpp',['-I',eigeninclude]);