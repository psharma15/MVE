%% MVE

function [B,d] = mve(A,b)

%  MVE finds the maximum volume ellipsoid bounded by polyhedron Ax <= b
    % x0 is the initial guess
    % Using MATLAB Optimization toolbox.
    % This is constrained nonlinear multivariable problem.
    % The ellipsoid is affine transform of a unit circle e = {B*u+d | ||u|| <=1}.
    % -----------------------------------------------------------------
    % Pragya Sharma, Cornell University, 08-25-2016
    % -----------------------------------------------------------------

[~,p] = size(A);

fun = @(x)log(det(inv(x(1:p,1:p))));

x0 = [ eye(p,p) zeros(p,1) ];
lb = ones(p,p+1) * -Inf;
ub = ones(p,p+1) * Inf;
     
A1=[];
b1=[];
Aeq=[];
beq=[];

x=fmincon(fun,x0,A1,b1,Aeq,beq,lb,ub,@(x)mveconstr(x,A,b,p));
B=x(1:p,1:p);
d=x(:,p+1);
end
