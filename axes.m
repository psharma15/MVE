%% AXES
function [a,b,varargout]=get_axes(B,d)

% AXES finds the axes for 2-D (ellipse) or 3-D (ellipsoid) input
  % B and d are the output of mve function
  % ----------------------------------------------------------------
  % Pragya Sharma, Cornell University, 08-25-2016
  % ----------------------------------------------------------------
  
% Center form: (x-d)'A(x-d) <= 1
invB=inv(B);
transinvB=ctranspose(invB);
A=transinvB*invB;

% Extraxting axes by SVD
[~,D,V]=svd(A);


a=1/sqrt(D(1,1));
b=1/sqrt(D(2,2));
if length(d) == 3
  varargout=1/sqrt(D(3,3));
end
