function [Ai] = SurfaceCubicInterpolator(A,n)
%A : matrix to interpolate
%n : points to interpolate over

A_size = size(A);
Ai = ones( A_size(1),n);

for i=1:A_size(1)
    Ai(i,:) = interp1(linspace(0,1,A_size(2)), A(i,:), linspace(0,1,n), 'PCHIP');
end