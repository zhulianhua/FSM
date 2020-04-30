function [ t, p ] = xyz_to_tp ( x, y, z )

%*****************************************************************************80
%
%% XYZ_TO_TP converts (X,Y,Z) to (Theta,Phi) coordinates on the unit sphere.
%
%  Modified:
%
%    09 September 2010
%
%  Author:
%
%    Dmitri Laikov
%
%  Reference:
%
%    Vyacheslav Lebedev, Dmitri Laikov,
%    A quadrature formula for the sphere of the 131st
%    algebraic order of accuracy,
%    Russian Academy of Sciences Doklady Mathematics,
%    Volume 59, Number 3, 1999, pages 477-481.
%
%  Parameters:
%
%    Input, real X, Y, Z, the Cartesian coordinates of a point
%    on the unit sphere.
%
%    Output, real T, P, the Theta and Phi coordinates of
%    the point.
%
  p = acos ( z );
  N = length(z);
  t = zeros(N,1);

  fact = sqrt ( x .* x + y .* y );

  for i = 1:N
      if (fact(i) > 0)
        t(i) = acos ( x(i) / fact(i));
      else
        t(i) = acos ( x(i));
        fprintf("t = %f, fact = %f\n", t(i), fact(i));
      end
    if ( y(i) < 0.0 )
        t(i) = - t(i);
    end
  end
  
  % Theata in [0, 2*pi];
  t = t + pi;
  
%
%  Convert to degrees.
%
%   t = t * 180.0 / pi;
%   p = p * 180.0 / pi;

  return
end
