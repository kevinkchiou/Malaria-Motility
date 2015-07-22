function out = circle_fit(positions)
%CIRCLE_FIT takes positional data in two-column or two-row format and finds the best-fit circle. The procedure is detailed in Materials and Methods of the primary manuscript
%	This subroutine uses a linear approach to fitting a circle to the input set of points. Circle fitting is non-linear when using the standard approach. However, with change of coordinates and measure it can be made linear, but not without cost. There is some sensitivity to error for tracks that subtend small amounts of solid angle or more than 2pi solid angle. These issues can generally be resolved by proper data window selection.

%change these warnings about rank deficiency into an error
w1=warning('error','MATLAB:rankDeficientMatrix');
w2=warning('error','MATLAB:singularMatrix');

sz=size(positions);
if sz(1)~=2 && sz(2) ~=2
    error('Error: position vector is incorrectly sized!\n');
end
%ensure column vector to do solution correctly
if sz(2)==2
    x=positions(:,1);y=positions(:,2);
    len=sz(1);
elseif sz(1)==2
    x=positions(1,:)';y=positions(2,:)';
    len=sz(2);
end
M = [2*x,2*y,ones(len,1)];
try
    resvec= M \ (x.^2 + y.^2);
    %resvec = [a;b;c^2] from the equation (x-a)^2 + (y-b)^2 = r^2 -->
    %2ax + 2by + c^2 = x^2 + y^2
catch err
    %fprintf('Catch: Returning NaN for radius parameters\n',err.message);
    out=[Inf;Inf;Inf];
    warning(w1);warning(w2);return;
end
x0=resvec(1);y0=resvec(2);
r=sqrt(resvec(3)+x0^2+y0^2);
out=[r;x0;y0];
warning(w1);warning(w2);
end

