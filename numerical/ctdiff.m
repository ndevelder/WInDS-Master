function dy=ctdiff(x,y)
%% dy=spdiff(x,y) -> Centered differentiation.
%
% Function computes the derivative of the function 'y' wrt 'x' numerically
% via centered difference method
% 
% ****Input(s)****
% x         Vector of independent values
% y         Array of dependent values such that f(x)=y
% 
% ****Output(s)****
% dy        Array of derivative values of y wrt x
% 
%
% This work is licensed under the Open Source Initiative BSD 2-Clause 
% License. To view a copy of this license, visit 
% http://www.opensource.org/licenses/BSD-2-Clause
% 
%
% Written by Thomas Sebastian (tommy.sebastian@gmail.com)
% Last edited March 7, 2011
%

%% Compute estimate of derivative
dx=abs(mean(diff(x)));
dy=zeros(size(y));

if numel(size(y))==2
    dy(1,:)=(-y(3,:)+4*y(2,:)-3*y(1,:))/(2*dx); %Forward finite-difference, O(h^2)
    dy(2,:)=(y(3,:)-y(1,:))/(2*dx); %Centered finite-difference, O(h^2)    
    dy(3:end-2,:)=(-y(5:end,:)+8*y(4:end-1,:)-8*y(2:end-3,:)+y(1:end-4,:))/(12*dx); %Centered finite-difference, O(h^4)
    dy(end-1,:)=(y(end,:)-y(end-2,:))/(2*dx); %Centered finite-difference, O(h^2)
    dy(end,:)=(3*y(end,:)-4*y(end-1,:)+y(end-2,:))/(2*dx); %Backward finite-difference, O(h^2)
else
    dy(:,:,1,:)=(-y(:,:,3,:)+4*y(:,:,2,:)-3*y(:,:,1,:))/(2*dx); %Forward finite-difference, O(h^2)
    dy(:,:,2,:)=(y(:,:,3,:)-y(:,:,1,:))/(2*dx); %Centered finite-difference, O(h^2)
    dy(:,:,3:end-2,:)=(-y(:,:,5:end,:)+8*y(:,:,4:end-1,:)-8*y(:,:,2:end-3,:)+y(:,:,1:end-4,:))/(12*dx); %Centered finite-difference, O(h^4)
    dy(:,:,end-1,:)=(y(:,:,end,:)-y(:,:,end-2,:))/(2*dx); %Centered finite-difference, O(h^2)
    dy(:,:,end,:)=(3*y(:,:,end,:)-4*y(:,:,end-1,:)+y(:,:,end-2,:))/(2*dx); %Backward finite-difference, O(h^2)
end





