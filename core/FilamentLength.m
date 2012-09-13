function [L]=FilamentLength(F1,F2)
%% L=FilamentLength(F1,F2,P,gamma,rc,d) -> Biot-Savart Law
%
% Function computes the length of vortex filaments
% 
% ****Input(s)****
% F1        Array containing first point of each vortex filament
% F2        Array containing second point of each vortex filament
% ****Output(s)****
% L         Filament length
% 
%
% This work is licensed under the Open Source Initiative BSD 2-Clause 
% License. To view a copy of this license, visit 
% http://www.opensource.org/licenses/BSD-2-Clause
% 
%
% Written by Thomas Sebastian (tommy.sebastian@gmail.com)
% Last edited May 24, 2011
%

%% Calculate Filament Length

%Filament start points
x1=F1(:,1,:,:);
y1=F1(:,2,:,:);
z1=F1(:,3,:,:);
clear F1

%Filament end points
x2=F2(:,1,:,:);
y2=F2(:,2,:,:);
z2=F2(:,3,:,:);
clear F2

x2x1=x2-x1;
y2y1=y2-y1;
z2z1=z2-z1;
L=x2x1.^2+y2y1.^2+z2z1.^2; %Length of vortex filament (NOTE: L is L^2, as rc is rc^2)

L=sqrt(L);

size(L);




