function [uind,L]=BiotSavart(F1,F2,P,gamma,rc,d,co,type)
%% uind=BiotSavart(F1,F2,P,gamma,rc,d,type) -> Biot-Savart Law
%
% Function computes the velocity contributions due to turbine motion and
% freestream flow in the inertial and blade coordinate systems.
% 
% ****Input(s)****
% F1        Array containing first point of each vortex filament
% F2        Array containing second point of each vortex filament
% P         Array containing points of interest (where induction is
%           computed)
% gamma     Array of vortex filament circulation strengths
% rc        Vortex core sizes (actually radius squared for code speed-up)
% d         Squared cut-off distance (if =0, then viscous correction used)
% co        Distance from wake nodes beyond which influence is negligible
% type      If 'length', then will only output filament length (for
%           filament stretching correction), if 'full', will compute
%           induction on all points of interest
% 
% ****Output(s)****
% uind      Array of induced velocity at each of the points P due to
%           contributions from filaments defined by F1 and F2
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

%% Relabel filament endpoint variables, preallocate memory
sp=size(P); %Size of 4D array containing induced velocity points
if length(sp)==2
    sp(3)=1;
end
if length(sp)==3
    sp(4)=1;
end
ns=sp(1);
nt=sp(3);
nb=sp(4);

uind=zeros(sp);

if strfind(d,'visc')
    n=str2double(d(5:end));
else
    n=0; 
end

gammaout = size(gamma);

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

if strcmp(type,'length') %If true, then only returns filament length
    L=sqrt(L);
    uind=zeros(size(P));
elseif strcmp(type,'full')

%% Begin looping over POIs
    for k=1:nb
        for j=1:nt
            for i=1:ns
                px=P(i,1,j,k);
                py=P(i,2,j,k);
                pz=P(i,3,j,k);                
                
                
%% Compute vector difference calculations
                pxx1=px-x1;
                pyy1=py-y1;
                pzz1=pz-z1;
                pxx2=px-x2;
                pyy2=py-y2;
                pzz2=pz-z2;
    
%% Compute distances between points on triangle (filament to POI)
                r1=sqrt(pxx1.^2+pyy1.^2+pzz1.^2);
                r2=sqrt(pxx2.^2+pyy2.^2+pzz2.^2);            
                r1dr2=pxx1.*pxx2+pyy1.*pyy2+pzz1.*pzz2;
                r1tr2=r1.*r2;
                
                if (n~=0)
                    Ldr12=(x2x1.*pxx1+y2y1.*pyy1+z2z1.*pzz1).^2;
                    Cnu=((r1.^2)-(Ldr12./L));
                    Cnu=Cnu.*(rc.^n+Cnu.^n).^(-1/n);
                    ubar=Cnu.*gamma/(4*pi).*(r1+r2)./(r1tr2.*(r1tr2+r1dr2));
                else
                    ubar=(gamma/(4*pi)).*(r1+r2)./(r1tr2.*(r1tr2+r1dr2)+(d*L));
                end
                
                ubar(isnan(ubar) | isinf(ubar) | (r1>co & r2>co))=0;
                
                uind(i,1,j,k)=sum(sum(sum(ubar.*(pyy1.*pzz2-pzz1.*pyy2),1),3),4);
                uind(i,2,j,k)=sum(sum(sum(ubar.*(pzz1.*pxx2-pxx1.*pzz2),1),3),4);
                uind(i,3,j,k)=sum(sum(sum(ubar.*(pxx1.*pyy2-pyy1.*pxx2),1),3),4);               
            end
            %display('Went through the J loop');
        end
    end
end