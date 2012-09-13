function [cl,cd,phi,aoa,a,ap]=BEM(airfoils,blade,turbine,fastout,vel)
%% [cl,cd,phi,aoa,a,ap]=BEM(airfoils,blade,turbine,fastout,vel) -> BEM theory.
%
% Function computes spanwise and rotor performance and loads via blade
% element momentum theory. Includes corrections for skewed flow and 
% heavily loaded rotors.
% 
% ****Input(s)****
% airfoils  Structure containing airfoil performance tables
% blade     Structure containing blade geometry
% turbine   Structure containing turbine geometry
% fastout   Structure containing time-dependent kinematics
% vel       Structure containing velocity components in inertial and blade
%           coordinate systems
% 
% ****Output(s)****
% cl        Spanwise lift coefficient
% cd        Spanwise drag coefficient
% phi       Spanwise inflow angle
% aoa       Spanwise angle of attack
% a         Spanwise axial induction factor
% ap        Spanwise tangential induction factor
% 
%
% This work is licensed under the Open Source Initiative BSD 2-Clause 
% License. To view a copy of this license, visit 
% http://www.opensource.org/licenses/BSD-2-Clause
% 
%
% Written by Thomas Sebastian (tommy.sebastian@gmail.com)
% Last edited January 15, 2011
%

%% Preallocate space for variables within loop
%Determine size of test vectors/arrays
ns=length(blade.RNodes);
nt=length(fastout.Time); %Number of timesteps
na=length(airfoils.Names);
RNodes=blade.RNodes;
Chord=blade.Chord;
NFoil=blade.NFoil;

a0=zeros(ns,nt); %Old (previous iteration) axial induction factor
ap0=zeros(ns,nt); %Old (previous iteration) tangential induction factor
phi=zeros(ns,nt); %Local inflow angle
aoa=zeros(ns,nt); %Local angle of attack
cl=zeros(ns,nt); %Local lift coefficient
cd=zeros(ns,nt); %Local drag coefficient
ct=zeros(ns,nt); %Local thrust coefficient
ftip=zeros(ns,nt); %Tip loss factor
fhub=zeros(ns,nt); %Hub loss factor
f=zeros(ns,nt); %Total loss corection factor
fiter=zeros(ns,nt); %Converence flag for gridpoints ('1' if converged, '9999' if not)

%% Define convergence criteria
tol=1e-6; %Convergence tolerance
da=ones(ns,nt); %Set initial value for axial induction factor residual equal to 1
dap=ones(ns,nt); %Set initial value for tangential induction factor residual equal to 1

ncv=find(da>tol | dap>tol); %Identify all nonconverged points (all initially)
miter=5000; %Maximum number of allowable iterations
wt=0.1; %Weighting factor on corrections to balance speed with stability (faster as you 
        %approach 1, but less stable)

%% Compute relevant velocity/angle components
Uinf=sqrt(sum(vel.relhub.^2,2));
Om=fastout.RotSpeed*(2*pi/60);

twst=-blade.AeroTwst*pi/180;
ptch=fastout.BldPitch1*pi/180;

rP=-fastout.PtfmPitch*pi/180; %Rotor pitch (vector, wrt time)
rY=(fastout.PtfmYaw+fastout.NacYaw)*pi/180; %Rotor yaw (vector, wrt time)
if sign(rP)==0
    sp=sign(rY);
elseif sign(rY)==0
    sp=sign(rP);
else
    sp=sign(rP).*sign(rY);
end
gamma=sp.*acos(cos(rP).*cos(rY)); %Total skew angle
psi=pi-atan2(cos(rP).*sin(rY),sin(rP)); %Total azimuthal angle of skew

%% Compute initial guesses of key variables
Om=repmat(Om',ns,1);
Uinf=repmat(Uinf',ns,1);
ptch=repmat(ptch',ns,1);
gamma=repmat(gamma',ns,1);
psi=repmat(psi',ns,1);
twst=repmat(twst,1,nt);
sigmap=repmat(turbine.NumBl.*Chord./(2.*pi.*RNodes),1,nt); %Local solidity
RNodes=repmat(RNodes,1,nt);
lambdar=Om.*RNodes./Uinf; %Local speed ratio


% Initial values for axial and tangential induction factors
a=real(0.25*(2+pi*lambdar.*sigmap-sqrt(4-4*pi*lambdar.*sigmap+pi*lambdar.^2.*sigmap.* ...
    (8*(twst+ptch)+pi*sigmap))));
ap=zeros(size(a));

%% Primary loop for BEM

for j=1:200
    
    % Save previous values of axial and tangential induction factors
    a0(ncv)=a(ncv);
    ap0(ncv)=ap(ncv);
    
    % Compute inflow angle and angle of attack
    phi(ncv)=atan2(Uinf(ncv).*(1-a(ncv)),Om(ncv).*RNodes(ncv).*(1+ap(ncv)));    
    aoa(ncv)=(phi(ncv)-(twst(ncv)+ptch(ncv)))*180/pi;
    
    % Interpolate over airfoil database for lift and drag coefficients
    for k=1:na
        cl(NFoil==k,:)=interp1(airfoils.profiles(k,1).AoA,airfoils.profiles(k,1).Cl, ...
            aoa(NFoil==k,:));
        cd(NFoil==k,:)=interp1(airfoils.profiles(k,1).AoA,airfoils.profiles(k,1).Cd, ...
            aoa(NFoil==k,:));
    end
    
    % Compute elemental thrust coefficient
    ct(ncv)=sigmap(ncv).*(1-a(ncv)).^2.*(cl(ncv).*cos(phi(ncv))+cd(ncv).*sin(phi(ncv)))./ ...
        sin(phi(ncv)).^2;
    
    % Compute loss correction factor due to tip and hub losses
    ftip(ncv)=2./pi.*acos(exp(-(turbine.NumBl.*(blade.TipRad-RNodes(ncv))./ ...
        (2.*RNodes(ncv).*sin(phi(ncv)))))); %Tip loss factor
    fhub(ncv)=2./pi.*acos(exp(-(turbine.NumBl.*(RNodes(ncv)-blade.HubRad)./ ...
        (2*blade.HubRad.*sin(phi(ncv)))))); %Hub loss factor
    f(ncv)=fhub(ncv).*ftip(ncv); %Total loss correction factor
    
    % Compute axial induction factor using conventional BEM theory
    a(ncv)=real((1+4.*f(ncv).*sin(phi(ncv)).^2./(sigmap(ncv).*(cl(ncv).*cos(phi(ncv))+ ...
        cd(ncv).*sin(phi(ncv))))).^-1);
    
    % Identify highly loaded gridpoints (requires use of modified Glauert correction for 
    % axial induction factor)
    ncvf=find(ct>0.96*f & (da>tol | dap>tol));
    
    % Compute axial induction factor using modified Glauert correction (on identified gridpoints)
    a(ncvf)=real((18.*f(ncvf)-20-3.*sqrt(ct(ncvf).*(50-36.*f(ncvf))+12.*f(ncvf).* ...
        (3.*f(ncvf)-4)))./(36.*f(ncvf)-50));
    
    % Compute tangential induction factor
    ap(ncv)=(4.*f(ncv).*cos(phi(ncv)).*sin(phi(ncv))./(sigmap(ncv).*(cl(ncv).*sin(phi(ncv)) ...
        -cd(ncv).*cos(phi(ncv))))-1).^-1;
    
    %  Apply skewed wake correction if flow is non-axial
    if abs(gamma)>1e-8;
        a(ncv)=a(ncv).*(1+15*pi/32.*RNodes(ncv)./blade.TipRad.*tan(0.5.*(0.6.*a(ncv)+1).* ...
            gamma(ncv)).*cos(psi(ncv)));
    end
    
    % Compute residuals
    da(ncv)=abs(a0(ncv)-a(ncv));
    dap(ncv)=abs(ap0(ncv)-ap(ncv));
    
    % Apply corrective weighting for convergence stability
    if wt>0
        a(ncv)=a0(ncv)+wt.*(a(ncv)-a0(ncv));
        ap(ncv)=ap0(ncv)+wt.*(ap(ncv)-ap0(ncv));
    end

    % Clear all gridpoint flags in preparation for next loop
    clear ncv ncvf ncvcl ida idap
    
    % Identify nonconverged gridpoints
    ncv=find(da>tol | dap>tol);
    
    % If all points meet convergence criteria, break loop
    if isempty(ncv)
        break
    end
    
    % If maximum allowable iterations has been reached, flag nonconverged gridpoints 
    % with '9999'
    if j==miter
        fiter(ncv)=9999;
    else
        fiter(ncv)=j;
    end
end