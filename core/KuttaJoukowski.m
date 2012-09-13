function [wake,perf,vel,j]=KuttaJoukowski(pos,vel,blade,turbine,wake,airfoils, ...
    user,perf,time,type,wind)
%% [wake,perf,vel]=KuttaJoukowski(pos,vel,blade,turbine,wake,airfoils,user,perf,time)
%  -> Kutta-Joukowski solver.
%
% Function computes the bound vortex filament strength via Kutta-Joukowski
% theorem, solving via fixed-point iteration or Brent's method
% 
% ****Input(s)****
% pos       Structure containing relevant positions
% vel       Structure containing velocity components in inertial and blade
%           coordinate systems
% blade     Structure containing blade geometry
% turbine   Structure containing turbine geometry
% wake      Structure containing wake node positions, filament strengths,
%           vortex core radii, and vortex Reynolds number
% airfoils  Structure containing airfoil performance tables
% user      Structure containing user-defined variables
% perf      Structure containing performance-related variables
% time      Index for current timestep
% type      If 'fixed', will use fixed-point iteration, if 'brent', will
%           use Brent's method
% 
% ****Output(s)****
% wake      Structure containing wake node positions, filament strengths
%           (updated), vortex core radii, and vortex Reynolds number
% perf      Structure containing performance-related variables (updated)
% vel       Structure containing velocity components in inertial and blade
%           coordinate systems
% 
%
% This work is licensed under the Open Source Initiative BSD 2-Clause 
% License. To view a copy of this license, visit 
% http://www.opensource.org/licenses/BSD-2-Clause
% 
%
% Written by Thomas Sebastian (tommy.sebastian@gmail.com)
% Last edited February 20, 2011
%

%% Check condition for fixed-point iteration or Brent's method
if strcmp(type,'fixed')
    j=0;
    dg=1;
    while max(max(abs(dg)))>user.tol & j<user.maxiter %#ok<AND2> %Fixed-point iteration
        gamma=wake.gamma.shed{time}(:,:,1,:);
        [dg,wake,perf,vel]=kj(gamma,vel,wake,pos,blade,turbine,perf,airfoils,time,user,wind);
        j=j+1;
    end
elseif strcmp(type,'brent')
    
%% Iteration via Brent's method
    %Preallocate for speed
    flag=zeros(6,1); %Space for logical values from conditional tests    
    na=length(airfoils.Names); %Number of airfoils
    nb=turbine.NumBl; %Number of blades    
    Vinf=sqrt(sum(vel.blade(:,:,time-1,:).^2,2)); %Magnitude of wind at the blade
    
    %Loop over airfoils + blades, interpolate wrt AoA to determine Cl and Cd
    %Adjust AoA +/-10-degrees to set upper/lower bounds for Brent's method
    aoa=perf.aoa(:,:,time-1,:);
    cla=perf.cl(:,:,time-1,:);
    clb=perf.cl(:,:,time-1,:);
    dalpha=1;
    for k=1:na
        for m=1:nb
            cla(blade.NFoil==k,1,1,m)=interp1(airfoils.profiles(k,1).AoA, ...
                airfoils.profiles(k,1).Cl,squeeze(aoa(blade.NFoil==k,1,1,m)-dalpha));
            clb(blade.NFoil==k,1,1,m)=interp1(airfoils.profiles(k,1).AoA, ...
                airfoils.profiles(k,1).Cl,squeeze(aoa(blade.NFoil==k,1,1,m)+dalpha));
        end
    end

    a=0.5*Vinf.*repmat(blade.Chord,[1 1 1 turbine.NumBl]).*cla;
    b=0.5*Vinf.*repmat(blade.Chord,[1 1 1 turbine.NumBl]).*clb;
    
    fa=kj(a,vel,wake,pos,blade,turbine,perf,airfoils,time,user,wind);
    fb=kj(b,vel,wake,pos,blade,turbine,perf,airfoils,time,user,wind);
    fs=ones(size(fb));
    
    %Check that bounds are opposite signs (soln must be between bounds)
    if any(fa(2:end-1,:,:,:).*fb(2:end-1,:,:,:)>0); 
        j=0;
        dg=1;
        while max(max(abs(dg)))>user.tol & j<user.maxiter %#ok<AND2> %Fixed-point iteration
            gamma=wake.gamma.shed{time}(:,:,1,:);
            [dg,wake,perf,vel]=kj(gamma,vel,wake,pos,blade,turbine,perf,airfoils,time,user,wind);
            j=j+1;
        end
        return
    end

    %If any values are zero (Cl=0, for example), then no sign... assign
    %+/-1 depending on the number of +/- values in bound
    if any(fa.*fb==0);
        if numel(fa<0)>numel(fa>0)
            fa(fa==0)=-1;
            fb(fb==0)=1;
        else
            fa(fa==0)=1;
            fb(fb==0)=-1;
        end
    end

    %Set |fb| < |fa|
    if abs(fa(mid(fa)))<abs(fb(mid(fb)));
        [b,a,fb,fa]=deal(a,b,fa,fb);
    end
    
    %Set initial values and conditions
    c=a;
    fc=fa;
    flag(1)=true;
    j=0;
 
%% Iterate until convergence or max. iterations reached
    while max(abs(fs))>user.tol & j<user.maxiter %#ok<AND2>
        flag(2)=all(all(fa~=fc)) && all(all(fb~=fc));
        if flag(2) %Inverse quadratic interpolation
            s=a.*fb.*fc./((fa-fb).*(fa-fc))+b.*fa.*fc./((fb-fa).*(fb-fc))+c.*fa.*fb./ ...
                ((fc-fa).*(fc-fb));
        else %Secant rule
            s=b-fb.*(b-a)./(fb-fa);
        end

        t1=0.25*(3*a+b);
        t2=b;
        if t2(mid(t2))<t1(mid(t1));
            [t2,t1]=deal(t1,t2);
        end
        
        %Conditional flags for method(s) used
        flag(3)=~(t1(mid(t1))<s(mid(s)) && s(mid(s))<t2(mid(t2)));
        flag(4)=flag(1) && abs(s(mid(s))-b(mid(b)))>=0.5*abs(b(mid(b))-c(mid(c)));
        flag(5)=~flag(1) && abs(s(mid(s))-b(mid(b)))>=0.5*abs(c(mid(c))-d(mid(d)));
        flag(6)=flag(1) && abs(b(mid(b))-c(mid(c)))<user.tol;
        flag(7)=~flag(1) && abs(c(mid(c))-d(mid(d)))<user.tol;
        
        if any(flag(3:7))
            s=0.5*(a+b); %Bisection method
            flag(1)=true;
        else
            flag(1)=false;
        end

        %Apply Kutta-Joukowski theorem to bound filament strength 's'
        [fs,wake,perf,vel]=kj(s,vel,wake,pos,blade,turbine,perf,airfoils,time,user,wind);
        s=wake.gamma.shed{time}(:,:,1,:);
        d=c;
        c=b;
        fc=fb;
        
        %Swap to set new bounds
        if any(fa.*fs<0)
            b=s;
            fb=fs;
        else
            a=s;
            fa=fs;
        end
        
        %Set |fb| < |fa|
        if abs(fa(mid(fa)))<abs(fb(mid(fb)));
            [b,a,fb,fa]=deal(a,b,fa,fb);
        end

        j=j+1;
    end
end
end


function [dg,wake,perf,vel]=kj(gamma,vel,wake,pos,blade,turbine,perf,airfoils,time,user,wind)
%% [dg,wake,perf,vel]=kj(gamma,vel,wake,pos,blade,turbine,perf,airfoils,time,user)
% -> Kutta-Joukowski theorem.
%
% Function computes the bound vortex filament strength via Kutta-Joukowski
% theorem, solving via fixed-point iteration or Brent's method
% 
% ****Input(s)****
% pos       Structure containing relevant positions
% vel       Structure containing velocity components in inertial and blade
%           coordinate systems
% blade     Structure containing blade geometry
% turbine   Structure containing turbine geometry
% wake      Structure containing wake node positions, filament strengths,
%           vortex core radii, and vortex Reynolds number
% airfoils  Structure containing airfoil performance tables
% user      Structure containing user-defined variables
% perf      Structure containing performance-related variables
% time      Index for current timestep
% 
% ****Output(s)****
% wake      Structure containing wake node positions, filament strengths
%           (updated), vortex core radii, and vortex Reynolds number
% perf      Structure containing performance-related variables (updated)
% vel       Structure containing velocity components in inertial and blade
%           coordinate systems
% 
% Written by Thomas Sebastian (tommy.sebastian@gmail.com)
% Last edited February 20, 2011
%

%% Preallocate for speed
%Determine size of test vectors/arrays
na=length(airfoils.Names); %Number of airfoils
nb=turbine.NumBl; %Number of blades
ns=length(blade.RNodes); %Number of shed nodes (stations)
cl=perf.cl(:,:,time-1,:);
cd=perf.cd(:,:,time-1,:);
vel.rot=zeros(size(vel.blade(:,:,time,:)));
wake.gamma.shed{time}(:,:,1,:)=gamma;

%% Compute induced velocity on lifting line due to shed and trailing filament induction
if strcmp(user.biotmex,'true')
[vel.uindb_shed tmpL]=BiotSavartMex(wake.domain{time}(1:end-1,:,:,:),wake.domain{time}(2:end,:,:,:), ...
    pos.bound(:,:,time,:),wake.gamma.shed{time},wake.rc_eff.shed{time},user.d,user.cmod,user.co,user.gpurun,user.gpu,user.gpuhw);
[vel.uindb_trail tmpL]=BiotSavartMex(wake.domain{time}(:,:,2:end,:),wake.domain{time}(:,:,1:end-1,:), ...
    pos.bound(:,:,time,:),wake.gamma.trail{time},wake.rc_eff.trail{time},user.d,user.cmod,user.co,user.gpurun,user.gpu,user.gpuhw);    
   %vel.uind_trail(isnan(vel.uind_trail)  | abs(vel.uind_trail)>100)=0;
   %vel.uind_shed(isnan(vel.uind_shed)  | abs(vel.uind_shed)>100)=0;
else
vel.uindb_shed=BiotSavart(wake.domain{time}(1:end-1,:,:,:),wake.domain{time}(2:end,:,:,:), ...
    pos.bound(:,:,time,:),wake.gamma.shed{time},wake.rc_eff.shed{time},user.d,user.co,'full');
vel.uindb_trail=BiotSavart(wake.domain{time}(:,:,2:end,:),wake.domain{time}(:,:,1:end-1,:), ...
    pos.bound(:,:,time,:),wake.gamma.trail{time},wake.rc_eff.trail{time},user.d,user.co,'full');
end 

uindb_shedsize = size(vel.uindb_shed);
uindb_trailsize = size(vel.uindb_trail);

if(strcmp(user.towerflag,'true'))
   vel.uind_tower = tower_shadow(squeeze(pos.bound(:,:,time,:)),wind.infty(time),[],[],user.towerupdown,turbine,time);
   towersize = size(vel.uind_tower);
   vel.uindb=vel.uindb_shed+vel.uindb_trail+vel.uind_tower; %Add in pre-calculated tower shadow uind
else
   vel.uindb=vel.uindb_shed+vel.uindb_trail;
end

%% Perform coordinate transformation on induced velocity (inertial to blade)
vel.rot(:,1,:,:)=pos.nodes.bxn(:,1,time,:).*vel.uindb(:,1,:,:)+pos.nodes.bxn(:,2,time,:).* ...
    vel.uindb(:,2,:,:)+pos.nodes.bxn(:,3,time,:).*vel.uindb(:,3,:,:);
vel.rot(:,2,:,:)=pos.nodes.byn(:,1,time,:).*vel.uindb(:,1,:,:)+pos.nodes.byn(:,2,time,:).* ...
    vel.uindb(:,2,:,:)+pos.nodes.byn(:,3,time,:).*vel.uindb(:,3,:,:);
vel.rot(:,3,:,:)=pos.nodes.bzn(:,1,time,:).*vel.uindb(:,1,:,:)+pos.nodes.bzn(:,2,time,:).* ...
vel.uindb(:,2,:,:)+pos.nodes.bzn(:,3,time,:).*vel.uindb(:,3,:,:);

%% Compute effective wind in blade coordinate system
vel.tot=vel.blade(:,:,time,:)+vel.rot;
u=vel.tot(:,1,:,:);
v=vel.tot(:,2,:,:);
w=vel.tot(:,3,:,:);

Vinf=sqrt(sum(vel.blade(:,:,time,:).^2,2));
Vtot=sqrt(sum(vel.tot.^2,2));

%% Compute angle of attack and sideslip angle
aoa=atan2(-v,u)*(180/pi);
beta=asind(w./Vtot);

%% Interpolate over airfoil data tables
for k=1:na
    for m=1:nb
        cl(blade.NFoil==k,1,1,m)=interp1(airfoils.profiles(k,1).AoA, ...
            airfoils.profiles(k,1).Cl,squeeze(aoa(blade.NFoil==k,1,1,m)));
        cd(blade.NFoil==k,1,1,m)=interp1(airfoils.profiles(k,1). ...
            AoA,airfoils.profiles(k,1).Cd,squeeze(aoa(blade.NFoil==k,1,1,m)));
    end
end

%Check for NaN values of Cl
if any(isnan(cl));
    error('Diverging soln!!!');
end

%% Compute bound vorticity via Kutta-Joukowski theorem
gamma=0.5*Vinf.*repmat(blade.Chord,[1 1 1 turbine.NumBl]).*cl;
dg=gamma-wake.gamma.shed{time}(:,:,1,:); %Change in bound vorticity between iterations

if strcmp(user.kjtype,'fixed')
    wake.gamma.shed{time}(:,:,1,:)=wake.gamma.shed{time}(:,:,1,:)+user.relax*dg;
else
    wake.gamma.shed{time}(:,:,1,:)=gamma;
end
wake.gamma.trail{time}(:,:,1,:)=diff([zeros(1,1,1,nb) ; wake.gamma.shed{time}(:,:,1,:) ; ...
    zeros(1,1,1,nb)],1);
wake.gamma.shed{time}(:,:,2:end,:)=diff(cat(3,cumsum(wake.gamma.trail{time} ...
    (1:end-1,:,:,:),1),zeros(ns,1,1,nb)),1,3);

dg=dg./(abs(gamma)+1);

%% Compute performance variables and coefficients
perf.cl(:,:,time,:)=cl;
perf.cd(:,:,time,:)=cd;
perf.aoa(:,:,time,:)=aoa;
perf.beta(:,:,time,:)=beta;

end
