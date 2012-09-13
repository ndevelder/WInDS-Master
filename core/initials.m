function [wake,vel,perf]=initials(pos,vel,blade,turbine,wind,airfoils,fastout,const,user)
%% [wake,vel,perf]=initials(pos,vel,blade,turbine,wind,airfoils,fastout,const,user)
%  -> Define initial values.
%
% Function preallocates memory for wake and response structures and
% variables and computes initial results for the first timestep.
% 
% ****Input(s)****
% pos       Structure containing relevant positions
% vel       Structure containing velocity components in inertial and blade
%           coordinate systems
% blade     Structure containing blade geometry
% turbine   Structure containing turbine geometry
% wind      Structure containing imported wind data
% airfoils  Structure containing airfoil performance tables
% fastout   Structure containing imported FAST-generated results
% const     Structure containing model and atmospheric constants
% user      Structure containing user-defined variables
% 
% ****Output(s)****
% wake      Structure containing wake node positions, filament strengths,
%           vortex core radii, and vortex Reynolds number
% vel       Structure containing velocity components in inertial and blade
%           coordinate systems, now including induced velocity
% perf      Structure containing performance-related variables
% 
%
% This work is licensed under the Open Source Initiative BSD 2-Clause 
% License. To view a copy of this license, visit 
% http://www.opensource.org/licenses/BSD-2-Clause
% 
%
% Written by Thomas Sebastian (tommy.sebastian@gmail.com)
% Last edited February 18, 2011
%

%% Preallocate for speed
%Determine size of test vectors/arrays
nt=length(fastout.Time); %Number of timesteps
nb=turbine.NumBl; %Number of blades
nst=length(blade.RTrail); %Number of trailing nodes (+1 number of station)
ns=length(blade.RNodes); %Number of shed nodes (stations)

wake.domain=cell(nt,1);
wake.domain(1:nt)={zeros([nst 3 nt+1 nb])};
vel.domain=cell(nt,1);
vel.domain(1:nt)={zeros([nst 3 nt+1 nb])};
vel.uind=cell(nt,1);
vel.uind(1:nt)={zeros([nst 3 nt+1 nb])};
vel.uindb=cell(nt,1);
vel.uindb(1:nt)={zeros([nst 3 nt+1 nb])};

wake.Re.shed=cell(nt,1);
wake.Re.shed(1:nt)={zeros([ns,1,nt+1,nb])};
wake.Re.trail=cell(nt,1);
wake.Re.trail(1:nt)={zeros([nst,1,nt,nb])};

wake.rc.shed=cell(nt,1);
wake.rc.shed(1:nt)={zeros([ns,1,nt+1,nb])};
wake.rc.trail=cell(nt,1);
wake.rc.trail(1:nt)={zeros([nst,1,nt,nb])};

wake.length.shed=cell(nt,1);
wake.length.shed(1:nt)={zeros([ns,1,nt+1,nb])};
wake.length.trail=cell(nt,1);
wake.length.trail(1:nt)={zeros([nst,1,nt+1,nb])};

wake.rc_eff.shed=cell(nt,1);
wake.rc_eff.shed(1:nt)={zeros([ns,1,nt+1,nb])};
wake.rc_eff.trail=cell(nt,1);
wake.rc_eff.trail(1:nt)={zeros([nst,1,nt,nb])};

wake.gamma.shed=cell(nt,1);
wake.gamma.shed(1:nt)={zeros([ns,1,nt+1,nb])};
wake.gamma.trail=cell(nt,1);
wake.gamma.trail(1:nt)={zeros([nst,1,nt+1,nb])};

perf.cl=zeros([ns,1,nt,nb]);
perf.cd=zeros([ns,1,nt,nb]);
perf.aoa=zeros([ns,1,nt,nb]);
perf.beta=zeros([ns,1,nt,nb]);

%% Substitute in initial values and truncate size of variables by timestep
for j=1:nt
    wake.domain{j}(:,:,1,:)=pos.quarter(:,:,j,:);
    wake.domain{j}(:,:,2,:)=pos.trail(:,:,j,:);
    wake.domain{j}(:,:,j+2:end,:)=[];
    
    vel.domain{j}(:,:,1:j+1,:)=repmat(wind.infty(j,:),[nst 1 j+1 nb]);
    vel.domain{j}(:,:,j+2:end,:)=[];
    vel.uind{j}(:,:,j+2:end,:)=[];  
    
    wake.Re.shed{j}(:,:,j+2:end,:)=[];  
    wake.Re.trail{j}(:,:,j+1:end,:)=[];

    wake.rc.shed{j}(:,:,j+2:end,:)=[];
    wake.rc.trail{j}(:,:,j+1:end,:)=[];
    
    wake.length.shed{j}(:,:,j+2:end,:)=[];
    wake.length.trail{j}(:,:,j+1:end,:)=[];  
    
    wake.rc_eff.shed{j}(:,:,j+2:end,:)=[];
    wake.rc_eff.trail{j}(:,:,j+1:end,:)=[];
    
    wake.gamma.shed{j}(:,:,j+2:end,:)=[];
    wake.gamma.trail{j}(:,:,j+1:end,:)=[];   
end

%% Define initial induced velocities via 1st-order methods
aoa=pos.aoag(:,1,1);
if strcmp(user.filename,'elliptical')
    cl=2*pi/(1+2/turbine.ellip.AR)*aoa*pi/180;
    perf.cl(:,1,1,1:nb)=repmat(cl,[1 1 1 nb]);
    perf.aoa(:,1,1,1:nb)=repmat(aoa,[1 1 1 nb]);
else
    [perf.bem.cl,perf.bem.cd,perf.bem.phi,perf.bem.aoa,perf.bem.a]=BEM(airfoils, ...
        blade,turbine,fastout,vel);
    perf.cl(:,1,1,1:nb)=repmat(perf.bem.cl(:,1),[1 1 1 nb]);
    perf.aoa(:,1,1,1:nb)=repmat(perf.bem.aoa(:,1),[1 1 1 nb]);
end

%% Define initial vortex strength
%Use Kutta-Joukowski theorem to define bound circulation strength
wake.gamma.shed{1}(:,:,1,:)=0.5*wind.inftyM(1).*repmat(blade.Chord,[1 1 1 nb]).* ...
    perf.cl(:,:,1,:);
%Compute spanwise change in bound filament to compute first set of trailing filaments
wake.gamma.trail{1}=diff([zeros(1,1,1,nb) ; wake.gamma.shed{1}(:,:,1,:) ; zeros(1,1,1,nb)],1);
%Shed filaments computed via spanwise summation of trailing filaments (ensure Kelvin's theorem 
%is satisfied)
wake.gamma.shed{1}(:,:,2:end,:)=diff(cat(3,cumsum(wake.gamma.trail{1}(1:end-1,:,:,:),1), ...
    zeros(ns,1,1,nb)),1,3);

%% Define initial vortex core size 
T0=2*pi*blade.TipRad./(12*fastout.TipSpdRat.*wind.inftyM);
wake.r0=sqrt(4*const.alpha*const.nu*const.delta*T0);

%% Modify core size using Ramasamy-Leishman model
wake=vcore(wake,const,fastout,user,1);

%% Compute induced velocity at all points in domain and convect points to next timestep
%Velocity induced by shed filaments on all nodes in wake
vel.uind_shed=BiotSavart(wake.domain{1}(1:end-1,:,:,:),wake.domain{1}(2:end,:,:,:), ...
    wake.domain{1},wake.gamma.shed{1},wake.rc_eff.shed{1},user.d,user.co,'full');
%Velocity induced by trailing filaments on all nodes in wake
vel.uind_trail=BiotSavart(wake.domain{1}(:,:,2:end,:),wake.domain{1}(:,:,1:end-1,:), ...
    wake.domain{1},wake.gamma.trail{1},wake.rc_eff.trail{1},user.d,user.co,'full');
%Sum the induced velocity contributions due to shed and trailing filaments
vel.uind{1}=vel.uind_shed+vel.uind_trail;
%Add the total induced velocity in the wake to the freestream velocity
vel.domain{1}=vel.domain{1}+vel.uind{1};
%Numerically convect wake nodes to time+1 via forward Euler
wake=fe(wake,vel,user,1);