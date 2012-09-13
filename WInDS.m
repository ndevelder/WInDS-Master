%% WInDS Driver -> Wake Induced Dynamics Simulator
% 
% Driver script to compute wind turbine performance via unsteady lifting
% line method.
%
% Uses FAST input and output files to define wind turbine geometry and 
% operating conditions. WInDS then predicts wind turbine performance due 
% to wake evolution via free vortex wake method and lifting-line theory. 
%
%
% ****Function(s)****
% constants          Load constants used by other functions
% elliptical         Generate geometry and variables for elliptical wing
% rotor              Generate geometry and variables for rotor
% input_import       Import FAST-formatted input files
% output_import      Import FAST-formatted output files
% input_mod          Modify inputs, remove discontinuities
% kinematics         Compute positions of blade stations
% velocity           Compute velocity contributions due to kinematics
% initials           Set initial conditions and preallocate memory
% performance        Compute performance and load values
% 
%
% This work is licensed under the Open Source Initiative BSD 2-Clause 
% License. To view a copy of this license, visit 
% http://www.opensource.org/licenses/BSD-2-Clause
% 
%
% Written by Thomas Sebastian (tommy.sebastian@gmail.com)
% Last edited December 16, 2011
%

%% Clear command window and workspace
clear all
close all
clc

%% !!!User-defined variables!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
user.t=[0 5 10]; %Initial t, final t, and frequency in Hz
user.filename='NRELrotor'; %Test case (elliptical, rotor type, or .fst file)
user.tol=1e-8; %Tolerance value for convergence of numerical methods
user.d=0.0001; %Core model for filaments (numerical values are the cutoff radius, 
                %'viscX' applied viscous model of index X) 
user.co=1000; %Distance from wake nodes beyond which influence is negligible
user.integ='pcc'; %Numerical integration scheme
user.ns=20; %Number of spanwise stations
user.maxiter=30; %Maximum number of iterations for Kutta-Joukowski theorem
user.roll='true'; %If 'true', will apply induction to all wake nodes
user.anim='true'; %If 'true', will generate animation of wake evolution
user.time=datestr(now ,'mm-dd-yyyy_HHMM'); %Date and time of code execution
user.kjtype='fixed'; %Use either fixed point or Brent's method for convergence (Brent is
                     %still a bit coarse)
user.relax=0.25; %Relaxation value for fixed-point iteration

%Variables related to MEX/GPU computing
user.biotmex='true';
user.gpu='false'; %if 'true', use the CUDA Mex function with GPU instead of Matlab Biot-Savart
user.gpuhw=[256 1 0]; %[Thread Count, Number of GPUs, GPU Device Indices]
user.gpurun='fuln';

%Variables for tower shadow inclusion
user.towerflag='false'; %Include induced velocity from tower effects (true or false)
user.towerupdown=1; %Upwind rotor=1 Downwind rotor=2
user.toweroverhang=5.0; %Tower overhang: Should be -5.0 for NREL 5MW Upwind, positive for all downwind rotors
user.TowerBaseDiam=6.0; %Base diameter of tower: 6.0 for NREL 5MW
user.TowerTopDiam=3.87; %Top diameter of tower: 3.87 for NREL 5MW

%%Variables for user.ellip.* used only if user.filename='elliptical'
user.ellip.b=10; %Elliptical wingspan
user.ellip.AR=6; %Elliptical wing aspect ratio (AR=b^2/S)
user.ellip.wind=[1 0 0]; %Wind velocity vector
user.ellip.pitch=[5 5 0]; %Pitch angle of elliptical wing (in degrees)
user.ellip.pitchrate=0; %Pitch rate of elliptical wing (in degrees)
user.ellip.yaw=0; %Yaw angle of elliptical wing (in degrees)

%%Variables for user.rotor.* used only if user.filename='rotor'
user.rotor.wind=[11.4 0 0]; %Wind velocity vector
user.rotor.tsr=7; %Tip speed ratio
user.rotor.casetype='static_rated';
user.rotor.pitch=0; %Pitch angle of rotor blade (in degrees)
user.rotor.yaw=0;
user.rotor.modes=[];%{'Surge' 0.72520 0.00740 -1.16256 -0.44205 0.07750 2.60940 13.60156 10};

if(isnumeric(user.d))
user.cmod = 'none';
else
  if strcmp(user.gpu,'true')
  user.cmod= user.d(1:4);
  user.d = str2double(user.d(5:end));
  end
end

addpath(genpath(fullfile(cd))); %Add directories to search path

%% Load constants (physical and derived)
[const]=constants;


%% Load test case (elliptical wing, rotor, or FAST-generated)
if strcmp(user.filename,'elliptical')
    [blade,turbine,platform,fastout,airfoils,wind]=elliptical(user);
elseif strcmp(user.filename,'NRELflat')
    [blade,turbine,platform,fastout,airfoils,wind]=NRELflat(user);
elseif strcmp(user.filename,'NRELrotor')
    [blade,turbine,platform,fastout,airfoils,wind]=NRELrotor(user);
elseif strcmp(user.filename,'FAST')
    [airfoils,blade,turbine,platform,wind]=input_import(user.filename);
    [fastout]=output_import(user.filename,user.t);
end

%% Compute positions of blade stations in inertial reference frame
[pos]=kinematics(blade,turbine,platform,fastout);

%% Compute velocities of blade stations due to external motions
[vel,pos]=velocity(pos,blade,turbine,wind,fastout);

%% Define initial values (wake strength, geometry, etc)
[wake,vel,perf]=initials(pos,vel,blade,turbine,wind,airfoils,fastout,const,user);

%% !!!PRIMARY LOOP OVER TIMESERIES!!!
%Determine size of test vectors/arrays
nt=length(fastout.Time); %Number of timesteps
nb=turbine.NumBl; %Number of blades
ns=length(blade.RNodes); %Number of shed nodes (stations)
tm=zeros(nt,1); %Preallocate memory for timer (time for each timestep)

for p=2:nt
    tic; %Begin timing this timestep
%Update shed and trailing filament strength
    %Bound filament for previous timestep becomes new bound filament
    wake.gamma.shed{p}(:,:,1,:)=wake.gamma.shed{p-1}(:,:,1,:);
    %Compute spanwise change in bound filament to compute first set of trailing filaments
    wake.gamma.trail{p}(:,:,1,:)=diff([zeros(1,1,1,nb) ; wake.gamma.shed{p}(:,:,1,:) ; ...
        zeros(1,1,1,nb)],1);
    %Previous set of trailing filaments becomes new set of trailing filaments
    wake.gamma.trail{p}(:,:,2:end,:)=wake.gamma.trail{p-1};
    %Shed filaments computed via spanwise summation of trailing filaments (ensure Kelvin's 
    %theorem is satisfied)
    wake.gamma.shed{p}(:,:,2:end,:)=diff(cat(3,cumsum(wake.gamma.trail{p}(1:end-1,:,:,:),1), ...
        zeros(ns,1,1,nb)),1,3); 
    
%Modify vortex core size via Ramasamy-Leishman model and include effect of filament stretching
%from previous timestep
    wake=vcore(wake,const,fastout,user,p);
    
%Compute induced velocity at all points 
    %Velocity induced by shed filaments on all nodes in wake
    if strcmp(user.roll,'true')
        if strcmp(user.biotmex,'true')
        [vel.uind_shed tmpL]=BiotSavartMex(wake.domain{p}(1:end-1,:,:,:),wake.domain{p}(2:end,:,:,:), ...
            wake.domain{p},wake.gamma.shed{p},wake.rc_eff.shed{p},user.d,user.cmod, user.co,user.gpurun,user.gpu,user.gpuhw);
        %Velocity induced by trailing filaments on all nodes in wake
        [vel.uind_trail tmpL]=BiotSavartMex(wake.domain{p}(:,:,2:end,:),wake.domain{p}(:,:,1:end-1,:), ...
            wake.domain{p},wake.gamma.trail{p},wake.rc_eff.trail{p},user.d,user.cmod, user.co,user.gpurun,user.gpu, user.gpuhw);    
        %vel.uind_trail(isnan(vel.uind_trail)  | abs(vel.uind_trail)>100)=0;
        %vel.uind_shed(isnan(vel.uind_shed)  | abs(vel.uind_shed)>100)=0;
        else
        vel.uind_shed=BiotSavart(wake.domain{p}(1:end-1,:,:,:),wake.domain{p}(2:end,:,:,:), ...
            wake.domain{p},wake.gamma.shed{p},wake.rc_eff.shed{p},user.d,user.co,'full');
        %Velocity induced by trailing filaments on all nodes in wake
        vel.uind_trail=BiotSavart(wake.domain{p}(:,:,2:end,:),wake.domain{p}(:,:,1:end-1,:), ...
            wake.domain{p},wake.gamma.trail{p},wake.rc_eff.trail{p},user.d,user.co,'full');
        end
        %Sum the induced velocity contributions due to shed and trailing filaments
        vel.uind{p}=vel.uind_shed+vel.uind_trail;
    end 
    %Add the total induced velocity in the wake to the freestream velocity
    vel.domain{p}=vel.domain{p}+vel.uind{p};    
    
%Numerically convect wake nodes to time+1
    if strcmp(user.integ,'fe') && p~=nt
        wake=fe(wake,vel,user,p); %Foward euler
    elseif strcmp(user.integ,'ab2') && p~=nt
        wake=ab2(wake,vel,user,p); %2nd-order Adams-Bashforth
    elseif strcmp(user.integ,'ab4') && p~=nt
        wake=ab4(wake,vel,user,p); %2nd-order Adams-Bashforth
    elseif strcmp(user.integ,'pcc') && p~=nt
        wake=pcc(wake,vel,const,fastout,user,p); %Predictor-corrector, central-difference
    end

%Compute strength of new bound vortex via Kutta-Joukowski theorem
    [wake,perf,vel,ctj]=KuttaJoukowski(pos,vel,blade,turbine,wake,airfoils,user,perf,p, ...
        user.kjtype,wind);
    
%Determine time spent on current timeloop and estimate time remaining
    tm(p-1)=toc; %Time spent on current loop
    if p>2
        pt=polyfit([0 ; (2:p)'],cumsum([0 ; tm(1:p-1)]),2);
        tr=polyval(pt,nt)-sum(tm(1:p-1)); %Extrapolate to determine time remaining
        clc; disp([num2str(ctj) ': ' num2str(p/nt*100) ...
            '% complete, estimated time remaining: ' num2str(tr/60) ' minutes'])
    end
end

%% Compute performance metrics
perform;

%% Tidy up the workspace
clear yn j nb nt wb1 vs vt pg nst ns tr
save(['savedsims\' user.time '_' user.filename '_' user.rotor.casetype '.mat'])

%% Generate wake figure
if strcmp(user.anim,'true')
    j=length(fastout.Time);
    wakeplot(pos,vel,turbine,blade,wake,fastout,j);
end