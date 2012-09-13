function wake=pcc(wake,vel,const,fastout,user,time)
%% wake=pcc(wake,vel,user,time) -> %Predictor-corrector w/ central-difference
%
% Function numerically convects wake nodes to time+1 via predictor (forward
% Euler), corrector (w/ central difference) method
% 
% ****Input(s)****
% wake      Structure containing wake node positions, filament strengths,
%           vortex core radii, and vortex Reynolds number
% vel       Structure containing velocity components in inertial and blade
%           coordinate systems
% const     Structure containing model and atmospheric constants
% fastout   Structure containing time-dependent kinematics
% user      Structure containing user-defined variables
% time      Index for current timestep
% 
% ****Output(s)****
% wake      Structure containing wake node positions (time+1), filament
%           strengths, vortex core radii, and vortex Reynolds number
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

%% Determine size of test vectors/arrays
nst=size(wake.domain{1},1); %Number of trailing nodes (+1 number of station)
ns=nst-1; %Number of shed nodes (stations)
nb=size(wake.domain{1},4); %Number of blades
dt=1/user.t(3); %Size of timestep

%% Integrate to time+1
wake=fe(wake,vel,user,time); %Use forward Euler as predictor

if time>2
%Update shed and trailing filament strength
    %Bound filament for previous timestep becomes new bound filament
    wake.gamma.shed{time+1}(:,:,1,:)=wake.gamma.shed{time}(:,:,1,:);
    %Compute spanwise change in bound filament to compute first set of trailing filaments
    wake.gamma.trail{time+1}(:,:,1,:)=diff([zeros(1,1,1,nb) ; wake.gamma.shed{time+1}(:,:,1,:) ; zeros(1,1,1,nb)],1);
    %Previous set of trailing filaments becomes new set of trailing filaments
    wake.gamma.trail{time+1}(:,:,2:end,:)=wake.gamma.trail{time};
    %Shed filaments computed via spanwise summation of trailing filaments (ensure Kelvin's theorem is satisfied)
    wake.gamma.shed{time+1}(:,:,2:end,:)=diff(cat(3,cumsum(wake.gamma.trail{time+1}(1:end-1,:,:,:),1),zeros(ns,1,1,nb)),1,3);
    
%Modify vortex core size via Ramasamy-Leishman model and include effect of filament stretching from previous timestep
    wake=vcore(wake,const,fastout,user,time+1);  

%Compute induced velocity at all points
    if strcmp(user.roll,'true')
        if strcmp(user.biotmex,'true')
        %Velocity induced by shed filaments on all nodes in wake
        [vel.uind_shed tmpL]=BiotSavartMex(wake.domain{time+1}(1:end-1,:,:,:),wake.domain{time+1}(2:end,:,:,:),wake.domain{time+1},wake.gamma.shed{time+1},wake.rc_eff.shed{time+1},user.d,user.cmod,user.co,'full',user.gpu,user.gpuhw);
        %Velocity induced by trailing filaments on all nodes in wake
        [vel.uind_trail tmpL]=BiotSavartMex(wake.domain{time+1}(:,:,2:end,:),wake.domain{time+1}(:,:,1:end-1,:),wake.domain{time+1},wake.gamma.trail{time+1},wake.rc_eff.trail{time+1},user.d,user.cmod,user.co,'full',user.gpu,user.gpuhw);
        else
        %Velocity induced by shed filaments on all nodes in wake
        vel.uind_shed=BiotSavart(wake.domain{time+1}(1:end-1,:,:,:),wake.domain{time+1}(2:end,:,:,:),wake.domain{time+1},wake.gamma.shed{time+1},wake.rc_eff.shed{time+1},user.d,user.co,'full');
        %Velocity induced by trailing filaments on all nodes in wake
        vel.uind_trail=BiotSavart(wake.domain{time+1}(:,:,2:end,:),wake.domain{time+1}(:,:,1:end-1,:),wake.domain{time+1},wake.gamma.trail{time+1},wake.rc_eff.trail{time+1},user.d,user.co,'full');
        end
        %Sum the induced velocity contributions due to shed and trailing filaments
        vel.uind{time+1}=vel.uind_shed+vel.uind_trail;
    end
    %Add the total induced velocity in the wake to the freestream velocity
    vel.domain{time+1}=vel.domain{time+1}+vel.uind{time+1};
    
%Correct integration estimate via central difference step
    wake.domain{time+1}(:,:,3:end,:)=wake.domain{time}(:,:,2:end,:)...
        +dt/2*(vel.domain{time+1}(:,:,3:end,:)+vel.domain{time}(:,:,2:end,:));
end