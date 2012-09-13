function [vel,pos]=velocity(pos,blade,turbine,wind,fastout)
%% [vel]=velocity(blade,turbine,wind,fastout) -> Turbine motion-derived and freestream
%  velocities.
%
% Function computes the velocity contributions due to turbine and platform
% motions and freestream flow in the inertial and blade coordinate systems.
% 
% ****Input(s)****
% pos       Structure containing relevant positions 
% blade     Structure containing blade geometry from FAST input file
% turbine   Structure containing turbine geometry from FAST input file
% wind      Structure containing imported wind data
% fastout   Structure containing imported FAST-generated results
% 
% ****Output(s)****
% vel       Structure containing velocity components in inertial and blade
%           coordinate systems
% pos       Structure containing relevant positions and angles
%
% This work is licensed under the Open Source Initiative BSD 2-Clause 
% License. To view a copy of this license, visit 
% http://www.opensource.org/licenses/BSD-2-Clause
% 
%
% Written by Thomas Sebastian (tommy.sebastian@gmail.com)
% Last edited June 7, 2010
%

%% Determine size of test vectors/arrays and preallocate memory
nt=length(fastout.Time); %Number of timesteps
nb=turbine.NumBl; %Number of blades
ns=length(blade.RNodes); %Number of shed nodes (stations)

%Preallocate for speed
vel.bound=zeros(ns,3,nt,nb);
vel.blade=zeros(ns,3,nt,nb);

%% Compute kinematically-derived inertial velocities using central differencing
vel.platform=ctdiff(fastout.Time,pos.platform);
vel.hub=ctdiff(fastout.Time,pos.hub);
vel.relhub=vel.platform+vel.hub+wind.infty;

vel.bound=ctdiff(fastout.Time,pos.bound);
for c1=1:nb
    for c2=1:ns %Loop over number of blades
        vel.bound(c2,:,:,c1)=-squeeze(vel.bound(c2,:,:,c1))+wind.infty';
    end
end

%% Determine velocity in BCS via coordinate transformation (inertial to blade)
pos.nodes.bxt=pos.trail-pos.quarter;
pos.nodes.bxt=pos.nodes.bxt./repmat(sqrt(sum(pos.nodes.bxt.^2,2)),[1 3 1]);
pos.nodes.bzt=diff(pos.quarter,1,1);
pos.nodes.bzt=pos.nodes.bzt./repmat(sqrt(sum(pos.nodes.bzt.^2,2)),[1 3 1]);
pos.nodes.bzt=cat(1,pos.nodes.bzt(1,:,:,:),pos.nodes.bzt);
pos.nodes.byt=cross(pos.nodes.bzt,pos.nodes.bxt,2);

pos.nodes.bxn=pos.end-pos.bound;
pos.nodes.bxn=pos.nodes.bxn./repmat(sqrt(sum(pos.nodes.bxn.^2,2)),[1 3 1]);
pos.nodes.bzn=diff(pos.bound,1,1);
pos.nodes.bzn=pos.nodes.bzn./repmat(sqrt(sum(pos.nodes.bzn.^2,2)),[1 3 1]);
pos.nodes.bzn=cat(1,pos.nodes.bzn(1,:,:,:),pos.nodes.bzn);
pos.nodes.byn=cross(pos.nodes.bzn,pos.nodes.bxn,2);

vel.blade(:,1,:,:)=pos.nodes.bxn(:,1,:,:).*vel.bound(:,1,:,:)+pos.nodes.bxn(:,2,:,:).* ...
    vel.bound(:,2,:,:)+pos.nodes.bxn(:,3,:,:).*vel.bound(:,3,:,:);
vel.blade(:,2,:,:)=pos.nodes.byn(:,1,:,:).*vel.bound(:,1,:,:)+pos.nodes.byn(:,2,:,:).* ...
   vel.bound(:,2,:,:)+pos.nodes.byn(:,3,:,:).*vel.bound(:,3,:,:);
vel.blade(:,3,:,:)=pos.nodes.bzn(:,1,:,:).*vel.bound(:,1,:,:)+pos.nodes.bzn(:,2,:,:).* ...
   vel.bound(:,2,:,:)+pos.nodes.bzn(:,3,:,:).*vel.bound(:,3,:,:);

%% Compute geometric total angle of attack (w/o induced velocity)
pos.aoag=atan2(-vel.blade(:,2,:,:),vel.blade(:,1,:,:))*180/pi; %Geometric Total AoA
pos.aoag(isnan(pos.aoag) | abs(pos.aoag)==180)=0;