function [pos]=kinematics(blade,turbine,platform,fastout)
%% [pos]=kinematics(blade,turbine,platform,fastout) 
%  -> Inertial position of rotor and blade stations.
%
% Function computes the station locations of each blade in the inertial
% coordinate system
% 
% ****Input(s)****
% blade     Structure containing blade geometry from FAST input file
% turbine   Structure containing turbine geometry from FAST input file
% platform  Structure containing platform geometry from FAST input file
% fastout   Structure containing imported FAST-generated results
% 
% ****Output(s)****
% pos       Structure containing relevant positions
% 
%
% This work is licensed under the Open Source Initiative BSD 2-Clause 
% License. To view a copy of this license, visit 
% http://www.opensource.org/licenses/BSD-2-Clause
% 
%
% Written by Thomas Sebastian (tommy.sebastian@gmail.com)
% Last edited March 25, 2010
%

%% Position of platform reference point in inertial coordinate system
pos.platform=[fastout.PtfmSurge fastout.PtfmSway fastout.PtfmHeave];

%% Position of rotor cone apex (hub) in inertial coordinate system
hx=turbine.OverHang*cosd(turbine.ShftTilt);
hy=0;
hz=platform.PtfmRef+turbine.TowerHt+turbine.Twr2Shft+turbine.OverHang*sind(turbine.ShftTilt);
hub_nominal=[hx hy hz]; %Coordinates of hub in ICS

%Rotation sequence for hub in ICS due to platform+nacelle motions
hub_rotseq=[fastout.PtfmYaw fastout.PtfmPitch fastout.PtfmRoll fastout.NacYaw];
hub_rotated=DCMRot(hub_nominal,hub_rotseq,[],'zyxz',0);
pos.hub=pos.platform+hub_rotated;

%% Position of spanwise stations and nodes in inertial coordinate system
nt=length(fastout.Time); %Number of timesteps
nb=turbine.NumBl; %Number of blades
nst=length(blade.RTrail); %Number of trailing nodes (+1 number of station)
ns=length(blade.RNodes); %Number of shed nodes (stations)

%Blade stations defined radially along z-axis
blade_lead=[-0.25*blade.ChordTrail zeros(nst,1) blade.RTrail];
blade_bound=[zeros(ns,1) zeros(ns,1) blade.RNodes];
blade_colloc=[0.25*blade.Chord zeros(ns,1) blade.RNodes];
blade_quarter=[zeros(nst,1) zeros(nst,1) blade.RTrail];
blade_trail=[0.75*blade.ChordTrail zeros(nst,1) blade.RTrail];
blade_end=[0.75*blade.Chord zeros(ns,1) blade.RNodes];

%Rotation sequence from rotor to inertial coordinate system
rotor_rotseq=[fastout.Azimuth turbine.ShftTilt*ones(nt,1) flipdim(hub_rotseq,2)];

%Preallocate for speed
pos.lead=zeros(nst,3,nt,nb);
pos.bound=zeros(ns,3,nt,nb);
pos.colloc=zeros(ns,3,nt,nb);
pos.quarter=zeros(nst,3,nt,nb);
pos.trail=zeros(nst,3,nt,nb);
pos.end=zeros(ns,3,nt,nb);
if strcmp(platform.Type,'EllipticalWing')
    pos.blade_rotseq=zeros(nt,9,nb);
else
    pos.blade_rotseq=zeros(nt,10,nb);
end

%Determine azimuth angle between blades, using # of blades
Azstep=360/nb;
Az=[0 cumsum(Azstep*ones(1,nb-1))];

if turbine.NumBl==2
    fastout.BldPitch(:,1)=fastout.BldPitch1;
    fastout.BldPitch(:,2)=fastout.BldPitch2;
elseif turbine.NumBl==3
    fastout.BldPitch(:,1)=fastout.BldPitch1;
    fastout.BldPitch(:,2)=fastout.BldPitch2;
    fastout.BldPitch(:,3)=fastout.BldPitch3;
end

if strcmp(platform.Type,'EllipticalWing')
    rseq='zzyxxyzxyz';
else
    rseq='zzzyxxyzxyz';
end

for c1=1:nb %Blade-specific rotation sequences
    if strcmp(platform.Type,'EllipticalWing')
        pos.blade_rotseq(:,:,c1)=[fastout.BldPitch(:,c1) turbine.PreCone(c1)*ones(nt,1) ...
            Az(c1)*ones(nt,1) rotor_rotseq];
    else
        pos.blade_rotseq(:,:,c1)=[90*ones(nt,1) fastout.BldPitch(:,c1) ...
            turbine.PreCone(c1)*ones(nt,1) Az(c1)*ones(nt,1) rotor_rotseq];
    end
    for c2=1:ns
        total_rotseq=[blade.AeroTwst(c2)*ones(nt,1) pos.blade_rotseq(:,:,c1)];
        pos.bound(c2,1:3,:,c1)=DCMRot(blade_bound(c2,:),total_rotseq,[],rseq,0)'+pos.hub';
        pos.end(c2,1:3,:,c1)=DCMRot(blade_end(c2,:),total_rotseq,[],rseq,0)'+pos.hub';
        pos.colloc(c2,1:3,:,c1)=DCMRot(blade_colloc(c2,:),total_rotseq,[],rseq,0)'+pos.hub';
    end
    for c2=1:nst
        total_rotseq=[blade.AeroTwstTrail(c2)*ones(nt,1) pos.blade_rotseq(:,:,c1)];
        pos.lead(c2,1:3,:,c1)=DCMRot(blade_lead(c2,:),total_rotseq,[],rseq,0)'+pos.hub';
        pos.quarter(c2,1:3,:,c1)=DCMRot(blade_quarter(c2,:),total_rotseq,[],rseq,0)'+pos.hub';
        pos.trail(c2,1:3,:,c1)=DCMRot(blade_trail(c2,:),total_rotseq,[],rseq,0)'+pos.hub';
    end
end   