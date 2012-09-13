function [blade,turbine,platform,fastout,airfoils,wind]=NRELflat(user)
%% [blade,turbine,platform,fastout,airfoils,wind]=rotor(user) -> Simple rotor test case.
%
% Function generates inputs used by WInDS to simulate a simple (flat plate 
% airfoil) version of the NREL 5-MW wind turbine rotor
% 
% ****Input(s)****
% user      Structure containing user-defined variables
% 
% ****Output(s)****
% blade     Structure containing blade geometry
% turbine   Structure containing turbine geometry
% platform  Structure containing platform geometry
% fastout   Structure containing time-dependent kinematics
% airfoils  Structure containing airfoil performance tables
% wind      Structure containing wind data
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

%% User-defined variables
U=user.rotor.wind;
lambda=user.rotor.tsr;
pitch=user.rotor.pitch; %#ok<NASGU>

%% NREL 5-MW rotor geometry
blade.TipRad=63;
blade.HubRad=1.5;
blade.RTrail=[blade.HubRad 2.86670 5.6 8.33330...
    11.75 15.85 19.95 24.05 28.15 32.25 36.35...
    40.45 44.55 48.65 52.75 56.1667 58.90...
    61.6333 blade.TipRad]';
blade.ChordTrail=[3.542 3.542 3.854 4.167 4.557...
    4.652 4.458 4.249 4.007 3.748 3.502 3.256...
    3.010 2.764 2.518 2.313 2.086 1.419...
    1.419/2]';
blade.AeroTwstTrail=-[13.308 13.308 13.308 13.308...
    13.308 11.48 10.162 9.011 7.795 6.544 5.361...
    4.188 3.125 2.319 1.526 0.863 0.370 0.106...
    0.106/2]';

blade.AR=blade.TipRad^2/trapz(blade.RTrail, blade.ChordTrail);
blade.RNodes=blade.RTrail(1:end-1)+diff(blade.RTrail)/2;
blade.AeroTwst=interp1(blade.RTrail,blade.AeroTwstTrail,blade.RNodes);
blade.Chord=interp1(blade.RTrail,blade.ChordTrail,blade.RNodes);
blade.NFoil=2*ones(size(blade.RNodes));
blade.NFoil(1:3)=1;
blade.DRNodes=abs(diff(blade.RTrail));

%% Compute dynamics (rotational speed/azimuth)
omega=lambda*norm(U)/blade.TipRad;

turbine.NumBl=3;
turbine.OverHang=0;
turbine.TowerHt=90;
turbine.Twr2Shft=0;
turbine.ShftTilt=5;
turbine.PreCone=-2.5*ones(1,turbine.NumBl);

platform.Type='NRELRotor';
platform.TwrDraft=0;
platform.PtfmCm=0;
platform.PtfmRef=0;
platform.PtfmDraft=0;
platform.PtfmDiam=0;

%% Generate fastout structure
time=(user.t(1):1/user.t(3):user.t(2))';
fastout.Time=time;
fastout.WindVxi=U(1)*ones(size(time));
fastout.WindVyi=U(2)*ones(size(time));
fastout.WindVzi=U(3)*ones(size(time));
fastout.Azimuth=omega*(180/pi)*time;
for j=1:turbine.NumBl
    eval(['fastout.BldPitch' num2str(j) '=pitch*ones(size(time));']);
end
fastout.NacYaw=user.rotor.yaw*ones(size(time));

fastout.PtfmSurge=zeros(size(time));
fastout.PtfmSway=zeros(size(time));
fastout.PtfmHeave=zeros(size(time));
fastout.PtfmRoll=zeros(size(time));
fastout.PtfmPitch=zeros(size(time));
fastout.PtfmYaw=zeros(size(time));

if iscell(user.rotor.modes)
    for k=1:size(user.rotor.modes,1)
        A1=user.rotor.modes{k,2};
        f1=user.rotor.modes{k,3};
        p1=user.rotor.modes{k,4};
        A2=user.rotor.modes{k,5};
        f2=user.rotor.modes{k,6};
        p2=user.rotor.modes{k,7};
        mn=user.rotor.modes{k,8};
        tmm=time;
        tlt=find(tmm<user.rotor.modes{k,9});
        tgt=find(tmm>=user.rotor.modes{k,9});
        tmm(tlt)=0; %#ok<FNDSB>
        tmm(tgt)=tmm(tgt)-tmm(tgt(1));
        if strcmpi(user.rotor.modes{k,1},'surge')    
            fastout.PtfmSurge=A1*sin(2*pi*f1*tmm+p1)+A2*sin(2*pi*f2*tmm+p2)+mn;
        elseif strcmpi(user.rotor.modes{k,1},'sway')
            fastout.PtfmSway=A1*sin(2*pi*f1*tmm+p1)+A2*sin(2*pi*f2*tmm+p2)+mn;
        elseif strcmpi(user.rotor.modes{k,1},'heave')
            fastout.PtfmHeave=A1*sin(2*pi*f1*tmm+p1)+A2*sin(2*pi*f2*tmm+p2)+mn;
        elseif strcmpi(user.rotor.modes{k,1},'roll')
            fastout.PtfmRoll=A1*sin(2*pi*f1*tmm+p1)+A2*sin(2*pi*f2*tmm+p2)+mn;
        elseif strcmpi(user.rotor.modes{k,1},'pitch')
            fastout.PtfmPitch=A1*sin(2*pi*f1*tmm+p1)+A2*sin(2*pi*f2*tmm+p2)+mn;
        elseif strcmpi(user.rotor.modes{k,1},'yaw')
            fastout.PtfmYaw=A1*sin(2*pi*f1*tmm+p1)+A2*sin(2*pi*f2*tmm+p2)+mn;
        end
    end
end

fastout.TipSpdRat=lambda*ones(size(time));
fastout.RotSpeed=omega*(30/pi)*ones(size(time));

%% Generate airfoil database (assume drag-less cylinder and flat plate)
aoa0=0;
airfoils.Names={'zero' ; 'analytical'};
%Airfoil properties for a drag-less cylinder (zero lift and drag)
airfoils.profiles(1,1).StallAoA=0;
airfoils.profiles(1,1).Cn0AoA=0;
airfoils.profiles(1,1).Lift0Cn=0;
airfoils.profiles(1,1).Cn0AoA=0;
airfoils.profiles(1,1).StallAoACn=0;
airfoils.profiles(1,1).StallAoANCn=0;
airfoils.profiles(1,1).CdminAoA=0;
airfoils.profiles(1,1).Cdmin=0;
airfoils.profiles(1,1).AoA=(-pi:2*pi/360:pi)'*180/pi;
airfoils.profiles(1,1).Cl=zeros(size(airfoils.profiles(1).AoA));
airfoils.profiles(1,1).Cd=zeros(size(airfoils.profiles(1).AoA));
airfoils.profiles(1,1).Cm=zeros(size(airfoils.profiles(1).AoA));
%Airfoil properties for a flat plate
airfoils.profiles(2,1).StallAoA=0;
airfoils.profiles(2,1).Cn0AoA=0;
airfoils.profiles(2,1).Lift0Cn=0;
airfoils.profiles(2,1).Cn0AoA=0;
airfoils.profiles(2,1).StallAoACn=0;
airfoils.profiles(2,1).StallAoANCn=0;
airfoils.profiles(2,1).CdminAoA=0;
airfoils.profiles(2,1).Cdmin=0;
airfoils.profiles(2,1).AoA=(-2*pi:2*pi/360:2*pi)'*180/pi;
airfoils.profiles(2,1).Cl=2*pi*(airfoils.profiles(2,1).AoA-aoa0)*pi/180;
airfoils.profiles(2,1).Cl(abs(airfoils.profiles(2,1).Cl)>3)=3;
airfoils.profiles(2,1).Cd=zeros(size(airfoils.profiles(2,1).AoA));
airfoils.profiles(2,1).Cm=zeros(size(airfoils.profiles(2,1).AoA));

%% Generate wind structure
wind.infty=[fastout.WindVxi fastout.WindVyi fastout.WindVzi];
wind.time=fastout.Time;
wind.inftyM=sqrt(sum(wind.infty.^2,2)); %Compute mean of freestream wind