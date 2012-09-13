function [blade,turbine,platform,fastout,airfoils,wind]=elliptical(user)
%% [blade,turbine,platform,fastout,airfoils,wind]=elliptical(user) -> Elliptical wing test case.
%
% Function generates inputs used by WInDS to simulate an elliptical wing
% test case
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
b=user.ellip.b;
AR=user.ellip.AR;
U=user.ellip.wind;
ptch=user.ellip.pitch(1);
yw=user.ellip.yaw;
n=user.ns;

%% Elliptical wing geometry
S=b^2/AR; %Planform area of wing
c=4*S/(pi*b); %Root chord of wing

%Generate station locations via cosine distribution (more points near wingtips)
theta=(-pi:pi/n:0)';
x=(b/2)*cos(theta);
t=acos(2*x/b);
y=c*sin(t);

blade.TipRad=b/2;
blade.HubRad=-b/2;
blade.RTrail=x;
blade.ChordTrail=y;
blade.AeroTwstTrail=zeros(size(blade.RTrail));

blade.RNodes=blade.RTrail(1:end-1)+diff(blade.RTrail)/2;
blade.AeroTwst=interp1(blade.RTrail,blade.AeroTwstTrail,blade.RNodes);
blade.Chord=interp1(blade.RTrail,blade.ChordTrail,blade.RNodes);
blade.NFoil=ones(size(blade.RNodes));
blade.DRNodes=abs(diff(blade.RTrail));

turbine.NumBl=1;
turbine.OverHang=0;
turbine.TowerHt=0;
turbine.Twr2Shft=0;
turbine.ShftTilt=0;
turbine.PreCone=0;

turbine.ellip.b=b;
turbine.ellip.S=S;
turbine.ellip.AR=AR;

platform.Type='EllipticalWing';
platform.TwrDraft=0;
platform.PtfmCm=0;
platform.PtfmRef=0;
platform.PtfmDraft=0;
platform.PtfmDiam=0;

%% Generate fastout structure (most values are zero)
time=(user.t(1):1/user.t(3):user.t(2))';
fastout.Time=time;
fastout.WindVxi=U(1)*ones(size(time));
fastout.WindVyi=U(2)*ones(size(time));
fastout.WindVzi=U(3)*ones(size(time));
fastout.Azimuth=-90*ones(size(time));

if user.ellip.pitchrate~=0
    ptch1=user.ellip.pitch(1);
    ptch2=user.ellip.pitch(2);
    ptime=user.ellip.pitch(3);
    prate=user.ellip.pitchrate;
    ptch=interp1([time(1) ptime ptime+abs(ptch2-ptch1)/prate time(end)],[ptch1 ptch1 ptch2 ptch2],time);
end
fastout.BldPitch=ptch.*ones(size(time));
fastout.NacYaw=0*ones(size(time));
fastout.PtfmSurge=0*ones(size(time));
fastout.PtfmSway=0*ones(size(time));
fastout.PtfmHeave=zeros(size(time));
fastout.PtfmRoll=0*ones(size(time));
fastout.PtfmPitch=0*ones(size(time));
fastout.PtfmYaw=yw*ones(size(time));
fastout.TipSpdRat=ones(size(time));
fastout.RotSpeed=zeros(size(time));

%% Generate airfoil database (assume flat plate)
aoa0=0;
airfoils.Names={'Analytical'};
airfoils.Profiles(1).StallAoA=0;
airfoils.Profiles(1).Cn0AoA=0;
airfoils.Profiles(1).Lift0Cn=0;
airfoils.Profiles(1).Cn0AoA=0;
airfoils.profiles(1).StallAoACn=0;
airfoils.profiles(1).StallAoANCn=0;
airfoils.profiles(1).CdminAoA=0;
airfoils.profiles(1).Cdmin=0;
airfoils.profiles(1).AoA=(-pi/2:pi/100:pi/2)'*180/pi;
airfoils.profiles(1).Cl=2*pi*(airfoils.profiles(1).AoA-aoa0)*pi/180;
airfoils.profiles(1).Cd=zeros(size(airfoils.profiles(1).AoA));
airfoils.profiles(1).Cm=zeros(size(airfoils.profiles(1).AoA));

%% Generate wind structure
wind.infty=[fastout.WindVxi fastout.WindVyi fastout.WindVzi];
wind.time=fastout.Time;
wind.inftyM=sqrt(sum(wind.infty.^2,2));