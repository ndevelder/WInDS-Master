%% input_mod -> Modify inputs to remove discontinuties, flip signs, etc.
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
% 

%% Nodes represent section midpoints -> compute endpoints for trailing filaments
blade.RNodes=[blade.HubRad ; blade.RNodes ; blade.TipRad]; %Extrapolate to hub/tip
blade.DRNodes=diff(blade.RNodes);
blade.AeroTwst=-[blade.AeroTwst(1) ; blade.AeroTwst ; blade.AeroTwst(end)];
blade.Chord=[blade.Chord(1) ; blade.Chord ; blade.Chord(end)/2]; %Assume 50% chord at tip (from last given node)
blade.NFoil=[blade.NFoil(1) ; blade.NFoil ; blade.NFoil(end)];

if user.n~=0
    theta=(pi/2:-pi/(2*(user.n-1)):0)';
    rn=(blade.RNodes(end)-blade.RNodes(1))*cos(theta)+blade.RNodes(1);
    %rn=(blade.RNodes(1):(blade.RNodes(end)-blade.RNodes(1))/(user.n+1):blade.RNodes(end))';
    itp1=interp1(blade.RNodes,[blade.AeroTwst blade.Chord],rn,'cubic');
    itp2=interp1(blade.RNodes,blade.NFoil,rn,'nearest');
    blade.RNodes=rn;
    blade.DRNodes=diff(blade.RNodes);
    blade.AeroTwst=itp1(:,1);
    blade.Chord=itp1(:,2);
    blade.NFoil=itp2;
    clear rn itp1 itp2
end

%Interpolate for trailing filaments node locations
blade.RTrail=[blade.HubRad ; blade.RNodes(2:end-2)+blade.DRNodes(2:end-1)/2 ; blade.TipRad];
blade.AeroTwstTrail=interp1(blade.RNodes,blade.AeroTwst,blade.RTrail,'cubic');
blade.ChordTrail=interp1(blade.RNodes,blade.Chord,blade.RTrail,'cubic');

%Remove extrapolated hub/tip points for section midpoints
blade.RNodes=blade.RNodes(2:end-1);
blade.DRNodes=blade.DRNodes(2:end-1);
blade.AeroTwst=blade.AeroTwst(2:end-1);
blade.Chord=blade.Chord(2:end-1);
blade.NFoil=blade.NFoil(2:end-1);

%% Unwrap azimuth angle to remove jump discontinuities (360->0)
fastout.Azimuth=unwrap(fastout.Azimuth*pi/180)*180/pi;
id=find(diff(fastout.Azimuth)>25); %Identify remaining discontinuities (method inspired by R. Pawlowicz)
fastout.Azimuth([id ; id+1])=NaN;
bd=isnan(fastout.Azimuth); %Find "bad" points (NaN)
gd=find(~bd); %Find "good" points
bd([1:(min(gd)-1) (max(gd)+1):end])=0;
fastout.Azimuth(bd)=interp1(gd,fastout.Azimuth(gd),find(bd),'cubic');

%% Change sign/dimensions of some of the geometric inputs
for c1=1:turbine.NumBl
    eval(['fastout.BldPitch(:,c1)=fastout.BldPitch' num2str(c1) ';']);
end

%% Clear transient variables
clear id bd gd c1