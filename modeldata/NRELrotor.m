function [blade,turbine,platform,fastout,airfoils,wind]=NRELrotor(user)
%% [blade,turbine,platform,fastout,airfoils,wind]=rotor(user) -> Simple rotor test case.
%
% Function generates inputs used by WInDS to simulate the NREL 5-MW wind 
% turbine rotor
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
blade.NFoil=[1 1 1 2 3 4 4 5 6 6 7 7 8 8 8 8 8 8 8]';
blade.NFoil=interp1(blade.RTrail,blade.NFoil,blade.RNodes,'nearest');
blade.DRNodes=abs(diff(blade.RTrail));

%% Compute dynamics (rotational speed/azimuth)
omega=lambda*norm(U)/blade.TipRad;

turbine.NumBl=3;
turbine.OverHang=user.toweroverhang;
turbine.TowerHt=90;
turbine.Twr2Shft=0;
turbine.TowerBaseDiam=6.0;
turbine.TowerTopDiam=3.87;
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
airfoils.Names={'Cylinder1' ; 'Cylinder2' ; 'DU21_A17' ; 'DU25_A17' ; 'DU30_A17' ; 'DU35_A17' ; 'DU40_A17' ; 'NACA64_A17'};
for c1=1:length(airfoils.Names)
    adata=importdata(['NREL_files\' char(airfoils.Names(c1)) '.dat'],' ',14);

    id=find(diff(adata.data(:,1))==0); %ID non-distinct values for AoA
    adata.data(id,:)=[]; %#ok<FNDSB>
    
    if length(adata.textdata)==9
        eval(['airfoils.profiles(' num2str(c1) ' ,1).StallAoA=sscanf(char(adata.textdata(3)),''%f'');'])
        eval(['airfoils.profiles(' num2str(c1) ' ,1).Cn0AoA=sscanf(char(adata.textdata(4)),''%f'');'])
        eval(['airfoils.profiles(' num2str(c1) ' ,1).Lift0Cn=sscanf(char(adata.textdata(5)),''%f'');'])
        eval(['airfoils.profiles(' num2str(c1) ' ,1).StallAoACn=sscanf(char(adata.textdata(6)),''%f'');'])
        eval(['airfoils.profiles(' num2str(c1) ' ,1).StallAoANCn=sscanf(char(adata.textdata(7)),''%f'');'])
        eval(['airfoils.profiles(' num2str(c1) ' ,1).CdminAoA=sscanf(char(adata.textdata(8)),''%f'');'])
        eval(['airfoils.profiles(' num2str(c1) ' ,1).Cdmin=sscanf(char(adata.textdata(9)),''%f'');'])
    elseif length(adata.textdata)==14
        eval(['airfoils.profiles(' num2str(c1) ' ,1).StallAoA=sscanf(char(adata.textdata(5)),''%f'');'])
        eval(['airfoils.profiles(' num2str(c1) ' ,1).Cn0AoA=sscanf(char(adata.textdata(9)),''%f'');'])
        eval(['airfoils.profiles(' num2str(c1) ' ,1).Lift0Cn=sscanf(char(adata.textdata(10)),''%f'');'])
        eval(['airfoils.profiles(' num2str(c1) ' ,1).StallAoACn=sscanf(char(adata.textdata(11)),''%f'');'])
        eval(['airfoils.profiles(' num2str(c1) ' ,1).StallAoANCn=sscanf(char(adata.textdata(12)),''%f'');'])
        eval(['airfoils.profiles(' num2str(c1) ' ,1).CdminAoA=sscanf(char(adata.textdata(13)),''%f'');'])
        eval(['airfoils.profiles(' num2str(c1) ' ,1).Cdmin=sscanf(char(adata.textdata(14)),''%f'');'])
    end        

    eval(['airfoils.profiles(' num2str(c1) ' ,1).AoA=adata.data(:,1);'])
    eval(['airfoils.profiles(' num2str(c1) ' ,1).Cl=adata.data(:,2);'])
    eval(['airfoils.profiles(' num2str(c1) ' ,1).Cd=adata.data(:,3);'])
    eval(['airfoils.profiles(' num2str(c1) ' ,1).Cm=adata.data(:,4);'])
    clear adata adata1 adata2 id
end

%% Generate wind structure
wind.infty=[fastout.WindVxi fastout.WindVyi fastout.WindVzi];
wind.time=fastout.Time;
wind.inftyM=sqrt(sum(wind.infty.^2,2)); %Compute mean of freestream wind