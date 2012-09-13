%% perform -> WInDS post-processor.
%
%
% This work is licensed under the Open Source Initiative BSD 2-Clause 
% License. To view a copy of this license, visit 
% http://www.opensource.org/licenses/BSD-2-Clause
% 
%
% Written by Thomas Sebastian (tommy.sebastian@gmail.com)
% Last edited December 17, 2010
%

%% Compute total lift coefficient
blade.S=trapz(blade.RTrail,blade.ChordTrail); %Compute planform surface area
perf.CL=1/blade.S*squeeze(trapz(blade.RNodes,repmat(blade.Chord,[1 1 nt nb]).*perf.cl)); %Compute total blade lift coefficient

%% Compute spanwise induction & forces from velocities in RCS
if turbine.NumBl>1 %Compute values if a rotor (i.e., not elliptical wing)
%Transform from blade coordinate system into rotor coordinate system (remove twist and blade pitch)
velX=zeros(size(vel.blade));
for c1=1:nb %Blade-specific rotation sequences
    for c2=1:ns
        total_rotseq=[-90*ones(nt,1) blade.AeroTwst(c2)*ones(nt,1) eval(['fastout.BldPitch' num2str(c1)])];
        velX(c2,1:3,:,c1)=-DCMRot(squeeze(vel.blade(c2,1:3,:,c1))',total_rotseq,[],'zzz',0)';
    end
end
perf.TSR=(velX(:,2,:,:)./velX(:,1,:,:)); %Compute local tip speed ratio (accounting for motion-induced velocities)

perf.phi=zeros(size(perf.aoa)); 
perf.a=zeros(size(perf.aoa));
rspd=zeros(size(perf.aoa));

for c1=1:nb
    for c2=1:ns
        for c3=1:nt
            perf.phi(c2,1,c3,c1)=perf.aoa(c2,1,c3,c1)-(eval(['fastout.BldPitch' num2str(c1) '(c3)'])'+blade.AeroTwst(c2)); %Compute inflow angle
            rspd(c2,1,c3,c1)=fastout.RotSpeed(c3);
        end
    end
end

perf.aap=perf.TSR./tand(perf.phi); %Compute a/a'
perf.a=perf.aap.*(velX(:,1,:,:)-tand(perf.phi).*velX(:,2,:,:))./(perf.aap.*velX(:,1,:,:)+tand(perf.phi).*velX(:,2,:,:)); %Compute axial induction factor
perf.ap=perf.a./perf.aap; %Compute tangential induction factor

perf.dCZ=perf.cl.*cosd(perf.phi)+perf.cd.*sind(perf.phi); %Compute normal force coefficient
perf.dCX=perf.cl.*sind(perf.phi)-perf.cd.*cosd(perf.phi); %Compute tangential force coefficient

perf.dT=0.5.*const.rho.*repmat(blade.Chord,[1 1 nt nb]).*perf.dCZ.*((velX(:,1,:,:).*(1-perf.a)).^2+(velX(:,2,:,:).*(1+perf.ap)).^2).*repmat(blade.DRNodes,[1 1 nt nb]); %Compute spanwise thrust, Newtons
perf.dQ=0.5.*const.rho.*repmat(blade.Chord,[1 1 nt nb]).*perf.dCX.*((velX(:,1,:,:).*(1-perf.a)).^2+(velX(:,2,:,:).*(1+perf.ap)).^2).*repmat(blade.RNodes,[1 1 nt nb]).*repmat(blade.DRNodes,[1 1 nt nb]); %Compute spanwise torque, Newton-meters
for c1=1:nb
    for c2=1:ns
        for c3=1:nt
            perf.dP(c2,1,c3,c1)=(pi/30)*fastout.RotSpeed(c3).*perf.dQ(c2,1,c3,c1); %Compute spanwise power, Watts
        end
    end
end

%% Compute total (rotor) performance values
perf.T=squeeze(sum(sum(perf.dT,1),4)); %Newtons
perf.Q=squeeze(sum(sum(perf.dQ,1),4)); %Newton-meters
perf.P=squeeze(sum(sum(perf.dP,1),4)); %Watts

clear c1 c2 c3 velX rspd
end