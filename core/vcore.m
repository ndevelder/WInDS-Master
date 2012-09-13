function wake=vcore(wake,const,fastout,user,time)
%% wake=vcore(wake,const,fastout,user,time) -> Vortex filament core size.
%
% Function computes the effective vortex filament core size using the
% Ramasamy-Leishman model and filament stretching.
% 
% ****Input(s)****
% wake      Structure containing wake node positions, filament strengths,
%           vortex core radii, and vortex Reynolds number
% const     Structure containing model and atmospheric constants
% fastout   Structure containing time-dependent kinematics
% user      Structure containing user-defined variables
% time      Index for current timestep
% 
% ****Output(s)****
% wake      Structure containing wake node positions, filament strengths,
%           vortex core radii (updated), and vortex Reynolds number
% 
%
% This work is licensed under the Open Source Initiative BSD 2-Clause 
% License. To view a copy of this license, visit 
% http://www.opensource.org/licenses/BSD-2-Clause
% 
%
% Written by Thomas Sebastian (tommy.sebastian@gmail.com)
% Last edited February 20, 2011
%

%% Compute vortex Re #
wake.Re.shed{time}=abs(wake.gamma.shed{time}/const.nu);
wake.Re.trail{time}=abs(wake.gamma.trail{time}/const.nu);
  
%% Modify coresize using Ramasamy-Leishman model   
wake.rc.shed{time}=(wake.r0(time).^2+4*const.alpha*const.nu*(1+const.a1* ...
    wake.Re.shed{time}).*fastout.Time(time));
wake.rc.trail{time}=(wake.r0(time).^2+4*const.alpha*const.nu*(1+const.a1* ...
    wake.Re.trail{time}).*fastout.Time(time));
    
wake.rc_eff.shed{time}=wake.rc.shed{time};
wake.rc_eff.trail{time}=wake.rc.trail{time};

%% Determine filament lengths, then apply filament stretching
if strcmp(user.roll,'true')
    [wake.length.shed{time}]=FilamentLength(wake.domain{time}(1:end-1,:,:,:), ...
   wake.domain{time}(2:end,:,:,:));
    [wake.length.trail{time}]=FilamentLength(wake.domain{time}(:,:,2:end,:), ...
   wake.domain{time}(:,:,1:end-1,:));     
    
    %Effective vortex filament core size due to filament stretching between
    %current time and time-1
    wake=filamentmod(wake,time);
end