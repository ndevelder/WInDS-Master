function wake=filamentmod(wake,time)
%% wake=filamentmod(wake,time) -> Core size due to filament stretching.
%
% Function computes the effective vortex filament core size due to filament
% stretching between timesteps.
% 
% ****Input(s)****
% wake      Structure containing wake node positions, filament strengths,
%           vortex core radii, and vortex Reynolds number
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

%% Apply filament stretching if time index > 3
if time>3
    trailnew=sqrt(wake.length.trail{time}(:,:,2:end-1,:));
    trailold=sqrt(wake.length.trail{time-1}(:,:,2:end,:));
    shednew=sqrt(wake.length.shed{time}(:,:,2:end-1,:));
    shedold=sqrt(wake.length.shed{time-1}(:,:,2:end,:));
    
%% Compute strain of trailing and shed filaments
    wake.strain.trail=(trailnew-trailold)./trailold;    
    wake.strain.shed=(shednew-shedold)./shedold;    
    
    %Equations modified as rc and re_eff are squared
    wake.rc_eff.trail{time}(:,:,2:end-1,:)=wake.rc.trail{time}(:,:,2:end-1,:).* ...
        (1./(1+wake.strain.trail));
    wake.rc_eff.shed{time}(:,:,2:end-1,:)=wake.rc.shed{time}(:,:,2:end-1,:).* ...
        (1./(1+wake.strain.shed));
end