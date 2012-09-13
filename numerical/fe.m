function wake=fe(wake,vel,user,time)
%% wake=ab4(wake,vel,user,time) -> Forward Euler numerical integration
%
% Function numerically convects wake nodes to time+1 via forward Euler
% numerical integration
% 
% ****Input(s)****
% wake      Structure containing wake node positions, filament strengths,
%           vortex core radii, and vortex Reynolds number
% vel       Structure containing velocity components in inertial and blade
%           coordinate systems
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

%% Integrate to time+1 via forward Euler
dt=1/user.t(3); %Size of timestep
wake.domain{time+1}(:,:,3:end,:)=wake.domain{time}(:,:,2:end,:)+dt.*vel.domain{time}(:,:,2:end,:); %Euler method