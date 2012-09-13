function wake=ab2(wake,vel,user,time)
%% wake=ab2(wake,vel,user,time) -> 2nd-order Adams-Bashforth numerical integration
%
% Function numerically convects wake nodes to time+1 via 2nd-order
% Adams-Bashforth numerical integration
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

%% Determine size of test vectors/arrays
nst=size(wake.domain{1},1); %Number of trailing nodes (+1 number of station)
nb=size(wake.domain{1},4); %Number of blades
dt=1/user.t(3); %Size of timestep

%% Integrate to time+1 
if time<2 %Use forward Euler
    wake=fe(wake,vel,dt,time);
else %Use 2nd-order Adams-Bashforth
    wake.domain{time+1}(:,:,3:end,:)=wake.domain{time}(:,:,2:end,:)...
        +dt/2*(3*vel.domain{time}(:,:,2:end,:)-cat(3,zeros(nst,3,1,nb),vel.domain{time-1}(:,:,2:end,:)));
end