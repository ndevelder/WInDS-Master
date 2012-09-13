function wake=ab4(wake,vel,user,time)
%% wake=ab4(wake,vel,user,time) -> 4th-order Adams-Bashforth numerical integration
%
% Function numerically convects wake nodes to time+1 via 4th-order
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
if time<4 %Use forward Euler
    wake=fe(wake,vel,user,time);
else %Use 4th-order Adams-Bashforth
    wake.domain{time+1}(:,:,3:end,:)=wake.domain{time}(:,:,2:end,:)...
        +dt/24*(55*vel.domain{time}(:,:,2:end,:)...
        -59*cat(3,zeros(nst,3,1,nb),vel.domain{time-1}(:,:,2:end,:))...
        +37*cat(3,zeros(nst,3,2,nb),vel.domain{time-2}(:,:,2:end,:))...
        -9*cat(3,zeros(nst,3,3,nb),vel.domain{time-3}(:,:,2:end,:)));
end