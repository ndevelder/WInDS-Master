clear all
close all
clc

%define geometry of panels. note this would happen in an initialization and
%not need to be called at each time step
Np = 50; % square root of number of panels
   
x_min = -500;
x_max = 1500;
dx = (x_max-x_min)/Np;
xp = [x_min:dx:x_max]';

y_min = -375;
y_max = 375;
dy = (y_max-y_min)/Np;
yp = [y_min:dy:y_max]';

zp = zeros(Np+1,1);

[Xp,Yp,Zp] = meshgrid(xp',yp',zp');  %define the nodes of the vortex panel grid

co = 0.05; % cut off distance fraction for vortex segments 
Gamma = zeros(Np^2,1); % panel strengths

% as an example, here is how you would couple it to a time marching code
% you could paste something like what's in the for loop inside the time
% marcher in winds
for j=1:Nt-1;
            
    Gamma = VP_Influence_Coefficient_Matrix_v7(V_ext,Xp,Yp,Zp,co);  % function to calculate B and n_vec
    % note, requires a matrix V_ext of size (Np^2,3), which containts the total external velocity (including free 
    % stream and induced from blades/wake/tower at the center point of
    % every panel. 
    
    % calculate induced velocity
    V_ind = VP_Induced_Velocity_LM(P,Gamma,Xp,Yp,Zp,co);  % function to calculate induced velocity at all control points P. P is a 
    % matrix of size (N,3) where N is the total number of control points, with each row corresponding to the position vector of a 
    % control point. V_ind has the same size as P, with three velocity components per point.
                                                                                                                
    % here you would add induced velocity to other induced velocities and
    % integrate forward in time

end

