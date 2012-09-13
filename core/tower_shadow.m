%% WInDS Tower Shadow Velocity Deficit Code
function [V_ind] = tower_shadow(P, Uinf_in, zinf, alpha, up_down, turbine,time)

% P are the node positions must be a Nx3 matrix where each row is the x,y,z position
% of a node in intertial coordinates
% P = P(station number, location-x-y-z, blade number)

% Uinf is the free stream velocity. if zinf and alpha are not empty then
% use the power law to find velocity at any height. if either is empty then
% assume a uniform profile


num_blade_stations = length(P);
num_dirs = 3;
num_blades = turbine.NumBl;


V_ind = zeros(num_blade_stations,num_dirs,1,num_blades); %Assumes 3 blades - change in future

%% define geometry
zhub = turbine.TowerHt;  %tower height (m) set in NRELRotor.m
zbase = 0;
 
D_tower = [turbine.TowerBaseDiam; turbine.TowerTopDiam]; % must have at least 2 entries
z_tower = [zbase; zhub];  % z must be the same length as D to work and first and last entry must be zbase and zhub


   


%% Loop over blade number

for nb=1:num_blades

    %% find only points below top of tower and above bottom of tower

    ind_TS = find(P(:,3,nb)<=zhub & P(:,3,nb)>=zbase); % only points that are between bottom and top of tower
    
    if(~isempty(ind_TS))
    P2 = P(ind_TS,:,:);

    N_P = size(P,1);
    N_P2 = size(P2,1);
    
    %% initialize variables
    V = zeros(size(P));
    r = zeros(N_P2,1);
    theta = zeros(N_P2,1);
    R = zeros(N_P2,1);
    k = zeros(N_P2,1);
    U_r = zeros(N_P2,1);
    U_t = zeros(N_P2,1);
    u_rel = zeros(N_P2,1);
    v_rel = zeros(N_P2,1);
    u = zeros(N_P,1);
    v = zeros(N_P,1);
    Uinf = Uinf_in*ones(N_P2,1);
    
    %% define wind profile
    if isempty(zinf)==0 && isempty(alpha)==0  % sheared wind flow   
        Uinf = Uinf*(P2(:,3,nb)/zinf)^alpha; 
    end

    %% position of all relevant nodes
    r(:,1) = sqrt(P2(:,1,nb).^2+P2(:,2,nb).^2);
    theta(:,1) = atan2(P2(:,2,nb),P2(:,1,nb));

    %% tower radius and doublet strength
    R(:,1) = interp1(z_tower,D_tower,P2(:,3,nb))/2; % tower radius
    k(:,1) = 2*pi*Uinf.*(R(:,1)).^2;  % doublet strength

    %% Potential Flow Solution

    % Induced Velocity in Polar Coordinates (r and theta components) due to doublet 
    U_r(:,1) = -(k(:,1).*cos(theta(:,1)))./(2*pi.*(r(:,1)).^2);
    U_t(:,1) = -(k(:,1).*sin(theta(:,1)))./(2*pi.*(r(:,1)).^2);

    % Induced Velocity in Cartesian Coordinates (x and y components)due to doublet 
    u_rel = U_r(:,1).*cos(theta(:,1))-U_t(:,1).*sin(theta(:,1));
    v_rel = U_r(:,1).*sin(theta(:,1))+U_t(:,1).*cos(theta(:,1));

    %% downwind tower shadow solution
    if up_down == 2
         
        %Get indices of P2(Y) less than 2.5x tower diameter and P2(X) greater than tower radius 
        ind_DW_TS = find(P2(:,2,nb)<2.5*R(:,1)*2 & P2(:,1,nb)>R(:,1));  % downwind points within the wake

        if isempty(ind_DW_TS) == 0       
            u_rel(ind_DW_TS,1,nb) = -0.3*Uinf(ind_DW_TS,1);
            v_rel(ind_DW_TS,1,nb) = 0;       
        end
    end

    %% eliminate any points inside cylinder
    ind_no_TS = find(r<R);  % only points that are outside of the tower

    N = length(ind_no_TS);

    if isempty(ind_no_TS) == 0
        u_rel(ind_no_TS,1) = zeros(N,1);
        v_rel(ind_no_TS,1) = zeros(N,1);
    end

    %% All velocities 
    u(ind_TS,1) = u_rel(:,1);
    v(ind_TS,1) = v_rel(:,1);


    %% final induced velocity matrix
    V_ind(:,1,1,nb) = u(:,1);
    V_ind(:,2,1,nb) = v(:,1);
    V_ind(:,3,1,nb) = zeros(num_blade_stations,1);

    else
    %vindsize = size( V_ind(:,1,1,nb))  
    V_ind(:,1,1,nb) = zeros(num_blade_stations,1);
    V_ind(:,2,1,nb) = zeros(num_blade_stations,1);
    V_ind(:,3,1,nb) = zeros(num_blade_stations,1);
        
    end
end
