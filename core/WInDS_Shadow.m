%% WInDS Tower Shadow Velocity Deficit Code
function [V_ind] = WInDS_Shadow(P, Uinf, zinf, alpha)

% P are the node positions must be a Nx3 matrix where each row is the x,y,z position
% of a node in intertial coordinates

% Uinf is the free stream velocity. if zinf and alpha are not empty then
% use the power law to find velocity at any height. if either is empty then
% assume a uniform profile

tic

%% define geometry
zhub = 90;                     %tower height (m) arbitrary value
zbase = 0;
 
D_tower = [6; 5.75; 5.5; 5.25; 5]; % must have at least 2 entries
z_tower = [zbase; 22.5; 45; 67.5; zhub];  % z must be the same length as D to work and first and last entry must be zbase and zhub

%% initialize variables
V = zeros(size(P));
r = zeros(size(P,1),1);
theta = zeros(size(P,1),1);
R = zeros(size(P,1),1);
k = zeros(size(P,1),1);
U_r = zeros(size(P,1),1);
U_t = zeros(size(P,1),1);
u = zeros(size(P,1),1);
v = zeros(size(P,1),1);

ind_TS = find(P(:,3)<=zhub & P(:,3)>=zbase);  % only points that are between bottom and top of tower

%% define wind profile
if isempty(zinf)==0 && isempty(alpha)==0  % sheared wind flow   
    Uinf = Uinf*(P(ind_TS,3)/zinf);       
end

%% position of all relevant nodes
r(ind_TS,1) = sqrt(P(ind_TS,1).^2+P(ind_TS,2).^2);
theta(ind_TS,1) = atan2(P(ind_TS,2),P(ind_TS,1));

%% tower radius and doublet strength
R(ind_TS,1) = interp1(z_tower,D_tower,P(ind_TS,3))/2; % tower radius
k(ind_TS,1) = 2*pi*Uinf.*(R(ind_TS,1)).^2;  % doublet strength

%% Induced Velocity in Polar Coordinates (r and theta components) due to doublet 
U_r(ind_TS,1) = -(k(ind_TS,1).*cos(theta(ind_TS,1)))./(2*pi.*(r(ind_TS,1)).^2);
U_t(ind_TS,1) = -(k(ind_TS,1).*sin(theta(ind_TS,1)))./(2*pi.*(r(ind_TS,1)).^2);

%% Induced Velocity in Cartesian Coordinates (x and y components)due to doublet 
u(ind_TS,1) = U_r(ind_TS,1).*cos(theta(ind_TS,1))-U_t(ind_TS,1).*sin(theta(ind_TS,1));
v(ind_TS,1) = U_r(ind_TS,1).*sin(theta(ind_TS,1))+U_t(ind_TS,1).*cos(theta(ind_TS,1));

%% eliminate any points inside cylinder
ind_no_TS = find(r<R);  % only points that are outside of the tower

if isempty(ind_no_TS) == 0
    u(ind_no_TS,1) = zeros(length(ind_no_TS),1);
    v(ind_no_TS,1) = zeros(length(ind_no_TS),1);
end

%% final induced velocity matrix
V_ind = [u v zeros(size(P,1),1)];

toc


