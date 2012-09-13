function a_dyn=PittPeters(t,a,time,ct,f,U,c)
%% a_dyn=PittPeters(t,a,time,ct,f,U,c) -> Pitt and Peters ODE.
%
% Function computes the time derivative of axial induction via Pitt and
% Peters dynamic inflow ODE, as implemented in GH Bladed
% 
%
% This work is licensed under the Open Source Initiative BSD 2-Clause 
% License. To view a copy of this license, visit 
% http://www.opensource.org/licenses/BSD-2-Clause
% 
%
% Written by Thomas Sebastian (tommy.sebastian@gmail.com)
% Last edited Sept. 29, 2010
%

%%
cti=interp1(time,ct',t)'; %Interpolate for thrust coefficient wrt time
ui=interp1(time,U',t)'; %Interpolate for velocity wrt time
fi=interp1(time,f',t)'; %Interpolate for loss factor wrt time

a_dyn=ui./c.*(cti./fi-4*a.*(1-a));