function [const]=constants
%% [const]=constants -> Constants used by WInDS.
%
% Function contains constant values associated with vortex models and with
% atmospheric conditions.
% 
% 
% ****Output(s)****
% const     Structure containing model and atmospheric constants
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

%% Constants associated with Ramasamy-Leishman vortex model
const.alpha=1.25643;
const.nu=1e-5;
const.delta=100;
const.a1=6.5e-5;

%% Atmospheric properties
const.rho=1.23; %Atmospheric density (Standard at MSL)