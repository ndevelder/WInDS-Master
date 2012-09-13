function [fastout]=output_import(filename,t)
%% [fastout]=output_import(filename,t) -> FAST-generated output importer.
%
% Function imports FAST output files and interpolates time-series data to 
% user-specifications.
%
% ****Input(s)****
% filename  String containing path to FAST input file (.fst)
% t         1x3 vector containing initial and final times and frequency
% 
% ****Output(s)****
% fastout   Structure containing imported FAST-generated results
% 
%
% This work is licensed under the Open Source Initiative BSD 2-Clause 
% License. To view a copy of this license, visit 
% http://www.opensource.org/licenses/BSD-2-Clause
% 
%
% Written by Thomas Sebastian (tommy.sebastian@gmail.com)
% Last edited February 23, 2010
%

%% Use FAST input file to ID other relevant files
data=importdata(filename,'\t'); %Import FAST input file

% Determine if Simulink-derived results or not
simq=sscanf(char(data(13)),'%i');
if simq==2 % Output file name based on use of Simulink or executable
    fastout.filename=[filename(1:end-4) '_SFunc.out'];
else
    fastout.filename=[filename(1:end-3) 'out'];
end

dt=sscanf(char(data(11)),'%f'); %Integration time step in FAST

%% Import FAST output
if exist(fastout.filename,'file')
    data=importdata(fastout.filename,'\t',7);
else
    fastout.filename=[fastout.filename(1:end-4) '_Sfunc.out'];
    data=importdata(fastout.filename,'\t',7);
    simq=2;
end
ovnames=genvarname(data(7,:)'); %Identify output variable names
odata=importdata(fastout.filename,'\t',7+1);
if simq~=2
    odata.data(:,1)=(odata.data(1,1):dt:odata.data(end,1))';
end

%% Interpolate to user-defined times
if t(3)==0 %If user-selected freq is zero, then use freq that the data is sampled at
    t(3)=1/(mean(diff(odata.data(1,:))));
end

if t(1)<odata.data(1,1);
    t(1)=odata.data(1,1);
    disp(['User selected initial time out-of-range, reset to ' num2str(t(1)) ' seconds.'])
    disp(' ')
end
if t(2)>odata.data(end,1);
    t(2)=odata.data(end,1);
    disp(['User selected final time out-of-range, reset to ' num2str(t(2)) ' seconds.'])
    disp(' ')
end
odatai=interp1(odata.data(:,1),odata.data(:,2:end),(t(1):1/t(3):t(2))');
odatai=[(t(1):1/t(3):t(2))' odatai]; %#ok<NASGU>
for c1=1:length(ovnames)
    eval(['fastout.' char(ovnames(c1)) '=odatai(:,c1);'])
end