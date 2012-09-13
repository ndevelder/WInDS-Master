clear
clc

addpath(genpath(fullfile(cd))); %Add directories to search path

P = ones(20,3,100,3);
F1 = ones(20,3,100,3);
F2 = ones(20,3,100,3);
gamma = ones(20,3,100,3);
rc = ones(20,3,100,3);
co = ones(20,3,100,3);
casetype = 'leng';
d = 'visc';
flag = 'mex';
P(1,1,1,1) = 5.00;
F1(1,1,1,1) = 5.00;
F1(1,2,1,1) = 5.00;
F1(1,3,1,1) = 5.00;
P(1,2,1,2) = 16.000;

if(flag == 'mex')
[uind L] = BiotSavartMex(F1,F2,P,gamma,rc,d,co,casetype);
end

if(flag == 'm')
[uind L] = BiotSavart(F1,F2,P,gamma,rc,d,co,casetype);
end