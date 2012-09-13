clc
clear all
close all

simlen = [2.0 5.0 10.0 20.0 30.0];

matsimlen = [2.0 5.0 10.0 20.0 30.0];

mexnogpusimlen = [2.0 5.0 10.0 20.0 30.0];

gpuwalltime = [12.783 37.589 106.271 379.250 886.314];

mexnogpuwalltime = [14.413 73.423 390.337 2596.392 8131.897];

matwalltime = [45.546 251.837 1094.687 6545.483 20056.238];

plot(simlen,gpuwalltime./60,'r');
hold on

plot(mexnogpusimlen,mexnogpuwalltime./60,'g');
hold on

plot(matsimlen,matwalltime./60,'b');

title('WInDS Computation Timings');
xlabel('Simulation Time (Seconds)');
ylabel('Computation Time (Minutes)');

l = legend ('MEX with GPU', 'MEX w/o GPU', 'MATLAB Only');

set(l,'Location', 'NorthWest');


figure(2)

gpuspeedupmat = matwalltime./gpuwalltime;

gpuspeedupmex = mexnogpuwalltime./gpuwalltime;

plot(simlen,gpuspeedupmat,'r');
hold on

plot(simlen,gpuspeedupmex,'g');
hold on

title('GPU Speedup Compared to MEX w/o GPU and MATLAB');
xlabel('Simulation Time (Seconds)');
ylabel('Times Faster (i.e. 2x)');

a = legend ('GPU vs. MATLAB', 'GPU vs. Mex Function');

set(a,'Location', 'NorthWest');