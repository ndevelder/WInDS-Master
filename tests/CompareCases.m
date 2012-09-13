clear;

Out1 = load('04-06-2012_1127_NRELrotor_static_rated.mat', '-mat', 'perf');
Out2 = load('04-06-2012_1158_NRELrotor_static_rated.mat', '-mat', 'perf');
Userout = load('04-06-2012_1127_NRELrotor_static_rated.mat', '-mat', 'user');
Bladeout = load('04-06-2012_1127_NRELrotor_static_rated.mat', '-mat', 'blade');


timestep = 1/Userout.user.t(3);
timegrid = 0:timestep:5;

tsr = Userout.user.rotor.tsr;
wind = Userout.user.rotor.wind(1);
radius = Bladeout.blade.TipRad;

%calculate two per rev time
omega = tsr*wind/radius;
maxval = max(Out1.perf.P);
opgrid = 0:100:maxval;
revpersec = omega/(2*pi);
secperrev = 1/revpersec;
twop = secperrev/2;

%Compare 2 cases
%timestep = (timelength.*step)./length(Out1.P)

diffPower = Out1.perf.P - Out2.perf.P;

plot(timegrid, diffPower,'--og');
hold on
plot(timegrid, Out1.perf.P)
hold on
plot(timegrid, Out2.perf.P,'r')
hold on
plot(twop,opgrid,'c')
