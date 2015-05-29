
var y c k i;
varexo e s;
parameters a d;

a=0.5;
d=0.9;

model;
y=k^a+e;
y=c+i;
i=s*y;
k=(1-d)*k(-1)+i;
end;

initval;
k=0.308641975308642;
y=0.555555555555556;
i=0.277777777777778;
c=0.277777777777778;
e=0;
s=0.5;
end;

shocks;
var e; stderr 0.01;
var s; stderr 0.01;
end;

stoch_simul(order=1, periods=1000);

conditional_forecast_paths;
var y;
periods 1,2;
values 0.6,0.7;
end;
conditional_forecast(parameter_set = calibration, controlled_varexo = (e), replic = 3000);
plot_conditional_forecast(periods = 10) y i;