function [Pc_interval , Pc_drain_max , simTimes] = PcInterval()

%This function calculates the interval for increasing thr capillary
%pressure

global pore_data throat_data

max_Pc_pore =  max(pore_data(:,28));
max_Pc_throat = max(throat_data(:,34));

min_Pc_pore = min(pore_data(:,28));
min_Pc_throat = min(throat_data(:,34));

max_Pc = max(max_Pc_pore,max_Pc_throat);
min_Pc = min(min_Pc_pore,min_Pc_throat);
Pc_interval = 35*(max_Pc - min_Pc)/1000;
P = 1.1;
Pc_drain_max = P*max_Pc;
simTimes = Pc_drain_max / Pc_interval;