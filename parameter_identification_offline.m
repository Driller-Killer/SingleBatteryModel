file_name = '03-11-17_08.47 25degC_5Pulse_HPPC_Pan18650PF.mat';
load(file_name);
% Sample time  ~= 0.1s
meas.Time = 0:0.1:0.1*length(meas.Time);
meas.Time(end) = [];
meas.Time = meas.Time';
addpath('./utilities/')
%%
ax1 = subplot(211);
plot(meas.Time, meas.Voltage);
title('terminal voltage(V)');
xlabel('t [s]');
ylabel('V_{t} [V]');

ax2 = subplot(212);
plot(meas.Time, meas.Current);
title('current(A)');
xlabel('t [s]');
ylabel('I [A]');
linkaxes([ax1, ax2], 'x');

clear ax1 ax2 file_name;
%%
min_ah = -1*min(meas.Ah);
SOC = (meas.Ah + min_ah)/min_ah; %minmax scaler

clear min_ah;

plot(meas.Time, SOC)
title( 'SOC ');
xlabel('t [s]');
ylabel('SOC');
%%
Current_flanks = flanks(meas.Current, 50);
plot(meas.Time, Current_flanks, meas.Time, meas.Voltage);
legend('flanks', 'meas.Voltage');
xlabel('t[s]');
%%
n = 1;
for i = 1:length(meas.Voltage)
    if(Current_flanks(i) == -1)
        OCV(n,1) = meas.Voltage(i);% pulse discharge current begins!
    end
    if(Current_flanks(i) == 1)
        OCV(n,2) = i; %pulse discharge current ends!
        n = n+1; % move on to next pulse discharge stage!
    end
end

n = 1;
for i = 1 : length(meas.Time) 
    Vocv(i,1) = OCV(n,1);
    if ( i == OCV(n,2))
        if(n<length(OCV(:,1)))
            n = n+1;
        end
    end
end

Vocv_table = [meas.Time,Vocv];
plot(meas.Time, meas.Voltage, Vocv_table(:,1), Vocv_table(:,2));
legend('meas.Voltage', 'Vocv');
title('Voltage vs Vocv');

clear i n Vocv
%%
n = 1;
for i = 1 : length(SOC)
    if(Current_flanks(i) == -1)
        SOC_table(n,1) = SOC(i);
        n = n+1;
    end
end
out_vector = fit_func(SOC_table(:,1), 9);
coefficient = ols(out_vector, OCV(:,1));

SOC_sample = 0:0.01:1;
out_vector = fit_func(SOC_sample', 9);
OCV_fit = out_vector * coefficient;

plot(SOC_table, OCV(:,1), '.',SOC_sample, OCV_fit, '.')
legend('experiment value', 'fit value')
title('fit vs experiment SOC-OCV relationship')

clear out_vector OCV_fit i n SOC_sample
%%
indexes = struct('start',0, 'end',0);

i=1;
for n = 1 : length(Current_flanks)

    if Current_flanks(n) == -1
        indexes(i).start = n;
        if i ~= 1
            indexes(i-1).end = n;
        end
         i = i+1;
    end
end

delta = 50;

for i = 1:67
    plot(meas.Time((indexes(i).start-delta):indexes(i).end),meas.Voltage((indexes(i).start-delta):indexes(i).end));
    pause;
end

clear i n 
%% RC fit, index start from 1
% i = 1;
% param.r1 = -param.r1;
% Current_buffer_table = [meas.Time((indexes(i).start-delta):indexes(i).end), meas.Current((indexes(i).start-delta):indexes(i).end)];
% Voltage_buffer_table =  [meas.Time((indexes(i).start-delta):indexes(i).end), meas.Voltage((indexes(i).start-delta):indexes(i).end)];
%sim("singleBatteryModel.slx")


i = 4;
delta = 100;
theta = 1500;
index_1 = (indexes(i).start-delta):indexes(i).end-theta;
index_2 = (indexes(i).start-delta-1):indexes(i).end-1-theta;
index_3 = (indexes(i).start-delta-2):indexes(i).end-2-theta;

uoc_1 = fit_func(SOC(index_1), 9) * coefficient;
uoc_2 = fit_func(SOC(index_2), 9) * coefficient;
uoc_3 = fit_func(SOC(index_3), 9) * coefficient;

ut_3 = meas.Voltage(index_3);%y
ut_2 = meas.Voltage(index_2);
ut_1 = meas.Voltage(index_1);

delta2 = ut_2 - uoc_2;%x2
delta3 = ut_3 - uoc_3;%x3

i_1 = -1*meas.Current(index_1);%x4
i_2 = -1*meas.Current(index_2);%x5
i_3 = -1*meas.Current(index_3);%x6

x = [delta2, delta3, i_1, i_2, i_3];
y = -uoc_1 + ut_1;

coefficient1 = ols(x,y);

plot(meas.Time(index_1),meas.Voltage(index_1),meas.Time(index_1), uoc_1+x*coefficient1)
legend('experiment value','fit value')
%%

% R0 value is obtained averaging the R0(i) calculated from every discharge pulse for each SOC 
% Value should be close to R0 = 0.0256 in our case.
T = 0.1;
param = decoder(coefficient1, T);
%%
time_buffer = meas.Time(index_1);
current_buffer = [time_buffer i_1];
voltage_buffer = [time_buffer uoc_1];
terminal_voltage_buffer = [time_buffer ut_1];












