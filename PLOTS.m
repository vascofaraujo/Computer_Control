%% Time Response
figure
plot(out.estimation_error.time, out.estimation_error.signals.values  , 'LineWidth', 1.0, 'Color', 'r')
hold on
plot(out.disturbance.time, out.disturbance.signals.values, 'LineWidth', 1.0, 'Color', 'b')
xlim([0 20])
ylim([0 20])
title('LQE time response to disturbance RE=1')
legend('E. error Var1','E. error Var2','E. error Var3','E. error Var4','disturbance')
xlabel('t [s]')
ylabel('y [rad]')
hold off

%% Frequency Response
figure
Fs = 1/Ts;
L = 500;
s = Fs*(0:(L/2))/L;
y_fft = fft(yout.signals.values);
P2 = abs(y_fft/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
plot(s, P1, 'LineWidth', 1.0, 'Color', 'r')
hold on
y_fft = fft(yin.signals.values);
P2 = abs(y_fft/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
plot(s, P1, 'LineWidth', 1.0, 'Color', 'b')
title('Frequency response to square wave')
xlabel('f [Hz]')
ylabel('Y(s)')
hold off