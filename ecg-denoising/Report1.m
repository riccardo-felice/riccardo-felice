clear; 
close all; 
clc;

%% 1. Generazione del segnale ECG sintetico
sfecg = 512;         % Frequenza di campionamento ECG [Hz]
sfint = 2*sfecg;     % Frequenza interna 
Nbeats = 10;         % Numero di battiti 
Anoise = 0;          % Rumore interno nullo 
hrmean = 60;         % Frequenza cardiaca media [bpm]
hrstd = 1;           % Deviazione standard della HR [bpm]
lfhfratio = 0.5;     % Rapporto LF/HF

% Genera il segnale "pulito" (artificiale)
[s, ipeaks] = ecgsyn(sfecg, Nbeats, Anoise, hrmean, hrstd, lfhfratio, sfint);

Nsample = length(s);
Tsample = 1/sfecg;
time = (0:1:Nsample-1)*Tsample;

figure;
plot(time, s, 'b');
xlabel('Tempo (s)');
ylabel('Ampiezza (mV)');
ylim([-0.6 1.25]);
xlim([0 5]);
title('Segnale ECG Sintetico');
grid on;
hold on;

% Salvo
%saveas(gcf, 'segnale_originale.png');

% Calcola il valore RMS del segnale "pulito"
rms_s = sqrt(mean(s.^2));

% Aggiunge rumore per ottenere SNR=1 (0 dB: std del rumore = rms del segnale)
noise_std = rms_s;
s_noisy = s+noise_std*randn(Nsample,1);

SNR = rms_s/noise_std % Rapporto segnale rumore
SNR = 20*log10(rms_s/noise_std); % conversione del rapporto segnale rumore in dB

figure;
plot(time, s_noisy, 'b');
xlabel('Tempo [s]'); ylabel('Ampiezza [mV]');
ylim([-2 2]);
xlim([0 5]);
title('Segnale ECG (con rumore)');
grid on;

% Salvo
%saveas(gcf, 'segnale_sintetico.png');

%% 2. Filtraggio lineare (low-pass Butterworth)

%s_noisy = s - mean(s_noisy);
% Parametri del filtro
BP = 40;           % Banda passante in Hz
BS = 1.1*BP;     % Banda di stop {1.1 2 5}
A = 100;           % Attenuazione {100 1000}
R = 0.01;          % Ripple (1%)
Rp = 20*log10(1+R);
Rs = 20*log10(A);

% Nyquist
Wp = BP/(sfecg/2);
Ws = BS/(sfecg/2);

% Filtro Butterworth con ordine limitato (a 8)
[n_filt, Wn] = buttord(Wp, Ws, Rp, Rs);
n_filt = min(n_filt, 8);  
[b, a] = butter(n_filt, Wn, 'low');

% Applica il filtro con filtfilt (evita delay)
y_filt = filtfilt(b, a, s_noisy);
y_filt = y_filt - mean(y_filt);

figure;
plot(time, s, 'b', 'LineWidth', 1.5); hold on;
plot(time, y_filt, 'r', 'LineWidth', 1.5);
xlabel('Tempo [s]'); ylabel('Ampiezza [mV]');
title('Filtraggio lineare (Butterworth) Bs=1.1*Bp, A=100');
xlim([0 5]); ylim([-0.75 1.5]);
legend('Segnale originale', 'Segnale filtrato');
grid on;

%% Plot FFT del segnale filtrato tramite Butterworth

% Rimuovo la media dal segnale filtrato
y_filt_zm = y_filt - mean(y_filt);

% Applico la finestra di Hann
w_filt = hann(length(y_filt_zm));
y_filt_windowed = y_filt_zm .* w_filt;

% Calcolo la FFT del segnale finestrato
Y_filt = fft(y_filt_windowed);

% Normalizzo la magnitudine utilizzando la somma della finestra
magY_filt = abs(Y_filt) / sum(w_filt);

% Definisco il vettore delle frequenze
frequency = (0:Nsample-1) * (1/(Nsample*Tsample));

% Plotto il grafico della FFT
figure;
plot(frequency, magY_filt, 'b', 'LineWidth', 1.5);
xlabel('Frequenza [Hz]');
ylabel('Magnitudine');
title('FFT del segnale filtrato (Butterworth)');
xlim([0 0.5*10^-3]);
grid on;
[pks_filt, locs_filt] = findpeaks(magY_filt, 'MinPeakHeight', 0.02);
hold on;
plot(frequency(locs_filt), pks_filt, 'r*', 'MarkerSize', 8);
hold off;

% Salvo
saveas(gcf, 'fft_filt.png');



%% 3. Denoising non lineare tramite FFT e soglia
s_noisy_zm = s_noisy - mean(s_noisy);

% Finestra di Hann per ridurre il spectral leakage
w = hann(length(Nsample));
s_windowed = s_noisy_zm.*w;

% Calcola la FFT del segnale windowed
Y = fft(s_windowed);

% Normalizzo la magnitudine usando sum(w)
magY = abs(Y)/sum(w);

% Definisco la soglia(es. moltiplifoc per 1, 3 o 10)
th = noise_std*sqrt(Nsample)*3; 

th_3 = noise_std*sqrt(Nsample)*3;

th_1 = noise_std*sqrt(Nsample);

% Tengo solo i coefficienti con magnitudine maggiore o uguale a th
Y_thresh = Y.*(magY >= th);

% Ricostruisco il segnale tramite ifft (solo la parte reale)
y_fft = real(ifft(Y_thresh));


[pks, locs] = findpeaks(magY, 'MinPeakHeight', 0.02);
pks = pks(2:end-1);
above_threshold = pks > th;
frequency = (0:1:Nsample-1)*(1/Nsample*Tsample);
Findex = locs(2:end-1);
figure;
plot(frequency(Findex), magY(Findex), 'b', 'LineWidth', 1.5); 
hold on;
yline(th, 'r', 'LineWidth', 1.5, 'Label', 'Soglia 10*');  % Soglia 10*
yline(th_3, 'r', 'LineWidth', 1.5, 'Label', 'Soglia 3*'); % Soglia 3*
yline(th_1, 'r', 'LineWidth', 1.5, 'Label', 'Soglia 1*'); % Soglia 1*
plot(frequency(Findex(above_threshold)), magY(Findex(above_threshold)), 'r*', 'MarkerSize', 8);
title("FFT del segnale filtrato (Non-Linear Denoising)");
xlabel('Frequenza [Hz]');
ylabel('Magnitudine');
xlim([0 0.25*10^-3]);
grid on;
hold off;

% Salvo
saveas(gcf, 'fft_3.png');

figure;
plot(time, s, 'b', 'LineWidth', 1.5); hold on;
plot(time, y_fft, 'r', 'LineWidth', 1.5);
xlabel('Tempo [s]'); ylabel('Ampiezza [mV]');
title('Denoising non lineare (FFT thresholding con finestra di Hann) (10*)');
legend('Segnale originale', 'Segnale ricostruito');
grid on;
xlim([0 5]); ylim([-1 1.5]);

% Salvo
%saveas(gcf, '10.png');

%% RMS bias, miglioramento SNR, MAE
close all

% relative RMS Bias
% ---- Caso 1: Butterworth filtrato ----
bias_filt = mean(y_filt - s);
rms_bias_filt = abs(bias_filt);
relative_rms_bias_filt = (rms_bias_filt^2) / (rms_s^2);

% ---- Caso 2: FFT thresholding ----
bias_fft = mean(y_fft - s);
rms_bias_fft = abs(bias_fft);
relative_rms_bias_fft = (rms_bias_fft^2) / (rms_s^2);

% SNR Improvement
SNR_noisy = 10*log10(sum(s.^2)/sum((s_noisy-s).^2));
SNR_filtered = 10*log10(sum(s.^2)/sum((s-y_filt).^2));
SNR_fft = 10*log10(sum(s.^2)/sum((s-y_fft).^2));

SNR_improvement_filtered = SNR_filtered-SNR_noisy;
SNR_improvement_fft = SNR_fft-SNR_noisy;

% MAE (Maximum Absolute Error)
MAE_filtered = max(abs(s-y_filt));
MAE_fft = max(abs(s-y_fft));

% Confronto
fprintf('Butterworth (Bs=1.1*Bp, A=100): RMS Bias = %.4f, SNR Improvement = %.2f dB, MAE = %.4f\n', relative_rms_bias_filt, SNR_improvement_filtered, MAE_filtered);

fprintf('Non-linear Denoising (1*): RMS Bias = %.4f, SNR Improvement = %.2f dB, MAE = %.4f\n', relative_rms_bias_fft, SNR_improvement_fft, MAE_fft);





