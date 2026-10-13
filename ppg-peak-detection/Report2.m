close all; 
clear all; 
clc

fs = 125;  % Hz


%% 1. Caricamento dati e separazione canali
 s3 = load('TROIKA Data Set-20250522\DATA_03_TYPE02.mat');
 ppg3 = s3.data(3,:);                        % solo riga 3
 acc3 = sqrt(sum(s3.data(4:6,:).^2,1));      % magnitudo accelerometro
 ecg3 = s3.data(2,:);  

s8 = load('TROIKA Data Set-20250522\DATA_08_TYPE02.mat');
ppg8 = s8.data(3,:);
acc8 = sqrt(sum(s8.data(4:6,:).^2,1));
ecg8 = s8.data(2,:);  
% acc8 = sqrt(s8.sig(4,:).^2 + s8.sig(5,:).^2 + s8.sig(6,:).^2); % PRESENTE NEL DOCUMENTO

%% 2. Pulizia PPG con filtro band-pass 0.5–8 Hz (filter)
bp = designfilt('bandpassiir','FilterOrder',4, ...
                'HalfPowerFrequency1',0.5,'HalfPowerFrequency2',8, ...
                'SampleRate',fs);
p3 = filter(bp, ppg3);
p8 = filter(bp, ppg8);

t3 = (0:length(p3)-1)/fs;
t8 = (0:length(p8)-1)/fs;

%% 3. Peak detection (vettori 1-D)

inds3_msptd   = msptd_beat_detector(p3, fs);
inds3_qppgfast =  qppgfast_beat_detector(p3(:), fs);

inds8_msptd   = msptd_beat_detector(p8, fs);
inds8_qppgfast = qppgfast_beat_detector(p8(:), fs);

%% 4. Plot sessione 3
subplot(2,1,1); hold on
  plot(t3, p3, 'k-','DisplayName','PPG filtrato');
  plot(t3(inds3_msptd), p3(inds3_msptd), 'ro','MarkerFaceColor','r','DisplayName','msptd');
  plot(t3(inds3_qppgfast), p3(inds3_qppgfast),'bd','MarkerFaceColor','b','DisplayName','qppgfast');
  xlabel('Tempo [s]'); ylabel('Ampiezza PPG');
  title('Sessione 3 — PPG con picchi rilevati');
  legend('Location','best'); grid on;
  %xlim([50 55]); % Imposta limite asse x

subplot(2,1,2);
  plot(t3, acc3);
  xlabel('Tempo [s]'); ylabel('Magnitudo Acc');
  title('Sessione 3 — Segnale accelerometro');
  %xlim([50 55]);
  grid on;

%% 5. Plot sessione 8
subplot(2,1,1); hold on
  plot(t8, p8, 'k-','DisplayName','PPG filtrato');
  plot(t8(inds8_msptd), p8(inds8_msptd), 'ro','MarkerFaceColor','r','DisplayName','msptd');
  plot(t8(inds8_qppgfast), p8(inds8_qppgfast),'bd','MarkerFaceColor','b','DisplayName','qppgfast');
  xlabel('Tempo [s]'); ylabel('Ampiezza PPG');
  title('Sessione 8 — PPG con picchi rilevati');
  legend('Location','best'); grid on;
  xlim([50 75]); % Imposta limite asse x

subplot(2,1,2);
  plot(t8, acc8);
  xlabel('Tempo [s]'); ylabel('Magnitudo Acc');
  title('Sessione 8 — Segnale accelerometro');
  grid on;
  xlim([50 75]); % Imposta limite asse x

