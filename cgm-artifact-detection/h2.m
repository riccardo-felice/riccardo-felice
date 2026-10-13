% RICCARDO FELICE - HOMEWORK H2_25
%% sezione 1
clear all; 
close all; 
clc;

load('sensor_data_1.mat'); 
ys = sensor_data.corrupted_CGM;
ys = ys(:);                
N = length(ys);            
Ts = 5;                    
ts = (0:N-1)'*Ts;        

PH = 20;         % prediction horizon (minuti): 20, 25, 30
mu = 0.99;       % forgettinng factor: 0.86, 0.9, 0.99
p  = 2;          % ordine del modello AR: 2, 4, 6
Soglia = 20;     % soglia di errore 

L = round(PH / Ts); % passi di predizione 

%% sezione 2: inizializzazione RLS
y_pred_aligned = NaN(N, 1); % predizione allineata al tempo t
residuals = NaN(N, 1);      

% inizializzazione RLS 
P = 100 * eye(p);           
ahat = (1/p) * ones(p, 1);  

%% sezione 3: identificazione modello AR attraverso RLS + predizione
for k = (p + 1) : (N - L)
    
    % RLS
    psi = ys(k-1 : -1 : k-p); 
    y_curr = ys(k);
   
    e_rls = y_curr - psi' * ahat;  % errore a priori
    
    % aggiornamento del gain vector e la matrice di covarianza P
    K = (P * psi) / (mu + psi' * P * psi);
    P = (1/mu) * (P - K * psi' * P); % utilizzo K appena calcolato
    
    % aggiornamento parametri del modello
    ahat = ahat + K * e_rls;
    
    % predizione a L passi in avanti
    % aggiornamento di psi
    psi_temp = [y_curr; psi(1:end-1)]; 
    
    for j = 1:L
        yp_step = psi_temp' * ahat;
        psi_temp = [yp_step; psi_temp(1:end-1)]; % aggiornamento stato con nuova predizione
    end
    
    y_pred_aligned(k + L) = yp_step; % salvataggio della predizione
end

%% sezione 4: fault detection
residuals = y_pred_aligned - ys; 
alarm_idx = find(residuals > Soglia); 
num_alarm = numel(alarm_idx) % quanti allarmi con questi parametri

%% sezione 5 
figure(1)
% CGM / AR
subplot(2,1,1); 
hold on; 
grid on;
plot(ts, ys, 'b-', 'LineWidth', 1.8);             % dati CGM
plot(ts, y_pred_aligned, 'r-', 'LineWidth', 2);   % predizione AR
xline(520, 'k-');
xline(620, 'k-');
title(['Modello AR - ordine p = ', num2str(p), '  - \mu = ', num2str(mu), '  -  PH = ', num2str(PH), ' min']);
%subtitle(['\mu = ', num2str(mu), '  –  PH = ', num2str(PH), ' min']);
xlabel('Tempo (minuti)');
ylabel('Glucosio (mg/dL)');
legend('Dati CGM', 'Predizione AR', 'Location', 'best');
xlim([0 ts(end)]);
%ylim([min(ys)-10 max(ys)+10]);
ylim([70 230]);

% residuo e fault detection
subplot(2,1,2); 
hold on; 
grid on;
res = plot(ts, residuals, 'k-', 'LineWidth', 1.5);
alarm = plot(ts(alarm_idx), residuals(alarm_idx), 'ro','MarkerFaceColor', 'r', 'MarkerSize', 5);
xline(520, 'k-');
xline(620, 'k-');
soglia = yline(Soglia, 'g-', 'LineWidth', 2);
xlabel('Tempo (minuti)');
ylabel('Errore (mg/dL)');
title('Residuo di predizione');
legend([res, alarm, soglia], 'Residuo', 'Allarme', 'Soglia');
xlim([0 ts(end)]);






