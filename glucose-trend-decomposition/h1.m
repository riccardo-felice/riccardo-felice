% RICCARDO FELICE - HOMEWORK H1_71
%% sezione 1
clear all; 
close all; 
clc;

load('smbg.mat');

ys = ys(:);             
N = length(ys);         
ts = (1:N)';            % asse temporale (campioni)
N_giorni = N / 3;         % ogni giorno ha 3 campioni

figure(1)
plot(ts, ys, 'r.', 'MarkerSize', 12)
title('Segnale SMBG grezzo')
xlabel('Campione'); 
ylabel('Glicemia'); 
grid on

%% sezione 2: stima del trend m(t)
w = 7;                  % (metà) ampiezza della finestra locale
m_hat = zeros(N,1);     

for k = 1:N
    % finestra centrata su k
    indice = (k-w):(k+w);
    indice = indice(indice>=1 & indice<=N); % qui impongo che la finestra non esca dai bordi
    
    % coordinate locali
    t_loc = ts(indice) - ts(k);
    y_loc = ys(indice);
    
    % pesi 
    dist_max = max(abs(t_loc)); % la distanza massima dal centro della finestra
    dist_max_norm = abs(t_loc) / dist_max; %  normalizzo le distanze
    W = (1 - dist_max_norm.^3).^3; % costruzione dei pesi
    
    % regressione locale pesata (minimi quadrati)
    G = [t_loc ones(length(t_loc),1)];
    W_radice = sqrt(W); % così integro l'errore correttamente nella regressione
    Gw = G .* W_radice;
    yw = y_loc .* W_radice;
    coefficienti = (Gw' * Gw) \ (Gw' * yw); % coefficienti della regressione
    
    % il trend 
    m_hat(k) = coefficienti(2); % valore del trend pari all'intercetta
end

% plot del trend
figure(2)
plot(ts, ys, 'r.', 'MarkerSize', 12); hold on;
plot(ts, m_hat, 'b-', 'LineWidth', 2);
title('Segnale SMBG e trend mensile m(t)')
xlabel('Campioni'); ylabel('Glicemia');
legend('y(t)','m(t)'); 
grid on;

%% sezione 3: parametri di interpolazione
m = 2;                % ordine di derivata per la regolarizzazione
passo = 0.1;
tv = (0:passo:3)'; 
nv = length(tv);      

% matrice di smoothing 
c = zeros(nv,1); c(1,1)=1; c(2,1)=-1; 
r = zeros(1,nv); r(1)=1; r(end)=0;
DELTA = toeplitz(c,r);
F = DELTA^m;
sigmaU = inv(F'*F);

% selezione dei punti realmente misurati
Gv = eye(nv); 
vt = zeros(nv,1);
vt(1) = 1; % t=0
vt(round(1/passo)+1) = 1;
vt(round(2/passo)+1) = 1; 
vt(round(3/passo)+1) = 1;  
G_local = Gv(vt==1,:);
 
% interpolazione
Interp = sigmaU * G_local' * inv(G_local * sigmaU * G_local');

%% sezione 4: estrazione componente ciclica g(t)
y_no_trend = ys - m_hat; % tolgo il trend al segnale
G_cicli = zeros(nv, N_giorni);

for k = 1:N_giorni
    % Seleziono 4 campioni (3 del giorno + 1)
    if k == 1
        indice = 1:4;       % Caso iniziale: non esiste il giorno prima, prendo il primo del giorno dopo
    elseif k == N_giorni
        indice = (N-3):N;   % Caso finale: 3 campioni di oggi + l'ultimo di ieri
    else
        indice = (3*k-3):(3*k); % Caso standard: ultimo di ieri + 3 di oggi
    end
    
    y_loc = y_no_trend(indice);
    
    % Proietto i 4 punti sulla curva 
    G_cicli(:,k) = Interp * y_loc;
end

%% sezione 5: bande di variabilità
G_ordinata = sort(G_cicli, 2);
perc_5 = G_ordinata(:, round(0.05*N_giorni));
perc_95 = G_ordinata(:, round(0.95*N_giorni));
g_media = mean(G_cicli, 2);

%% sezione 6: plot componente giornaliera
perimetro_tv = [tv; flipud(tv)];
perimetro_p5 = [perc_5; flipud(perc_95)];

figure(3); hold on;
h1 = fill(perimetro_tv, perimetro_p5, 'y', 'FaceAlpha', 0.5, 'EdgeColor', 'none'); 
h2 = plot(tv, G_cicli, 'Color', [0.7 0.7 0.7]);
h3 = plot(tv, g_media, 'Color', 'b', 'LineStyle', '-', 'LineWidth', 3); 
title('Componente giornaliera g(t) e banda di variabilità')
xlabel('Campioni giornalieri (0-3)'); 
ylabel('Glicemia');
legend([h1, h2(1), h3], 'Banda 5–95%', 'Singoli giorni', 'Media g(t)', 'Location','best');
grid on;



