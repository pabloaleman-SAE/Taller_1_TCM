% =========================================================================
% Taller 1: Identificación Gráfica de Sistemas FOTD
% Código corregido y simplificado con las 3 gráficas originales
% =========================================================================

clear; clc; close all;

%% 1. CARGA DE DATOS
datos = readmatrix('datos_motor.csv', 'NumHeaderLines', 1);
t = datos(:, 2); 
u = datos(:, 3); 
y = datos(:, 4); 

%% 2. ENCONTRAR EL INICIO DEL ESCALÓN
indice_escalon = 1;
for k = 1:length(u)
    if u(k) > u(1) 
        indice_escalon = k;
        break; 
    end
end
t_escalon = t(indice_escalon);

%% 3. CÁLCULO DE VALORES INICIALES, FINALES Y GANANCIA (K)
u_inicial = mean(u(1:indice_escalon-1));
u_final = mean(u(end-10:end));
y_inicial = mean(y(1:indice_escalon-1));
y_final = mean(y(end-10:end));

delta_u = u_final - u_inicial;
delta_y = y_final - y_inicial;
K = delta_y / delta_u; 

%% 4. BUSCAR TIEMPOS A LOS PORCENTAJES CLAVE (28.3% y 63.2%)
nivel_283 = y_inicial + 0.283 * delta_y;
nivel_632 = y_inicial + 0.632 * delta_y;


t_283 = 0;
for k = indice_escalon:length(y)
    if y(k) >= nivel_283
        t_283 = t(k-1) + (nivel_283 - y(k-1)) * (t(k) - t(k-1)) / (y(k) - y(k-1));
        break;
    end
end

t_632 = 0;
for k = indice_escalon:length(y)
    if y(k) >= nivel_632
        t_632 = t(k-1) + (nivel_632 - y(k-1)) * (t(k) - t(k-1)) / (y(k) - y(k-1));
        break;
    end
end

%% 5. MÉTODO ZIEGLER-NICHOLS (Recta Tangente)
dy = diff(y);
dt = diff(t);
pendiente = dy ./ dt;

% Buscamos la pendiente máxima
[m_max, idx_max] = max(pendiente(indice_escalon:end));
idx_max = idx_max + indice_escalon - 1;
t_max = t(idx_max); 
y_max = y(idx_max); 

% Ecuación de la recta despejada
t_cruce_abajo = (y_inicial - y_max) / m_max + t_max;
t_cruce_arriba = (y_final - y_max) / m_max + t_max;

theta_zc = t_cruce_abajo - t_escalon;
tau_zc = t_cruce_arriba - t_cruce_abajo;

%% 6. MÉTODO DE MILLER (63.2%)
theta_miller = max(theta_zc, 0); % aprobecha el retardo de ZN
tau_miller = (t_632 - t_escalon) - theta_miller;

%% 7. MÉTODO ANALÍTICO (Doble Punto)
tau_analitico = 1.5 * (t_632 - t_283);
theta_analitico = (t_632 - t_escalon) - tau_analitico;

%% 8. DIBUJA LAS GRÁFICAS DE LOS MODELOS
G_zc = tf(K, [tau_zc 1], 'InputDelay', theta_zc);
G_miller = tf(K, [tau_miller 1], 'InputDelay', theta_miller);
G_analitico = tf(K, [tau_analitico 1], 'InputDelay', max(theta_analitico, 0));

y_zc = lsim(G_zc, u - u_inicial, t) + y_inicial;
y_miller = lsim(G_miller, u - u_inicial, t) + y_inicial;
y_analitico = lsim(G_analitico, u - u_inicial, t) + y_inicial;

% --- GRÁFICA 1: Señales medidas ---
figure('Name', 'Senales medidas');
subplot(2, 1, 1);
plot(t, u, 'LineWidth', 1.5);
grid on; xlabel('Tiempo (s)'); ylabel('u(t)');
title('Senal de entrada');

subplot(2, 1, 2);
plot(t, y, 'k.', 'MarkerSize', 10); hold on;
plot(t, y_inicial + 0*t, '--b', 'LineWidth', 1.5);
plot(t, y_final + 0*t, '--r', 'LineWidth', 1.5);
grid on; xlabel('Tiempo (s)'); ylabel('y(t)');
title('Respuesta medida del sistema');
legend('Datos medidos', 'Valor inicial', 'Valor final', 'Location', 'southeast');

% --- GRÁFICA 2: Comparación de métodos ---
figure('Name', 'Comparacion de metodos graficos');
plot(t, y, 'k.', 'MarkerSize', 10); hold on;
plot(t, y_zc, 'b', 'LineWidth', 1.5);
plot(t, y_miller, 'r--', 'LineWidth', 1.5);
plot(t, y_analitico, 'g-.', 'LineWidth', 1.5);
plot(t_283, nivel_283, 'go', 'MarkerFaceColor', 'g');
plot(t_632, nivel_632, 'mo', 'MarkerFaceColor', 'm');
grid on; xlabel('Tiempo (s)'); ylabel('Respuesta y(t)');
title('Respuesta del sistema segun cada metodo');
legend('Datos medidos', 'ZC/tangente', 'Metodo 63.2 %', ...
       'Doble punto', '28.3 %', '63.2 %', 'Location', 'southeast');

% --- GRÁFICA 3: Construcción de los métodos ---
figure('Name', 'Construccion de los metodos');
plot(t, y, 'k.', 'MarkerSize', 10); hold on;

% Creamos la recta tangente para graficarla
y_tangente = m_max * (t - t_max) + y_max;
plot(t, y_tangente, 'b', 'LineWidth', 1.5);

plot(t_283, nivel_283, 'go', 'MarkerFaceColor', 'g');
plot(t_632, nivel_632, 'mo', 'MarkerFaceColor', 'm');
plot([t_cruce_abajo, t_cruce_arriba], [y_inicial, y_final], 'bs', 'MarkerFaceColor', 'b');

% Ajustamos la vista para que la recta larga no arruine el zoom de la gráfica
ylim([y_inicial-0.2 y_final+0.2]);

grid on; xlabel('Tiempo (s)'); ylabel('Respuesta y(t)');
title('Puntos caracteristicos y recta de maxima pendiente');
legend('Datos medidos', 'Recta de maxima pendiente', '28.3 %', ...
       '63.2 %', 'Intersecciones de la tangente', 'Location', 'southeast');