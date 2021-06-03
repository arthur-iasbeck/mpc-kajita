%#ok<*UNRCH>
clc; clear; close all;

% Variáveis de controle
plot_ref = 0;
print_resolution = '-r300';

% Variáveis associadas ao modelo
g = 9.81;               % Aceleração gravitacional
h = 0.24;               % Altura do centro de massa
T = 0.01;               % Tempo de amostragem

% Variáveis associadas à função custo
alpha = 1e-6;           % Peso da sobre-aceleração
gamma = 1;              % Peso do rastreamento do ZMP de referência

% Variáveis associadas à simulação
N = 300;                % Janela de predição
S = 50;                 % Tempo do passo/Tempo de amostragem
n_step = 12;            % Número de passos na simulação
len_step_x = 0.05;      % Tamanho do passo em x em metros
len_step_y = 0.1;       % Tamanho do passo em y em metros
w_foot = 0.06;          % Comprimento do pé
h_foot = 0.03;          % Largura do pé

% Definição do modelo
A = [1, T, T^2/2;
    0, 1, T;
    0, 0, 1];

B = [T^3/6;
    T^2/2;
    T];

A_z = [1, 0, -h/g];

% NOTA: Equações associadas ao modelo
% x_s = [x dx d2x]
% dx_s = A*x_s + B*d3x
% z_x = A_z*x_s

% Computação da matriz Pzu
P_zu = zeros(N, N);
for i = 1:N
    m = i;
    for j = 1:N
        if j > i
            break
        end
        factor = 1 + 3*(m - 1) + 3*(m - 1)^2;
        m = m - 1;
        P_zu(i, j) = factor*T^3/6 - h*T/g;
    end
end

% Computação da matriz Pzs
P_zs = zeros(N, 3);
for i = 1:N
    P_zs(i, :) = [1, i*T, i^2*T^2/2 - h/g];
end

% Construção da referência do ZMP
Z_rx_full = zeros(n_step*S + 1, 1);
Z_ry_full = len_step_y*ones(n_step*S + 1, 1);
i_0 = 2;
i_f = i_0 + S - 1;
zmp_y = len_step_y;
for i = 1:n_step
    Z_rx_full(i_0:i_f) = len_step_x*(i-1)*ones(S,1);
    Z_ry_full(i_0:i_f) = zmp_y*ones(S,1);
    i_0 = i_f + 1;
    i_f = i_0 + S - 1;
    if zmp_y == len_step_y
        zmp_y = 0;
    elseif zmp_y == 0
        zmp_y = len_step_y;
    end
end

% Representação da referência ao longo do tempo
figure; hold on;
t = (0:length(Z_rx_full)-1)*T;
plot(t, Z_rx_full);
plot(t, Z_ry_full);
legend({'z_x', 'z_y'}, 'Location', 'northwest');

% Definição da posição dos passos
step_pos_x = zeros(n_step, 1);
step_pos_y = zeros(n_step, 1);

for i = 1:n_step+1
    step_pos_x(i) = (i-1)*len_step_x;
    step_pos_y(i) = mod(i,2)*len_step_y;
end

% Inicialização de variáveis utilizadas na simulação
x_s = zeros(3, 1);
y_s = zeros(3, 1);
d3x = zeros(1, 1);
d3y = zeros(1, 1);
z_x = zeros(1, 1);
z_y = zeros(1, 1);

% Definição do estado inicial
x_s(:, 1) = [step_pos_x(1) 0 0]';
y_s(:, 1) = [step_pos_y(1) 0 0]';
z_x(:, 1) = step_pos_x(1);
z_y(:, 1) = step_pos_y(1);

% Loop de controle
k = 1;
while 1
    % Verificação da interrupção do laço
    if k+N-1 > length(Z_rx_full)
        break;
    end
    
    % Avaliação da referência do ZMP
    Z_rx = Z_rx_full(k:k+N-1);
    Z_ry = Z_ry_full(k:k+N-1);
    
    % Computação dos próximos esforços de controle
    d3X = (P_zu'*P_zu + (alpha/gamma)*eye(N))\...
        (P_zu'*(Z_rx - P_zs*x_s(:,k)));
    
    d3Y = (P_zu'*P_zu + (alpha/gamma)*eye(N))\...
        (P_zu'*(Z_ry - P_zs*y_s(:,k)));
    
    % Definição do esforço de controle a ser aplicado
    d3x(k) = d3X(1);
    d3y(k) = d3Y(1);
    
    % Evolução do modelo
    x_s(:, k+1) = A*x_s(:,k) + B*d3x(k);
    y_s(:, k+1) = A*y_s(:,k) + B*d3y(k);
    z_x(k+1) = A_z*x_s(:,k);
    z_y(k+1) = A_z*y_s(:,k);
    
    k = k + 1;
end

x_s(:,k) = [];
y_s(:,k) = [];
z_x(k) = [];
z_y(k) = [];
k = k - 1;


%% Representação gráfica dos resultados
k_max = k;

close all;
t = (0:k_max-1)*T;

% Rastreamento da referência do ZMP em x
figure;
subplot(2,1,1);
hold on;
stairs(t, z_x, '-', 'MarkerSize', 15, 'LineWidth', 2);
stairs(t, Z_rx_full(1:k_max), '-', 'MarkerSize', 15, 'LineWidth', 2);
legend({'z_x', 'z_{rx}'}, 'Location', 'eastoutside');
grid on; axis tight;
ax = gca; ax.FontSize = 13;
xlabel('t [s]');
ylabel('x [m]');
set(gcf, 'Position', get(0, 'Screensize'));

% Rastreamento da referência do ZMP em y
subplot(2,1,2);
hold on;
stairs(t, z_y, '-', 'MarkerSize', 15, 'LineWidth', 2);
stairs(t, Z_ry_full(1:k_max), '-', 'MarkerSize', 15, 'LineWidth', 2);
legend({'z_y', 'z_{ry}'}, 'Location', 'eastoutside');
grid on; axis tight;
ax = gca; ax.FontSize = 13;
xlabel('t [s]');
ylabel('y [m]');
set(gcf, 'Position', get(0, 'Screensize'));
print('ZMP_t', '-dpng', print_resolution);

% Evolução da velocidade do centro de massa
figure; hold on;
hold on;
stairs(t, x_s(2,:), '-', 'MarkerSize', 15, 'LineWidth', 2);
stairs(t, y_s(2,:), '-', 'MarkerSize', 15, 'LineWidth', 2);
legend({'v_x', 'v_y'}, 'Location', 'eastoutside');
grid on; axis tight;
ax = gca; ax.FontSize = 13;
xlabel('t [s]');
ylabel('v [m/s]');
set(gcf, 'Position', get(0, 'Screensize'));
print('CoM_v', '-dpng', print_resolution);

% Evolução da aceleração do centro de massa
figure; hold on;
hold on;
stairs(t, x_s(3,:), '-', 'MarkerSize', 15, 'LineWidth', 2);
stairs(t, y_s(3,:), '-', 'MarkerSize', 15, 'LineWidth', 2);
legend({'a_x', 'a_y'}, 'Location', 'eastoutside');
grid on; axis tight;
ax = gca; ax.FontSize = 13;
xlabel('t [s]');
ylabel('a [m/s²]');
set(gcf, 'Position', get(0, 'Screensize'));
print('CoM_a', '-dpng', print_resolution);

% Evolução da sobre-aceleração do centro de massa
figure; hold on;
hold on;
stairs(t, d3x, '-', 'MarkerSize', 15, 'LineWidth', 2);
stairs(t, d3y, '-', 'MarkerSize', 15, 'LineWidth', 2);
legend({'j_x', 'j_y'}, 'Location', 'eastoutside');
grid on; axis tight;
ax = gca; ax.FontSize = 13;
xlabel('t [s]');
ylabel('j [m/s³]');
set(gcf, 'Position', get(0, 'Screensize'));
print('CoM_j', '-dpng', print_resolution);

% Representação dos passos
figure; hold on; axis equal; grid on;
plot(step_pos_x, step_pos_y, 'k.', 'MarkerSize', 20),
for i = 1:length(step_pos_x)
    x_foot = step_pos_x(i) - w_foot/2;
    y_foot = step_pos_y(i) - h_foot/2;
    rectangle('Position', [x_foot, y_foot, w_foot, h_foot]);
end

% Evolução do ZMP e do CoM
plot_obj(1) = plot(x_s(1,:), y_s(1,:), 'LineWidth', 2);
plot_obj(2) = plot(z_x, z_y, 'LineWidth', 2);
plot_obj(3) = plot(Z_rx_full, Z_ry_full, '--', 'LineWidth', 2, 'MarkerSize', 13);
legend(plot_obj, {'CoM', 'ZMP', 'ZMP_r'}, 'Location', 'eastoutside');
ax = gca; ax.FontSize = 13;
xlabel('x [m]');
ylabel('y [m]');
set(gcf, 'Position', get(0, 'Screensize'));
xlim([min(z_x) - w_foot/2, max(z_x)])
print('ZMP_CoM', '-dpng', print_resolution);

fprintf('Simulação finalizada \n');


