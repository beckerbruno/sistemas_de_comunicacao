function ta2_ldpc_qam_ofdm()
% TA2 - Sistemas de Comunicacao 2026/1
% Projeto pratico: LDPC com 16-QAM e OFDM
% PUCRS - Escola Politecnica
% Trilha 1: LDPC didatico sem toolbox (matriz pequena 3x6)
%
% Objetivo: Comparar sistema QAM-OFDM sem FEC com sistema QAM-OFDM com FEC LDPC
% Entregaveis:
%   1. Codigo MATLAB comentado (este arquivo)
%   2. Graficos de BER versus SNR
%   3. Constelacoes QAM antes e depois do canal
%   4. Espectro aproximado do sinal OFDM
%   5. Tabela comparando sistema sem FEC e com LDPC
%   6. Relatorio (arquivo LaTeX separado)

%% ============================================================
% PARAMETROS DO SISTEMA
%% ============================================================

% Parametros LDPC (matriz didatica do exemplo do Capitulo 7)
% Codigo (6,3): n=6 bits codificados, k=3 bits de informacao, taxa R=1/2
H = [1 1 1 1 0 0; ...  % Equacao de paridade 1
     0 1 0 1 0 1; ...  % Equacao de paridade 2
     1 1 0 0 0 1];     % Equacao de paridade 3

% Calcular matriz geradora G via metodo direto
% H e m x n, queremos G de dimensao k x n onde H*G' = 0 (mod 2)
[m, n] = size(H);
k = n - m;  % Numero de bits de informacao = 3

% Para matriz H pequena (3x6), calcular G por forca bruta sistematica
% G deve satisfazer: para cada linha g de G, H * g' = 0 (mod 2)
% E G deve ter posto k = 3 (linhas linearmente independentes)

% Algoritmo: encontrar todas as solucoes de H*x = 0 e escolher k independentes
% Para codigo (6,3), o espaco nulo tem dimensao 3

% Metodo: eliminacao gaussiana em H para identificar variaveis livres
% Vamos colocar H na forma escalonada reduzida

H_rref = mod(H, 2);  % Garantir que e binaria

% Eliminacao gaussiana para forma escalonada
pivot_cols = [];
row = 1;
for col = 1:n
    % Procurar pivot
    pivot_row = 0;
    for r = row:m
        if H_rref(r, col) == 1
            pivot_row = r;
            break;
        end
    end
    
    if pivot_row > 0
        % Trocar linhas se necessario
        if pivot_row ~= row
            temp = H_rref(row, :);
            H_rref(row, :) = H_rref(pivot_row, :);
            H_rref(pivot_row, :) = temp;
        end
        
        % Eliminar abaixo e acima
        for r = 1:m
            if r ~= row && H_rref(r, col) == 1
                H_rref(r, :) = mod(H_rref(r, :) + H_rref(row, :), 2);
            end
        end
        
        pivot_cols = [pivot_cols, col];
        row = row + 1;
        if row > m
            break;
        end
    end
end

% Colunas livres (nao sao pivots)
free_cols = setdiff(1:n, pivot_cols);

% Construir G: para cada variavel livre, criar uma solucao basica
% G tem k linhas (uma para cada variavel livre) e n colunas
G = zeros(k, n);

for i = 1:k
    % i-esima linha de G: variavel livre i = 1, outras livres = 0
    x = zeros(1, n);
    x(free_cols(i)) = 1;  % Variavel livre
    
    % Calcular variaveis de pivot
    for j = 1:length(pivot_cols)
        pc = pivot_cols(j);
        % x(pc) = -sum(H_rref(j, free_cols) * x(free_cols)) mod 2
        x(pc) = mod(-H_rref(j, :) * x', 2);
    end
    
    G(i, :) = x;
end

% Debug: mostrar matrizes calculadas
fprintf('\n=== DEBUG Matrizes LDPC ===\n');
fprintf('H original:\n');
disp(H);
fprintf('H_rref (escalonada):\n');
disp(H_rref);
fprintf('Colunas pivot: '); disp(pivot_cols);
fprintf('Colunas livres: '); disp(free_cols);
fprintf('Matriz G calculada:\n');
disp(G);

% Salvar free_cols como variavel global para uso na decodificacao
free_cols_info = free_cols;

% Verificacao final
syndrome_check = mod(H * G', 2);
fprintf('H * G^T (deve ser zero):\n');
disp(syndrome_check);

if any(syndrome_check(:))
    error('Matrizes H e G inconsistentes! H*G^T != 0');
end

fprintf('Verificacao LDPC: H*G^T = 0 (mod 2) - OK\n');
fprintf('Matriz H (%dx%d): codigo (%d,%d) com taxa R=%.2f\n\n', m, n, n, k, k/n);

% Parametros do codigo LDPC
n_ldpc = 6;          % Comprimento do codeword
k_ldpc = 3;          % Bits de informacao
R_ldpc = k_ldpc/n_ldpc;  % Taxa do codigo = 1/2
max_iter = 30;       % Maximo de iteracoes do decodificador BP

% Parametros OFDM
N_FFT = 64;          % Tamanho da FFT
N_data = 48;         % Numero de subportadoras de dados (exclui DC e bordas)
cp_len = 16;         % Comprimento do prefixo ciclico (1/4 do simbolo)

% Parametros de modulacao
M = 16;              % 16-QAM
m_bits = log2(M);    % 4 bits por simbolo QAM

% Parametros de simulacao
SNR_dB = 0:2:14;     % Valores de SNR em dB
N_sym = 1000;        % Numero de simbolos OFDM por ponto de SNR

%% ============================================================
% CONFIGURACAO DA CONSTELACAO 16-QAM
%% ============================================================

% Amplitudes normalizadas para energia media unitaria
% 16-QAM com amplitudes +/-1, +/-3 em cada dimensao
amp = [-3, -1, 1, 3] / sqrt(10);  % Normalizacao para P_avg = 1

% Constelacao completa 16-QAM com mapeamento natural
% Bits (b1,b2) -> I, Bits (b3,b4) -> Q
constellation = zeros(1, M);
for idx = 0:M-1
    b1 = bitget(idx, 4);
    b2 = bitget(idx, 3);
    b3 = bitget(idx, 2);
    b4 = bitget(idx, 1);
    
    % Mapeamento: 00->-3, 01->-1, 10->+1, 11->+3
    % bi2de manual: valor = b1*2 + b2 (left-msb)
    I = amp(b1*2 + b2 + 1);
    Q = amp(b3*2 + b4 + 1);
    
    constellation(idx + 1) = I + 1j*Q;
end

% Tabela de bits para cada ponto da constelacao (para demodulacao)
bit_table = zeros(M, m_bits);
for i = 1:M
    bit_table(i, :) = bitget(i-1, m_bits:-1:1);
end

fprintf('Constelacao 16-QAM configurada com energia media = %.4f\n', ...
    mean(abs(constellation).^2));

%% ============================================================
% DIMENSIONAMENTO DO SISTEMA
%% ============================================================

% Bits por simbolo OFDM (subportadoras de dados * bits por QAM)
bits_per_ofdm = N_data * m_bits;

% Numero de blocos LDPC que cabem em um simbolo OFDM
% Cada bloco LDPC produz n_ldpc bits codificados
n_blocks_per_ofdm = floor(bits_per_ofdm / n_ldpc);

% Ajuste para garantir compatibilidade
bits_per_ofdm_ldpc = n_blocks_per_ofdm * n_ldpc;
info_bits_per_ofdm = n_blocks_per_ofdm * k_ldpc;

fprintf('\n=== Dimensionamento do Sistema ===\n');
fprintf('Bits por simbolo OFDM (16-QAM, %d subportadoras): %d\n', N_data, bits_per_ofdm);
fprintf('Blocos LDPC (n=%d) por simbolo OFDM: %d\n', n_ldpc, n_blocks_per_ofdm);
fprintf('Bits de informacao por simbolo (LDPC): %d\n', info_bits_per_ofdm);
fprintf('Bits codificados por simbolo (LDPC): %d\n', bits_per_ofdm_ldpc);
fprintf('Taxa efetiva: %.2f bits/simbolo QAM\n', info_bits_per_ofdm / N_data);
fprintf('Codigo LDPC: (%d,%d), taxa R = %.2f\n\n', n_ldpc, k_ldpc, R_ldpc);

%% ============================================================
% INICIALIZACAO DOS RESULTADOS
%% ============================================================

BER_uncoded = zeros(size(SNR_dB));
BER_ldpc = zeros(size(SNR_dB));

% Armazenar constelacoes para plot (SNR medio = 8 dB)
snr_const_idx = find(SNR_dB == 8);
if isempty(snr_const_idx)
    snr_const_idx = ceil(length(SNR_dB)/2);
end

tx_constellation = [];
rx_uncoded_const = [];
rx_ldpc_const = [];

%% ============================================================
% LOOP DE SIMULACAO
%% ============================================================

fprintf('Iniciando simulacao...\n');
fprintf('SNR (dB) | BER sem FEC | BER hard-LDPC| BER BP-LDPC  | Ganho (dB)\n');
fprintf('---------|-------------|--------------|--------------|------------\n');

for idx_snr = 1:length(SNR_dB)
    snr_db = SNR_dB(idx_snr);
    snr_lin = 10^(snr_db/10);
    
    % ========================================================
    % ETAPA 1: SISTEMA SEM FEC (referencia)
    % ========================================================
    
    num_err_uncoded = 0;
    total_bits_uncoded = 0;
    
    for sym = 1:N_sym
        % Gerar bits aleatorios
        bits = randi([0 1], 1, bits_per_ofdm);
        
        % Modulacao 16-QAM
        qam_symbols = qam16_mod(bits, constellation, m_bits);
        
        % OFDM modulacao
        ofdm_tx = ofdm_mod(qam_symbols, N_FFT, N_data, cp_len);
        
        % Canal AWGN - calculo da variancia do ruido
        % Es/N0 = SNR por simbolo. Para OFDM, cada subportadora eh um simbolo QAM
        % Potencia media do sinal OFDM deve ser normalizada
        signal_power = mean(abs(ofdm_tx).^2);
        noise_var = signal_power / snr_lin;
        noise_std = sqrt(noise_var/2);  % Divisao por 2 para componente complexa
        
        % Adicionar ruido AWGN
        noise = noise_std * (randn(size(ofdm_tx)) + 1j*randn(size(ofdm_tx)));
        ofdm_rx = ofdm_tx + noise;
        
        % OFDM demodulacao
        qam_rx = ofdm_demod(ofdm_rx, N_FFT, N_data, cp_len);
        
        % Demodulacao 16-QAM (hard decision)
        bits_rx = qam16_demod(qam_rx, constellation, bit_table, m_bits);
        
        % Contagem de erros
        num_err_uncoded = num_err_uncoded + sum(bits ~= bits_rx);
        total_bits_uncoded = total_bits_uncoded + length(bits);
        
        % Salvar constelacao para plot (apenas para ultimo simbolo no SNR de referencia)
        if idx_snr == snr_const_idx && sym == N_sym
            tx_constellation = qam_symbols;
            rx_uncoded_const = qam_rx;
        end
    end
    
    BER_uncoded(idx_snr) = num_err_uncoded / total_bits_uncoded;
    
    % ========================================================
    % ETAPA 2: SISTEMA COM LDPC
    % ========================================================
    
    num_err_ldpc = 0;
    num_err_hard = 0;  % DEBUG: para hard-decision sem BP
    total_info_bits = 0;
    
    for sym = 1:N_sym
        % Gerar bits de informacao
        info_bits = randi([0 1], 1, info_bits_per_ofdm);
        
        % Codificacao LDPC: processar em blocos de k_ldpc bits
        coded_bits = zeros(1, bits_per_ofdm_ldpc);
        for blk = 1:n_blocks_per_ofdm
            u = info_bits((blk-1)*k_ldpc + 1 : blk*k_ldpc);
            c = mod(u * G, 2);  % Codificacao matricial
            coded_bits((blk-1)*n_ldpc + 1 : blk*n_ldpc) = c;
        end
        
        % Modulacao 16-QAM
        qam_symbols = qam16_mod(coded_bits, constellation, m_bits);
        
        % Padding para preencher todas as subportadoras (se necessario)
        if length(qam_symbols) < N_data
            qam_symbols = [qam_symbols, zeros(1, N_data - length(qam_symbols))];
        end
        
        % OFDM modulacao
        ofdm_tx = ofdm_mod(qam_symbols, N_FFT, N_data, cp_len);
        
        % Canal AWGN (mesmo SNR que sistema sem FEC)
        signal_power = mean(abs(ofdm_tx).^2);
        noise_var = signal_power / snr_lin;
        noise_std = sqrt(noise_var/2);
        
        noise = noise_std * (randn(size(ofdm_tx)) + 1j*randn(size(ofdm_tx)));
        ofdm_rx = ofdm_tx + noise;
        
        % OFDM demodulacao
        qam_rx = ofdm_demod(ofdm_rx, N_FFT, N_data, cp_len);
        
        % Extrair apenas os simbolos correspondentes aos bits codificados
        qam_rx = qam_rx(1:(bits_per_ofdm_ldpc/m_bits));
        
        % Calcular LLRs para decodificacao LDPC
        % LLR = log(P(bit=0|y)/P(bit=1|y))
        % A variancia do ruido por dimensao (real/imag) apos FFT
        % No tempo: noise_var = signal_power / snr_lin
        % Na frequencia (apos FFT normalizada): mesma variancia por dimensao
        sigma2 = noise_var / 2;  % Variancia por dimensao real (I ou Q)
        llr_ch = compute_llr_16qam(qam_rx, constellation, bit_table, sigma2);
        
        % Decodificacao LDPC (iterativa - Belief Propagation)
        decoded_bits = zeros(1, info_bits_per_ofdm);
        
        for blk = 1:n_blocks_per_ofdm
            llr_block = llr_ch((blk-1)*n_ldpc + 1 : blk*n_ldpc);
            
            % Decodificador BP (Sum-Product)
            [decoded_block, ~] = ldpc_decode_bp(llr_block, H, max_iter);
            
            % Extrair bits de informacao
            % As colunas livres (free_cols_info) correspondem aos bits de info
            decoded_bits((blk-1)*k_ldpc + 1 : blk*k_ldpc) = decoded_block(free_cols_info);
        end
        
        % Contagem de erros (apenas nos bits de informacao)
        num_err_ldpc = num_err_ldpc + sum(info_bits ~= decoded_bits);
        total_info_bits = total_info_bits + length(info_bits);
        
        % DEBUG: Contagem hard-decision sem BP (comparacao)
        hard_bits = zeros(1, info_bits_per_ofdm);
        for blk = 1:n_blocks_per_ofdm
            llr_block = llr_ch((blk-1)*n_ldpc + 1 : blk*n_ldpc);
            % Hard decision direto dos LLRs: 1 se LLR < 0, 0 se LLR > 0
            hard_block = (llr_block < 0)';
            hard_bits((blk-1)*k_ldpc + 1 : blk*k_ldpc) = hard_block(free_cols_info);
        end
        err_hard = sum(info_bits ~= hard_bits);
        num_err_hard = num_err_hard + err_hard;
        
        % Salvar constelacao para plot
        if idx_snr == snr_const_idx && sym == N_sym
            rx_ldpc_const = qam_rx(1:min(length(qam_rx), N_data));
        end
    end
    
    BER_ldpc(idx_snr) = num_err_ldpc / total_info_bits;
    
    % Calcular ganho aproximado (para exibicao)
    % Ganho simples: melhoria multiplicativa na BER
    if BER_uncoded(idx_snr) > 0 && BER_ldpc(idx_snr) > 0
        ganho_factor = BER_uncoded(idx_snr) / BER_ldpc(idx_snr);
        ganho = 10*log10(ganho_factor);  % Em dB aproximado
    else
        ganho = 0;
    end
    
    ber_hard_debug = num_err_hard / total_info_bits;
    fprintf('%7.1f  |  %.4e  |   %.4e   |   %.4e  |  %6.2f\n', ...
        snr_db, BER_uncoded(idx_snr), ber_hard_debug, BER_ldpc(idx_snr), ganho);
end

fprintf('\nSimulacao concluida!\n');

%% ============================================================
% GRAFICOS E RESULTADOS
%% ============================================================

% Criar diretorio para figuras
if ~exist('figuras', 'dir')
    mkdir('figuras');
end

% 1. Grafico de BER versus SNR
figure('Position', [100 100 800 600]);
semilogy(SNR_dB, BER_uncoded, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName', 'Sem FEC (16-QAM)');
hold on;
semilogy(SNR_dB, BER_ldpc, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName', 'Com LDPC (R=1/2)');
grid on;
xlabel('E_s/N_0 (dB)', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title('Curvas BER versus SNR - Sistema QAM-OFDM', 'FontSize', 14);
legend('Location', 'southwest', 'FontSize', 10);
axis([min(SNR_dB) max(SNR_dB) 1e-4 1]);

% Salvar figura
saveas(gcf, 'figuras/ber_vs_snr.pdf');
saveas(gcf, 'figuras/ber_vs_snr.png');

% 2. Constelacoes QAM
figure('Position', [100 100 1200 400]);

% Transmissao
subplot(1, 3, 1);
plot(real(tx_constellation), imag(tx_constellation), 'bo', 'MarkerSize', 6);
grid on;
xlabel('Parte Real', 'FontSize', 10);
ylabel('Parte Imaginaria', 'FontSize', 10);
title('Constelacao 16-QAM - Transmissao', 'FontSize', 12);
axis([-1 1 -1 1]);
axis square;

% Recepcao sem FEC
subplot(1, 3, 2);
plot(real(rx_uncoded_const), imag(rx_uncoded_const), 'r.', 'MarkerSize', 4);
hold on;
plot(real(constellation), imag(constellation), 'k+', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
xlabel('Parte Real', 'FontSize', 10);
ylabel('Parte Imaginaria', 'FontSize', 10);
title(sprintf('Recepcao sem FEC (SNR = %d dB)', SNR_dB(snr_const_idx)), 'FontSize', 12);
axis([-1 1 -1 1]);
axis square;
legend('Simbolos recebidos', 'Simbolos ideais', 'Location', 'best');

% Recepcao com LDPC
subplot(1, 3, 3);
plot(real(rx_ldpc_const), imag(rx_ldpc_const), 'g.', 'MarkerSize', 4);
hold on;
plot(real(constellation), imag(constellation), 'k+', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
xlabel('Parte Real', 'FontSize', 10);
ylabel('Parte Imaginaria', 'FontSize', 10);
title(sprintf('Recepcao com LDPC (SNR = %d dB)', SNR_dB(snr_const_idx)), 'FontSize', 12);
axis([-1 1 -1 1]);
axis square;
legend('Simbolos recebidos', 'Simbolos ideais', 'Location', 'best');

saveas(gcf, 'figuras/constelacoes_qam.pdf');
saveas(gcf, 'figuras/constelacoes_qam.png');

% 3. Espectro do sinal OFDM
figure('Position', [100 100 800 400]);

% Gerar sinal OFDM para analise espectral
num_spectral_sym = 100;
ofdm_signal = [];
for i = 1:num_spectral_sym
    bits = randi([0 1], 1, bits_per_ofdm);
    qam_sym = qam16_mod(bits, constellation, m_bits);
    ofdm_frame = ofdm_mod(qam_sym, N_FFT, N_data, cp_len);
    ofdm_signal = [ofdm_signal, ofdm_frame];
end

% Calcular PSD usando FFT (sem dependencia de toolbox)
N_psd = 4096;
% Janela retangular (ou implementar hanning manualmente se desejado)
L = length(ofdm_signal);
w = 0.5 * (1 - cos(2*pi*(0:L-1)/(L-1)));  % Janela Hanning manual
X = fft(ofdm_signal .* w, N_psd);
pxx = abs(X).^2 / (sum(w.^2));
pxx = fftshift(pxx);  % Centralizar
f = (-N_psd/2:N_psd/2-1) / N_psd;  % Frequencias normalizadas

plot(f, 10*log10(pxx), 'b', 'LineWidth', 1.5);
grid on;
xlabel('Frequencia normalizada (f/f_s)', 'FontSize', 12);
ylabel('Densidade espectral de potencia (dB)', 'FontSize', 12);
title('Espectro aproximado do sinal OFDM', 'FontSize', 14);
xlim([-0.5 0.5]);
ylim([-80 20]);

% Marcar banda util (subportadoras de dados)
banda_util = N_data/N_FFT;
hold on;
xline(-banda_util/2, 'r--', 'Banda util');
xline(banda_util/2, 'r--');

saveas(gcf, 'figuras/espectro_ofdm.pdf');
saveas(gcf, 'figuras/espectro_ofdm.png');

%% ============================================================
% TABELA COMPARATIVA
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('TABELA COMPARATIVA: Sistema sem FEC vs. Sistema com LDPC\n');
fprintf('============================================================\n');
fprintf('SNR (dB) | BER sem FEC | BER com LDPC | Melhoria (x vezes)\n');
fprintf('---------|-------------|--------------|--------------------\n');

for i = 1:length(SNR_dB)
    if BER_ldpc(i) > 0
        melhoria = BER_uncoded(i) / BER_ldpc(i);
    else
        melhoria = Inf;
    end
    fprintf('%7.1f  |  %.4e  |   %.4e   |    %8.2f\n', ...
        SNR_dB(i), BER_uncoded(i), BER_ldpc(i), melhoria);
end

fprintf('============================================================\n');

% Salvar resultados em arquivo
fid = fopen('resultados_ber.txt', 'w');
fprintf(fid, 'Resultados da Simulacao TA2\n');
fprintf(fid, 'LDPC (%d,%d), 16-QAM, OFDM com %d subportadoras\n\n', ...
    n_ldpc, k_ldpc, N_data);
fprintf(fid, 'SNR_dB = [%s]\n', num2str(SNR_dB));
fprintf(fid, 'BER_sem_FEC = [%s]\n', num2str(BER_uncoded, '%.6e '));
fprintf(fid, 'BER_com_LDPC = [%s]\n', num2str(BER_ldpc, '%.6e '));
fclose(fid);

fprintf('\nResultados salvos em:\n');
fprintf('  - figuras/ber_vs_snr.pdf (e .png)\n');
fprintf('  - figuras/constelacoes_qam.pdf (e .png)\n');
fprintf('  - figuras/espectro_ofdm.pdf (e .png)\n');
fprintf('  - resultados_ber.txt\n');

end

%% ============================================================
% FUNCOES AUXILIARES
%% ============================================================

function symbols = qam16_mod(bits, constellation, m_bits)
% Modulacao 16-QAM
% bits: vetor de bits (0/1)
% constellation: pontos da constelacao
% m_bits: bits por simbolo (4 para 16-QAM)

    num_symbols = floor(length(bits) / m_bits);
    symbols = zeros(1, num_symbols);
    
    for i = 1:num_symbols
        bit_chunk = bits((i-1)*m_bits + 1 : i*m_bits);
        % bi2de manual (left-msb): converte vetor de bits para decimal
        idx = 1;
        for b = 1:m_bits
            idx = idx + bit_chunk(b) * 2^(m_bits - b);
        end
        symbols(i) = constellation(idx);
    end
end

function bits = qam16_demod(symbols, constellation, bit_table, m_bits)
% Demodulacao 16-QAM (hard decision)
% symbols: simbolos recebidos
% constellation: pontos da constelacao
% bit_table: tabela de mapeamento bit->simbolo

    num_symbols = length(symbols);
    bits = zeros(1, num_symbols * m_bits);
    
    for i = 1:num_symbols
        % Encontrar ponto mais proximo
        dist = abs(symbols(i) - constellation);
        [~, idx] = min(dist);
        
        % Mapear para bits
        bits((i-1)*m_bits + 1 : i*m_bits) = bit_table(idx, :);
    end
end

function llr = compute_llr_16qam(symbols, constellation, bit_table, sigma2)
% Calcular LLRs para 16-QAM
% LLR(b) = log(P(b=0|y)/P(b=1|y))
% Aproximacao: log-sum-exp para estabilidade numerica

    m_bits = size(bit_table, 2);
    num_symbols = length(symbols);
    llr = zeros(num_symbols * m_bits, 1);
    
    for i = 1:num_symbols
        y = symbols(i);
        
        for b = 1:m_bits
            % Calcular probabilidades para bit b = 0 e b = 1
            idx_0 = find(bit_table(:, b) == 0);
            idx_1 = find(bit_table(:, b) == 1);
            
            % P(b=0|y) proporcional a sum_{s in S_0} exp(-|y-s|^2/(2*sigma2))
            d2_0 = abs(y - constellation(idx_0)).^2;
            d2_1 = abs(y - constellation(idx_1)).^2;
            
            % LLR = log(sum(exp(-d2_0/2sigma))) - log(sum(exp(-d2_1/2sigma)))
            llr_val = logsumexp(-d2_0/(2*sigma2)) - logsumexp(-d2_1/(2*sigma2));
            
            llr((i-1)*m_bits + b) = llr_val;
        end
    end
end

function s = logsumexp(x)
% Calculo estavel de log(sum(exp(x)))
    m = max(x);
    s = m + log(sum(exp(x - m)));
end

function [decoded, converged] = ldpc_decode_bp(llr_ch, H, max_iter)
% Decodificador LDPC usando Belief Propagation (Sum-Product)
% llr_ch: LLRs do canal (vetor coluna)
% H: matriz de paridade
% max_iter: maximo de iteracoes

    [m, n] = size(H);
    
    % Inicializar mensagens variavel -> check
    L_v2c = repmat(llr_ch, 1, m) .* H';
    
    % Iteracoes BP
    for iter = 1:max_iter
        % Check node update: mensagens check -> variavel
        L_c2v = zeros(n, m);
        
        for j = 1:m
            % Variaveis conectadas ao check node j
            var_idx = find(H(j, :) == 1);
            
            for i = var_idx
                % Produto de tanh das mensagens de outras variaveis
                other_vars = setdiff(var_idx, i);
                if ~isempty(other_vars)
                    tanh_prod = prod(tanh(L_v2c(other_vars, j) / 2));
                    % Limitar para estabilidade numerica
                    tanh_prod = max(min(tanh_prod, 0.9999), -0.9999);
                    L_c2v(i, j) = 2 * atanh(tanh_prod);
                end
            end
        end
        
        % Variavel node update e decisao
        L_total = llr_ch + sum(L_c2v, 2);
        
        % Decisao hard
        decoded = (L_total < 0)';
        
        % Verificar se atingiu codeword valido (sindrome = 0)
        syndrome = mod(H * decoded', 2);
        if all(syndrome == 0)
            converged = true;
            return;
        end
        
        % Atualizar mensagens variavel -> check para proxima iteracao
        for i = 1:n
            check_idx = find(H(:, i) == 1);
            for j = check_idx'
                L_v2c(i, j) = llr_ch(i) + sum(L_c2v(i, setdiff(check_idx, j)));
            end
        end
    end
    
    converged = false;
end

function ofdm_signal = ofdm_mod(qam_symbols, N_FFT, N_data, cp_len)
% Modulacao OFDM
% qam_symbols: simbolos QAM (vetor linha)
% N_FFT: tamanho da FFT
% N_data: numero de subportadoras de dados
% cp_len: comprimento do prefixo ciclico

    % Montar vetor de frequencias (subportadoras centrais)
    ofdm_frame = zeros(1, N_FFT);
    
    % Posicionar simbolos nas subportadoras de dados (excluindo DC e bordas)
    start_idx = floor((N_FFT - N_data) / 2) + 1;
    end_idx = start_idx + N_data - 1;
    
    num_sym = min(length(qam_symbols), N_data);
    ofdm_frame(start_idx:start_idx+num_sym-1) = qam_symbols(1:num_sym);
    
    % IFFT (normalizada)
    time_signal = ifft(ofdm_frame) * sqrt(N_FFT);
    
    % Adicionar prefixo ciclico
    ofdm_signal = [time_signal(end-cp_len+1:end), time_signal];
end

function qam_symbols = ofdm_demod(ofdm_signal, N_FFT, N_data, cp_len)
% Demodulacao OFDM
% ofdm_signal: sinal OFDM recebido (com CP)
% N_FFT: tamanho da FFT
% N_data: numero de subportadoras de dados
% cp_len: comprimento do prefixo ciclico

    % Remover prefixo ciclico
    time_signal = ofdm_signal(cp_len+1:end);
    
    % FFT (normalizada)
    freq_signal = fft(time_signal) / sqrt(N_FFT);
    
    % Extrair subportadoras de dados
    start_idx = floor((N_FFT - N_data) / 2) + 1;
    end_idx = start_idx + N_data - 1;
    
    qam_symbols = freq_signal(start_idx:end_idx);
end
