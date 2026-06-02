# TA2 2026/1 - Sistemas de Comunicação
## Projeto Prático MATLAB: LDPC com 16-QAM e OFDM

### Autor
**Bruno Becker Silva**  
Escola Politécnica - PUCRS  
Disciplina: 445AH - Sistemas de Comunicação  
Professor: Flavio Eduardo Soares e Silva

---

## Descrição do Projeto

Implementação completa de uma cadeia digital de comunicação integrando:
- **Codificação LDPC** (código (6,3) didático sem toolbox)
- **Modulação 16-QAM**
- **OFDM** com 64 subportadoras (48 de dados)
- **Canal AWGN**
- **Decodificação LDPC** via Belief Propagation

Comparação entre sistema **sem FEC** e sistema **com LDPC**.

---

## Arquivos Entregues

| Arquivo | Descrição |
|---------|-----------|
| `ta2_ldpc_qam_ofdm.m` | Código MATLAB comentado com simulação completa |
| `ta2_relatorio.tex` | Relatório em LaTeX (IEEEtran) |
| `figuras/` | Diretório onde as figuras são salvas |
| `resultados_ber.txt` | Arquivo com resultados numéricos |

---

## Requisitos Mínimos Atendidos

- [x] Modulação 16-QAM
- [x] Sistema OFDM com 64 subportadoras (48 de dados)
- [x] Canal AWGN
- [x] Cálculo de BER
- [x] Codificação LDPC (matriz H 3×6 do Capítulo 7)
- [x] Decodificação LDPC (Belief Propagation)
- [x] Comparação sem codificação vs. com LDPC

---

## Execução

### No MATLAB:

```matlab
% Navegar até o diretório do projeto
cd /Users/brunobks/Faculdade/comunicacao

% Executar a simulação
ta2_ldpc_qam_ofdm()
```

### Saída Esperada:

1. **Verificação LDPC:** Confirmação de que H*G^T = 0 (mod 2)
2. **Dimensionamento:** Informações sobre bits por símbolo OFDM
3. **Progresso da simulação:** Tabela com BER em cada ponto de SNR
4. **Figuras geradas:**
   - `figuras/ber_vs_snr.png` - Curvas BER versus SNR
   - `figuras/constelacoes_qam.png` - Constelações antes/depois do canal
   - `figuras/espectro_ofdm.png` - Espectro do sinal OFDM
5. **Resultados numéricos:** Arquivo `resultados_ber.txt`

---

## Compilação do Relatório LaTeX

```bash
# Usando pdflatex
pdflatex ta2_relatorio.tex
pdflatex ta2_relatorio.tex  # Segunda passada para referências

# Ou usar latexmk
latexmk -pdf ta2_relatorio.tex
```

**Observação:** Certifique-se de que as figuras em `figuras/` foram geradas pelo MATLAB antes de compilar o relatório.

---

## Estrutura do Código MATLAB

### Funções Principais:

1. **`ta2_ldpc_qam_ofdm()`** - Função principal com simulação
2. **`qam16_mod()`** - Modulador 16-QAM
3. **`qam16_demod()`** - Demodulador 16-QAM (hard decision)
4. **`compute_llr_16qam()`** - Cálculo de LLRs para decodificação
5. **`ldpc_decode_bp()`** - Decodificador Belief Propagation
6. **`ofdm_mod()`** - Modulador OFDM (IFFT + CP)
7. **`ofdm_demod()`** - Demodulador OFDM (remoção CP + FFT)
8. **`logsumexp()`** - Função auxiliar para estabilidade numérica

### Parâmetros Configuráveis:

```matlab
% No início do arquivo
N_sym = 1000;        % Número de símbolos OFDM por ponto de SNR
max_iter = 30;       % Máximo de iterações do decodificador BP
SNR_dB = 0:2:14;     % Pontos de SNR para simulação
```

---

## Respostas às Perguntas do Relatório

As 10 perguntas da tarefa são respondidas no relatório LaTeX:

1. **Diferença detecção/correção:** Seção II-A
2. **BER e PER:** Seção II-B
3. **Função da matriz H:** Seção II-C
4. **Significado da síndrome:** Seção II-D
5. **Grafo de Tanner:** Seção II-E
6. **Decisão rígida/suave:** Seção II-F
7. **Prefixo cíclico no OFDM:** Seção II-G
8. **SNR na cadeia:** Seção V-A
9. **Melhoria da BER com LDPC:** Seção V-B
10. **Custo da melhoria:** Seção V-C

---

## Notas de Implementação

### Matriz LDPC (Trilha 1):

O código utiliza a matriz didática fornecida:
```
H = [1 1 1 1 0 0
     0 1 0 1 0 1
     1 1 0 0 0 1]
```

Esta é uma matriz pequena (3×6) que permite:
- Visualização didática do processo
- Simulação rápida
- Demonstração dos conceitos sem dependência de toolboxes

### Limitações do Código Didático:

- Código curto (n=6) limita o ganho de codificação
- Para aplicações práticas, códigos LDPC com n > 1000 são recomendados
- A curva de BER pode apresentar variância estatística devido ao tamanho limitado

---

## Referências

- Gallager, R. G. (1962). Low-density parity-check codes.
- MacKay, D. J. C. (1999). Good error-correcting codes based on very sparse matrices.
- Richardson & Urbanke (2008). Modern Coding Theory.
- Proakis & Salehi (2008). Digital Communications.

---

## Contato

Dúvidas sobre o código ou relatório: consultar o professor da disciplina.
