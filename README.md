# Qualidade do Ar: Análise e Engenharia de Features

Primeiro trabalho de **Redes Neuronais Artificiais 2026/2027**, Faculdade de Ciências Exatas e da Engenharia, Universidade da Madeira.

## Objetivo

Explorar, limpar, transformar e criar features a partir do dataset UCI Air Quality, usando o **CO(GT)** (concentração horária de monóxido de carbono) como variável alvo e as restantes variáveis como features de entrada. O foco está em compreender e documentar cada passo: pré-processamento, transformação, seleção e engenharia de features.

## Grupo

|      Nome       | Número de aluno |
|-----------------|-----------------|
| Rodrigo Soares  |     2075526     |
| Leonardo Caneco |     2048726     |


## Dataset

- **Fonte:** [UCI Machine Learning Repository: Air Quality](https://archive.ics.uci.edu/dataset/360/air+quality)
- **Conteúdo:** médias horárias das leituras de um dispositivo multissensor e de um analisador de referência, instalados numa cidade italiana (9357 horas, de março de 2004 a abril de 2005).
- **Variável alvo:** `CO(GT)`
- **Entradas:** as restantes variáveis de sensores, poluentes e meteorologia, mais a data e a hora.

O CSV original fica em `data/Raw/` e nunca é alterado.

## Estrutura do repositório

```
├── README.md
├── aqlib.m              Function library (um só ficheiro, struct de function handles)
├── main.m               Script do projeto: executa os passos 1 a 5 usando a aqlib
├── data/
│   ├── Raw/        AirQualityUCI.csv original
│   └── Processed/  Dados limpos e features criadas
├── Figures/             Gráficos exportados para o relatório e os slides
```

Ficheiros em `data/Processed/` (todos gerados pelo `main.m`):

| Ficheiro | Conteúdo |
|---|---|
| `airquality_clean.csv` | Dados limpos (passo 2), 7344 horas, unidades originais |
| `airquality_scaled.csv` | O mesmo, entradas em z-score |
| `scaling_params.csv` | Média e desvio padrão usados no z-score do passo 2 |
| `airquality_selected.csv` | As 7 entradas selecionadas (passo 4), unidades originais |
| `airquality_selected_scaled.csv` | O mesmo, entradas em z-score |
| `airquality_engineered.csv` | Dataset final (passo 5): 7 entradas + 3 features novas, 7329 horas, unidades originais |
| `airquality_engineered_scaled.csv` | O mesmo, entradas em z-score (o alvo fica em mg/m³) |
| `scaling_params_engineered.csv` | Média e desvio padrão de cada entrada do dataset final, para escalar dados novos |

## Requisitos

- MATLAB R2026b (versão usada no desenvolvimento e nos testes)
- Statistics and Machine Learning Toolbox

## Como executar

Abrir a pasta do repositório no MATLAB como pasta atual e correr o `main.m`. Cada secção do script corresponde a um passo do trabalho e pode ser executada isoladamente (Run Section).

## Passos do trabalho

- [x] **1. Importação e exploração dos dados:** descrição do dataset, valores em falta, distribuição da variável alvo, estatísticas de cada variável, correlação com a variável alvo, tabela resumo.
- [x] **2. Limpeza e pré-processamento:** identificação de valores em falta e anómalos, imputação, clipping, escalamento, comparações antes/depois.
- [x] **3. Análise exploratória (EDA):** gráficos que relacionam as entradas com a variável alvo, assimetria, outliers, relações lineares, redundância.
- [x] **4. Seleção de features:** relevância (Pearson, Spearman, informação mútua), redundância (agrupamento hierárquico e PCA), informação em contexto (RReliefF contra uma coluna de ruído) e estabilidade no tempo. Subconjunto final: 7 de 11 entradas (PT08.S1, PT08.S2, PT08.S4, PT08.S5, T, RH, AH).
- [x] **5. Engenharia de features:** 4 features novas calculadas numa grelha horária (hora do dia, média de 3 h, desvio à média de 24 h e variação horária do PT08.S2), avaliadas com os critérios do passo 4. Ficam 3; o dataset final tem 10 entradas e 7329 horas.


## Referência

S. De Vito, E. Massera, M. Piga, L. Martinotto, G. Di Francia, "On field calibration of an electronic nose for benzene estimation in an urban pollution monitoring scenario," *Sensors and Actuators B: Chemical*, vol. 129, no. 2, pp. 750–757, 2008.