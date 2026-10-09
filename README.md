# Qualidade do Ar: Análise e Engenharia de Features

Primeiro trabalho de **Redes Neuronais Artificiais 2026/2027**, Faculdade de Ciências Exatas e da Engenharia, Universidade da Madeira.

## Objetivo

Explorar, limpar, transformar e criar features a partir do dataset UCI Air Quality, usando o **CO(GT)** (concentração horária de monóxido de carbono) como variável alvo e as restantes variáveis como features de entrada. O foco está em compreender e documentar cada passo: pré-processamento, transformação, seleção e engenharia de features.

## Grupo

|     Nome       | Número de aluno |
| Rodrigo Soares |     2075526     |
|                |                 |


## Dataset

- **Fonte:** [UCI Machine Learning Repository: Air Quality](https://archive.ics.uci.edu/dataset/360/air+quality)
- **Conteúdo:** médias horárias das leituras de um dispositivo multissensor e de um analisador de referência, instalados numa cidade italiana (março de 2004 a fevereiro de 2005).
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

## Requisitos

- MATLAB R2016b ou mais recente
- Statistics and Machine Learning Toolbox

## Como executar

Abrir a pasta do repositório no MATLAB como pasta atual e correr o `main.m`. Cada secção do script corresponde a um passo do trabalho e pode ser executada isoladamente (Run Section).

## Passos do trabalho

- [ ] **1. Importação e exploração dos dados:** descrição do dataset, valores em falta, distribuição da variável alvo, estatísticas de cada variável, correlação com a variável alvo, tabela resumo.
- [ ] **2. Limpeza e pré-processamento:** identificação de valores em falta e anómalos, imputação, clipping, escalamento, comparações antes/depois.
- [ ] **3. Análise exploratória (EDA):** gráficos que relacionam as entradas com a variável alvo, assimetria, outliers, relações lineares, redundância.
- [ ] **4. Seleção de features:** análise de redundância, critérios de seleção, subconjunto final de features justificado.
- [ ] **5. Engenharia de features:** pelo menos duas features novas, visualização, avaliação da sua utilidade.


## Referência

S. De Vito, E. Massera, M. Piga, L. Martinotto, G. Di Francia, "On field calibration of an electronic nose for benzene estimation in an urban pollution monitoring scenario," *Sensors and Actuators B: Chemical*, vol. 129, no. 2, pp. 750–757, 2008.