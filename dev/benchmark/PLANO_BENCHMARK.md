# Plano de benchmark do BLASTr (rascunho, 2026-11-18)

Objetivo para o artigo: determinar a configuração ótima de busca e mostrar **se** (e quando)
a busca paralelizada e dividida do BLASTr é melhor que um único `blastn` com todas as queries.

## 1. Ponto de partida

**Como o BLASTr busca hoje** (`R/parallel_blast.R` → `R/blast_cmd.R`):

- Faz um `blastn` por sequência única. Cada um é disparado por `condathis::run_bin()`
  dentro de um worker `mirai`.
- Parâmetros: `-max_hsps 1`, `-perc_identity 80`, `-qcov_hsp_perc 80`, `-num_alignments 4`
  e `-mt_mode 2`, com `total_cores` processos de `num_threads` threads cada.
- Cada processo abre e varre o banco inteiro: o core_nt tem 299 GB em 92 volumes.

**O que já foi medido numa corrida de COI do mtbrcdR** (mtbrcdR/CLAUDE.md §12, 2026-11-18), com o
banco maior que o cache de páginas livre:

- 30 workers a ~26% de CPU;
- I/O pressure de ~45% e ~240 MB/s de leitura;
- ~15 min por query de 66 pb, o que dá ~120 queries/h.

**Hipótese:** o custo é dominado pela varredura do banco, que se repete N vezes. O
`blastn` com várias queries concatena as queries em lotes internos (variável `BATCH_SIZE`)
e varre o banco **uma vez por lote**. Então lotes de queries deveriam amortizar o I/O.

**Servidor** (`nproc` etc., 2026-11-18):

| Item | Valor |
|---|---|
| CPU | 2× Xeon Gold 6248: 40 núcleos físicos, 80 threads, 2 nós NUMA |
| RAM | 503 GB, com ~416 GB em cache de páginas |
| Disco `/data` | RAID de discos rotacionais (PERC H730P), 47 TB |
| BLAST | sistema 2.12.0+; o BLASTr usa o env conda (conferir a versão; `-mt_mode` exige ≥ 2.12) |
| Ferramentas | `fincore` e `pidstat` disponíveis; cgroup v1 |
| Carga na hora da medição | load ~65, com 92 `blastn` rodando (de dois usuários) |

**Conclusão: o benchmark precisa de uma janela com o servidor exclusivo.**

## 2. Perguntas → experimentos

| # | Pergunta | Fatores | Resposta (gráfico) |
|---|---|---|---|
| Q1 | Lote vs. 1 query/processo vs. 1 `blastn` único | estratégia A/B/C | tempo total e bytes lidos por estratégia |
| Q2 | Como escala com núcleos? Existe ótimo P × T × tamanho de lote? | P, T, k | heatmap P×T por k; curva de speedup/eficiência |
| Q3 | Peso de 3 a 10 hits por query | `num_alignments` ∈ {1, 3, 5, 10, 25, 50} | tempo, CPU e RAM relativos a 1 hit |
| Q4 | CPU, RAM e disco | métricas (§4) | fronteira de Pareto tempo × CPU-h × RAM |
| Q5 | Como fica em outros computadores? | RAM/banco, núcleos, tipo de disco | tempo normalizado vs. (tamanho do banco / RAM disponível) |
| Q6 | O resultado é o mesmo em todas as configurações? | todas | % de hits idênticos (top-k, staxid, identidade, e-value) |
| Q7 | Efeito da composição das queries | tamanho, marcador, com/sem hit | tempo por query vs. tamanho, por marcador |

**Estratégias:**

- **A (atual):** 1 query por `blastn`; P processos × T threads.
- **B (único):** 1 `blastn` com todas as queries e T threads, testando `-mt_mode` 0 (auto),
  1 (threads dividem as queries) e 2 (threads dividem o banco).
- **C (híbrida, candidata):** as N queries divididas em lotes de tamanho k e distribuídas
  em P processos × T threads.

Para referência, também medir **D: `blastn` puro no shell**, sem R, mirai nem condathis.
Isso isola o overhead do pacote (o `micromamba run` custa por chamada) do custo do BLAST.

## 3. Conjunto de queries (dataset representativo)

**Fontes já localizadas** (só sequências e metadados técnicos; nenhum nome de amostra ou cliente):

| Fonte | Marcadores | Seqs | Observação |
|---|---|---|---|
| fixtures do golden do mtbrcdR (duas corridas) | 12S p12SV5inv, 16S p16SMam1, rbcL, trnL | 1.497 com hit (+ sem hit no rds) | hits do core_nt anterior já conhecidos |
| análise de COI do mtbrcdR (swarm) | COI (CoqFish), merged + concat | 15.826 merged (50–712 pb), 6.052 concat | concat com `N`×10 no meio |
| base de uma análise de plantas do mtbrcdR (sqlite) | ITS2, rbcL (plantas) | 2.529 | — |
| bases do curatR (`db/*.data.sqlite`) | quatro projetos (114 mil, 2.057, 1.658 e 300 seqs) | — | colunas de pseudo-score/ID: estratos com/sem hit |

Os caminhos e os nomes das fontes ficam numa tabela local, fora do git (`BLASTR_BENCH_SOURCES`;
as variáveis de ambiente dos scripts estão no início de `R/bench_lib.R`).

**Amostragem estratificada:**

- Estratos: marcador × faixa de tamanho (≤100, 101–200, 201–400, >400 pb) × classe de
  hit pela ID já conhecida (≥98% pseudo-score, 80–98%, sem hit) × origem (merged/concat).
- Pool de ~1.000 queries, de onde saem subconjuntos **aninhados**:
  - `S100`: 100 queries, balanceado, usado em Q1–Q4;
  - `S500` e `S2000`: usados no escalonamento em N.
- Um conjunto "experimento real" completo (uma corrida inteira, ~3 mil seqs), para validar a
  configuração escolhida de ponta a ponta.
- Congelar em `dev/benchmark/data/queries.fasta` + `queries_meta.tsv`
  (id, marcador, tamanho, classe, origem, fonte), com md5.

## 4. Métricas e instrumentação

Cada rodada gera uma linha em `results/runs.tsv` e uma série temporal em `results/ts/<run_id>.tsv`.

- **Tempo:** wall clock total, tempo por query (estratégia A) ou por lote (C).
- **CPU:** user + sys de todos os `blastn`, lidos de `/proc/<pid>/stat` ou com `pidstat -u`.
  Eficiência = CPU-s / (wall × núcleos alocados).
- **RAM:** PSS (`/proc/<pid>/smaps_rollup`), porque o banco é mapeado em memória (mmap) e o
  RSS conta as páginas compartilhadas várias vezes. Registrar também a RAM livre e o cache.
- **Disco:** `read_bytes` de `/proc/<pid>/io`; `sda` em `/proc/diskstats` (bytes e %util);
  I/O pressure em `/proc/pressure/io`, se o kernel 5.4 expuser.
- **Cache:** `fincore` nos volumes do banco antes e depois da rodada (% do banco em RAM).
- **Contexto:** host, CPU, RAM, disco, versões (BLAST, R, BLASTr), tamanho do banco, load
  no início e a ordem da rodada.

**Estado do cache é um fator, não ruído.** Cada configuração roda:

- **fria:** banco despejado do cache. Usar `vmtouch -e`/`posix_fadvise(DONTNEED)` nos
  arquivos do banco, ou `drop_caches` se houver sudo. **Isso afeta os outros BLAST do servidor.**
- **quente:** banco pré-carregado, com `cat`/`vmtouch -t`, se couber.

**Repetições:** 3 por configuração, em ordem aleatória (para não confundir a configuração
com o tempo ou o aquecimento do cache).

## 5. Fases (em ordem; cada uma decide o grid da próxima)

0. **Infra e piloto (banco local de 12S + 20 queries, minutos):**
   - harness e coletor de métricas;
   - checar que A, B e C dão os mesmos hits;
   - medir o overhead fixo por chamada (R + mirai + condathis) contra um banco minúsculo.
1. **Triagem de estratégia (core_nt, `S100`, 3 hits, fria e quente):**
   - A com P ∈ {8, 30, 60} e T = 1;
   - B com T ∈ {8, 30, 60} × `mt_mode` {0, 1, 2};
   - C com k ∈ {10, 25, 50} e P×T ≈ 60.
   Pergunta-chave: quantas vezes o banco é lido em cada estratégia?
2. **Superfície de resposta da melhor estratégia:**
   - P ∈ {1, 2, 4, 8, 16, 32}, T ∈ {1, 2, 4, 8, 16} (P×T ≤ 80, mais 1–2 pontos de sobrecarga);
   - k ∈ {1, 10, 25, 50, 100, N/P}.
   - Fixar NUMA (`numactl`) num subconjunto, para ver se importa.
3. **Hits por query:** `num_alignments` ∈ {1, 3, 5, 10, 25, 50} na melhor configuração e no
   modo atual. Inclui o caso concat (sem `-max_hsps`, ~50 alinhamentos).
4. **Escalonamento em N:** `S100` → `S500` → `S2000` → uma corrida completa. Gráfico de
   queries/h × N em log-log: a vantagem do lote deve crescer com N.
5. **Outras máquinas:**
   - simular no servidor: menos núcleos com `taskset`; menos RAM ocupando memória com um
     processo que segura páginas (o cgroup v1 não limita o cache de forma confiável);
     banco menor (o banco local de 12S, ou core_nt filtrado por táxon);
   - medir em 1–2 máquinas reais (notebook/workstation com SSD) com o mesmo `S100` e um banco
     que caiba e outro que não caiba na RAM.
   - Modelo para o artigo: tempo ≈ f(varreduras do banco × tamanho do banco / vazão do disco
     ou da RAM, CPU por query). Assim a recomendação vale além deste servidor.

## 6. Gráficos-alvo (rascunho das figuras)

1. Barras/pontos: tempo total e GB lidos por estratégia (A/B/C), faceta fria/quente.
2. Heatmap P × T (cor = wall time), uma faceta por k; ótimo marcado.
3. Speedup e eficiência paralela vs. núcleos totais (curva de Amdahl ajustada).
4. Queries/h vs. N (log-log), uma linha por estratégia.
5. Custo relativo vs. `num_alignments`, com tempo, CPU-h e RAM.
6. Fronteira de Pareto: wall time × CPU-h, com tamanho do ponto = pico de PSS.
7. Tempo por query vs. tamanho da query, colorido por marcador e com/sem hit.
8. Tempo normalizado vs. (banco / RAM disponível), uma cor por máquina.
9. Concordância de resultados entre estratégias (tabela ou tile).

## 7. Mudanças no pacote que o benchmark puxa

- **`parallel_blast(batch_size = )`** (estratégia C): FASTA com cabeçalhos `>q<i>`, para
  remapear para a sequência; cache por lote; retry por lote. Implementar o protótipo
  **depois** da fase 1 confirmar o ganho, para que o artigo meça "pelo pacote".
- **Corrigir o `run_blast()`:** hoje ele grava um vetor de sequências com `cat()`, o que
  vira **uma única query concatenada por espaços**.
- **Conferir `-num_alignments` × `-max_target_seqs` com `-outfmt 6`:** confirmar que vêm de
  fato k hits por query, e não mais nem menos.
- Registrar no `BLASTr_metadata` os parâmetros de paralelismo e a duração.

## 8. Decisões pendentes (Heron)

1. **Janela exclusiva no servidor.** Hoje há dois usuários com BLAST no core_nt, rodando há
   mais de 2 dias. As rodadas frias despejam o cache que eles usam.
2. **Sudo** para `drop_caches`, ou instalar `vmtouch`.
3. **Banco do artigo:** congelar o `core_nt_nov26`, que vai mudar.
4. **Dados publicáveis:** as seqs de projetos de cliente podem ir para o
   repositório e o material suplementar do artigo?
5. **Quais "outros computadores"** estão disponíveis para a fase 5.
