# Hélice — exploração genética local

Aplicativo em português, executado no seu computador, para importar genótipos e variantes humanas e examinar evidências públicas. **Versão exploratória funcional, sem validação clínica.** Não promete prever todas as doenças, fazer diagnóstico ou selecionar medicamentos automaticamente.

## Iniciar no Windows

1. Abra `iniciar.cmd` nesta pasta. O navegador abre em **http://127.0.0.1:8765**.
2. Clique em **Importar DNA** e escolha seu TXT/CSV ou VCF.
3. Prefira detectar referência e fita pelo cabeçalho. Se faltarem, use as informações do laboratório original. Não adivinhe o assembly ou a orientação.
4. Consulte os painéis e as fontes. **Relatório HTML** salva um documento independente, que pode ser aberto e impresso em PDF pelo navegador. JSON e CSV estão em **Método e fontes**.
5. Para encerrar, feche o processo com **Ctrl+C** na janela do servidor. Fechar apenas a aba não encerra o processo.

O iniciador usa o Python do ambiente Codex quando disponível e, como alternativa, `py -3`. Em outro computador é necessário **Python 3.11 ou superior**. O núcleo médico e os parsers usam a biblioteca padrão. A projeção NHGRI K=21 e a reamostragem requerem **NumPy >= 1.26**, já disponível no Python deste ambiente. Em outra instalação: `python -m pip install -r requirements.txt`.

O repositório contém código, documentação, coordenadas públicas dos marcadores e exemplos sintéticos. Bases ClinVar e referências populacionais precisam ser baixadas localmente pelos utilitários. Dados pessoais, resultados e downloads não são versionados.

Alternativa no terminal:

```powershell
python server.py
```

Para abrir já com um arquivo, preservando o original:

```powershell
python server.py --input "C:\caminho\dados.txt"
```

Opções: `--port 8766`, `--no-browser`, `--build GRCh37`, `--build GRCh38`, `--strand +`, `--sample NOME`. Assembly e fita manuais devem ser confirmados pelo laboratório. O software rejeita contradição explícita com o cabeçalho.

## O que está implementado

| Área | Funcionalidade disponível | Limite principal |
|---|---|---|
| Genética médica | Correspondência local com todas as variantes simples suportadas do VCF ClinVar instalado; pesquisa de genes/condições, classificação, conflitos, estrelas e fonte. Painel inicial HFE, F5, F2, HBB. | Não transforma associação em diagnóstico, não resolve herança, penetrância ou fase. Não identifica variantes novas como patogênicas. |
| Câncer e desenvolvimento | O motor ClinVar inclui genes como BRCA1, BRCA2, genes de reparo, cardiopatias e distúrbios do desenvolvimento, quando a variante e os alelos estão efetivamente no arquivo. | Um chip não cobre integralmente esses genes. Resultado sem correspondências não é teste negativo. |
| Farmacogenética | 11 marcadores em CYP2C19, CYP2C9, SLCO1B1, ABCG2, VKORC1, DPYD e NUDT15, com contexto CPIC e fontes. | Não calcula diplótipos completos, atividade total, dose, tratamento ideal, HLA ou CNVs de CYP2D6. |
| Nutrição | Lactase, MTHFR/folato e associação GC/vitamina D; distinção entre traço, associação e exames bioquímicos. | Não diagnostica deficiência, doença celíaca ou sensibilidade não celíaca ao glúten. Não recomenda dieta ou suplementos por SNP. |
| Performance | ACTN3 R577X e seu contexto funcional. | Não prevê talento ou modalidade, nem prescreve treinamento. |
| Ancestralidade | Painel público NHGRI GRCh37 K=21 instalado, poda empírica de LD, projeção de mistura e relatório com 200 reamostragens de blocos. Aceita outros painéis. | Componentes do modelo não são países, tribos ou povos históricos. Faixas condicionadas ao painel não incluem vieses de referência ou de genotipagem. |
| Auditoria | SHA-256, referência, chamadas ausentes, duplicatas/conflitos, cromossomos, filtros e busca por rsID/coordenadas. | A taxa de chamadas mede o arquivo importado; não a cobertura do genoma nem a acurácia do laboratório. |

## Formatos escolhidos

### TXT, CSV e TSV de genotipagem

Formato 23andMe, com exemplo sintético:

```text
# reference=GRCh37
# strand=+
# rsid chromosome position genotype
rs4988235  2  136608646  AG
```

Aceita tabulação, vírgula, ponto e vírgula ou espaços no formato simples de quatro colunas. CSV com cabeçalho aceita `rsid,chromosome,position,genotype`, ou `allele1,allele2` no lugar de `genotype`. Existem aliases básicos para colunas em português. Não é um importador de qualquer planilha comercial arbitrária.

A referência e a orientação declaradas no cabeçalho são lidas pelo programa, **mas um cabeçalho inserido por conversor não comprova a preservação das informações do laboratório original**. Confira a documentação de origem antes de interpretar um arquivo real.

São interpretados alelos A/C/G/T. `--`, `00`, `NC`, `NN` e outras chamadas ausentes são excluídas. Códigos `I/D`, `II`, `DD` e indels sem sequências explícitas não são convertidos por suposição. Cromossomos 23/24/26 são X/Y/MT; código 25 (XY) não é confundido com MT.

### VCF e VCF.GZ

VCF com `GT`, uma amostra selecionada e referência GRCh37 ou GRCh38 declarada ou informada. Suporta chamadas diploides, chamadas haploides em X/Y/MT, genótipos multialélicos e indicação de fase dentro do registro. Em VCF com várias amostras é obrigatório selecionar uma pelo nome.

- FILTER diferente de `PASS` ou `.` é excluído.
- DP < 10 ou GQ < 20 é excluído, quando esses valores existem.
- DP/GQ ausentes ou FILTER `.` são contabilizados como qualidade incompleta; não são apresentados como qualidade clinicamente validada.
- Chamadas parciais, GT inválido, alelos simbólicos chamados e ploidia autossômica não diploide são excluídos.
- Omissão de uma posição no VCF não é convertida em `0/0`.
- Indels são comparados pela representação exata de REF/ALT e posição. **Não há normalização/left alignment com FASTA**; representações equivalentes diferentes podem não corresponder.

Upload pela interface: até 512 MB por arquivo. Arquivos maiores podem ser lidos com `--input`; o limite prático é a memória disponível. O programa mantém chamadas e índices em memória. Para um genoma completo, alguns GB de RAM podem ser necessários. Não é um pipeline de streaming de baixo consumo.

### Por que não FASTQ, BAM ou CRAM diretamente?

Esses arquivos exigem etapas especializadas de alinhamento, controle de qualidade e chamada de variantes antes da interpretação. Esta versão começa em **VCF**, formato adequado para integrar resultados de pipelines de exoma/genoma, ou em genótipos de chip. Não implementa basecalling, alinhamento ou variant calling.

## ClinVar: referência ampla de variantes

As bases são baixadas da fonte oficial e indexadas na instalação local. O número de registros, a edição declarada, os registros excluídos e o SHA-256 aparecem na auditoria. Nenhuma base ClinVar é incluída no Git.

Em **Bases de referência**, selecione GRCh37/38 e baixe o VCF público do NCBI. É possível importar uma cópia local oficial. O assembly do VCF é validado antes da substituição do índice.

```powershell
python tools/download_clinvar.py GRCh37
python tools/download_clinvar.py GRCh38
```

URLs fixas:

- [ClinVar GRCh37](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz)
- [ClinVar GRCh38](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz)
- [Documentação dos formatos ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/)

O download é grande e o SQLite pode ocupar mais de 1 GB. O índice guarda SHA-256, origem, data do arquivo, data de importação e assembly. Atualização é manual; a aplicação não verifica automaticamente novas classificações.

Correspondência:

1. Assembly e posição precisam concordar.
2. No VCF pessoal, REF precisa coincidir e o ALT do ClinVar precisa estar realmente chamado no GT.
3. No chip, ambos os alelos precisam ser compatíveis com a variante simples da base; rsID divergente, quando conhecido dos dois lados, impede o match.
4. Só o alelo correspondente à classificação registrada é usado; nunca se supõe que todo ALT é patogênico.
5. O app conserva classificações benignas, incertas, conflitantes e outras associações; P/LP é uma categoria da base, não um diagnóstico.

As estrelas refletem o **nível de revisão da classificação**, não a qualidade da chamada pessoal. Variantes raras em dados de chips podem ser falsos positivos. As classificações agregadas não substituem a revisão da condição específica. Não se assume que duas variantes estejam em trans, que uma variante seja causal ou que o modo de herança seja compatível.

O VCF público não representa todas as alterações do ClinVar. Existem limitações para variantes estruturais, haplótipos e alelos que coincidem com bases da referência. Nenhum resultado de ausência exclui doença. Detalhamento em memória é limitado a 20 mil associações não prioritárias; correspondências P/LP e conflitantes são preservadas e eventuais omissões são contadas no relatório.

## Preparar e utilizar o painel NHGRI

A versão 1.1 fornece um utilitário para baixar e preparar o painel oficial de projeção do NHGRI/CRGGH, com frequências A2 de 21 componentes, de Baker, Rotimi & Shriner (2017). O download original possui 19.075 SNPs em GRCh37. O painel foi combinado por posição e alelos com os genótipos públicos de 3.528 indivíduos do recurso de 2014, exclusivamente para controle empírico de LD.

O preparo exclui pares A/T e C/G quando existentes, variantes com taxa de chamada de referência <95%, MAF <1% e pares com r² >0,1 dentro de 1 Mb. A seleção gulosa em ordem genômica resultou em **15.625 marcadores públicos**. Esta poda agrupada não remove necessariamente LD de longo alcance ou LD específico de uma população. O software não imputa os genótipos pessoais.

Primeiro prepare a referência pública, sem selecionar qualquer arquivo pessoal:

```powershell
python tools/nhgri_ancestry.py --prepare-only
```

Para gerar um relatório pessoal com 200 reamostragens:

```powershell
python tools/nhgri_ancestry.py --input "C:\caminho\dados.txt" --output reports/ancestralidade_genera.json --bootstrap 200
```

O número de SNPs compatíveis e de blocos é determinado para cada arquivo. As faixas usam blocos físicos de 5 Mb e semente fixa, mantendo as frequências de referência constantes. Os percentis 2,5–97,5 indicam estabilidade condicional; não constituem validação externa e não abrangem vieses das frequências ou erros sistemáticos de genotipagem.

Os arquivos pessoais são gerados em `reports/` e não entram no Git. O JSON inclui componentes, agrupamentos, auditoria, hashes dos downloads públicos, ajuste sem cada cromossomo e controle sintético interno. A análise salva é reutilizada somente quando versão e SHA-256 do arquivo e do painel coincidirem. Outro arquivo GRCh37 é projetado na referência instalada, sem salvar relatório automaticamente.

O script baixa apenas recursos públicos. Os genótipos pessoais nunca são enviados. Um snapshot público com hash inesperado é rejeitado para revisão de versão. A referência processada está em `data/ancestry/nhgri_k21_grch37.json`; os ZIPs originais permanecem nessa pasta para auditoria.

Fontes: [NHGRI/CRGGH — painel e orientação A2](https://www.genome.gov/about-nhgri/Center-for-Research-on-Genomics-and-Global-Health/Research), [Baker et al. (2017)](https://doi.org/10.1038/s41598-017-01837-7), [formato PLINK BED](https://www.cog-genomics.org/plink/1.9/formats#bed).

## Ancestralidade com outro painel

O modelo usa frequências populacionais fixas e maximiza a verossimilhança binomial de mistura. Exige:

- Mesmo assembly e fita positiva.
- 2 a 32 componentes de referência; acima de 12, NumPy é necessário. O preparador de VCF próprio mantém o limite conservador de 12 populações.
- Pelo menos 500 SNPs autossômicos informativos em 5 cromossomos.
- Frequências válidas, alelos compatíveis, posições sem duplicatas e populações distinguíveis.
- Fonte conhecida e declaração explícita de poda de LD (`ld_pruned: true`). A aplicação não consegue auditar sozinha se essa declaração é correta.

Esquema JSON (um registro ilustrativo, insuficiente para rodar o modelo):

```json
{
  "assembly": "GRCh37",
  "strand": "+",
  "source": "Descreva a versão, populações, filtros e fonte do painel real",
  "ld_pruned": true,
  "populations": ["POP_A", "POP_B"],
  "markers": [{
    "chrom": "1", "pos": 10000,
    "other_allele": "A", "effect_allele": "G",
    "frequencies": {"POP_A": 0.1, "POP_B": 0.8}
  }]
}
```

Para preparar a referência, obtenha um VCF público consentido, por exemplo do [IGSR / 1000 Genomes](https://www.internationalgenome.org/), selecione populações, harmonize alelos e aplique filtros e poda de LD com ferramentas especializadas. Evite tratar populações já miscigenadas como ancestrais homogêneos ou usar marcadores de saúde selecionados como painel ancestral.

TSV de amostras:

```text
sample	population
SAMPLE_1	POP_A
SAMPLE_2	POP_B
```

```powershell
python tools/prepare_ancestry_panel.py --vcf referencias_filtradas_ld.vcf.gz --samples amostras.tsv --assembly GRCh37 --source "IGSR; versão e filtros usados" --confirm-ld-pruned --output painel.json
```

O preparador exige por padrão 40 alelos chamados por população e locus. O modelo rejeita cobertura insuficiente e referências praticamente idênticas. Percentuais são exploratórios, sem intervalos de confiança ou validação externa: populações omitidas, LD residual e seleção do painel podem enviesar fortemente o resultado. Não determina etnia, nacionalidade, haplogrupos ou proporções de povos históricos.

## Privacidade

- Servidor vinculado somente a **127.0.0.1**; proteção de Host/Origin e token local nas APIs.
- Sem contas, telemetria, fontes externas na interface ou serviços remotos de interpretação.
- Não envia genótipos ao NCBI, Ensembl ou outros serviços. Downloads consultam somente referências públicas.
- Upload usa temporário local, removido em sucesso ou falha; uma interrupção abrupta pode deixar esse temporário na pasta temporária do sistema.
- Genótipos são mantidos na memória da sessão e não são salvos automaticamente em banco ou relatório. A memória/paginação do sistema operacional não é apagada de forma segura pelo app.
- **Limpar dados da sessão** libera os dados da aplicação; não apaga o original nem relatórios já exportados.
- As referências públicas ficam em `data/`. O aplicativo não modifica a sincronização ou os backups das pastas do usuário.
- Relatórios HTML/JSON/CSV contêm dados genéticos sensíveis. Eles são gerados somente mediante exportação solicitada pela interface ou CLI.
- As fontes externas são abertas somente quando você clica. Sem dependências de CDNs.

## Testes e exemplos

```powershell
python -m unittest discover -s tests -v
node --check web/app.js
node tests/test_frontend.cjs
python server.py --input examples/exemplo_sintetico.txt
python server.py --input examples/exemplo_sintetico.vcf
python tools/analyze.py examples/exemplo_sintetico.vcf --output reports/exemplo.json
```

Exemplos e testes são **sintéticos**. Não representam uma pessoa. Os testes verificam assembly, fita, alelos, variantes multialélicas, chamadas ausentes, duplicatas, seleção de amostra, qualidade, correspondência REF/ALT, proteção contra classificação de referência homozigota como ALT observado e mistura ancestral sintética conhecida. Isso testa implementação; **não constitui validação clínica ou populacional**.

## Organização

```text
server.py                  HTTP local, jobs e sessão em memória
genetics.py                parsers e painel de marcadores
clinical.py                download, índice SQLite e match ClinVar
ancestry.py                modelo exploratório com painel fornecido
catalog.py                 textos e alelos de interesse do painel inicial
data/markers.json          snapshot de coordenadas públicas verificadas
web/                       interface sem dependências externas
tools/                     atualização de referências, CLI e preparo de painel
tests/                     testes com dados sintéticos
examples/                  arquivos sintéticos para demonstração
```

Atualizar coordenadas públicas dos 20 marcadores:

```powershell
python tools/fetch_markers.py
```

Isso consulta Ensembl GRCh37/38, exige mapeamento único na fita positiva e checa compatibilidade do alelo de interesse. O significado biológico do alelo é definido no catálogo e nas fontes; não é inferido automaticamente do alelo alternativo da API.

## Fontes principais

- [ClinVar — acesso, representação e limitações dos formatos](https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/)
- [CPIC — diretrizes farmacogenéticas](https://cpicpgx.org/guidelines/)
- [MedlinePlus — interpretação de testes de consumo](https://medlineplus.gov/genetics/understanding/dtcgenetictesting/dtcresults/)
- [Tandy-Connor et al. — falsos positivos e confirmação clínica](https://pubmed.ncbi.nlm.nih.gov/29565420/)
- [NIDDK — diagnóstico de doença celíaca](https://www.niddk.nih.gov/health-information/digestive-diseases/celiac-disease/diagnosis)
- [CDC — MTHFR e folato](https://www.cdc.gov/folic-acid/data-research/mthfr/index.html)
- [Consenso sobre testes genéticos e desempenho esportivo](https://pmc.ncbi.nlm.nih.gov/articles/PMC4680136/)

As fontes específicas de cada marcador e as coordenadas públicas ficam acessíveis no aplicativo e no relatório exportado.
