# Estado da validação

Este documento separa verificação de implementação de evidência de acurácia.

## O que foi verificado

- Parsers TXT/CSV/VCF: referência, fita, chamadas ausentes/parciais, duplicatas, ploidia, genótipos multialélicos e seleção de amostra.
- Matching ClinVar com fixtures sintéticas: posição, REF/ALT, referência homozigota, incompatibilidade de alelos, classificação conflitante e atualização de base sem perda da anterior.
- Projeção de mistura: soluções conhecidas em 21 componentes, inversão de alelos, fronteira do simplex, equivalência de pesos e reamostragem reproduzível.
- Leitor PLINK BED: códigos A2, ausências, assinatura e dimensões.
- Poda empírica de LD: remoção de pares correlacionados e exclusão de SNPs palindrômicos.
- Cache pessoal: exige SHA-256 do arquivo de entrada, SHA-256 do painel e versão compatíveis.
- HTTP local: Host/Origin, token, upload sintético, falha de importação preservando a sessão e limpeza.
- Exportação HTML/CSV: escape de conteúdo, URLs, restrições de conteúdo, fórmulas CSV e inclusão das faixas de reamostragem.

Na preparação inicial do repositório, passaram 50 testes Python e os testes JavaScript de exportação. A integração contínua executa a suíte sem baixar bases ou acessar dados pessoais. A contagem não é uma métrica de acurácia biológica.

## O que ainda não foi demonstrado

1. **Concordância em amostras externas independentes.** Não foi comparado um conjunto cego de indivíduos com um método estabelecido e referência previamente definida.
2. **Acurácia em chips comerciais.** Não existe estudo de sensibilidade, especificidade ou falsos positivos deste software por variante, gene e plataforma.
3. **Calibração dos componentes de ancestralidade.** Recuperar proporções em simulação gerada pelo próprio modelo não valida o modelo para pessoas reais.
4. **Cobertura real das faixas de reamostragem.** Os percentis condicionam nas frequências fixas; não incorporam vieses de referência, escolhas de K, incerteza de frequência ou erros sistemáticos do conversor.
5. **Validação clínica ou farmacogenética.** Não implementa diagnóstico, penetrância, recomendações terapêuticas completas ou diplótipos abrangentes.
6. **Novidade científica.** A integração de parsers, ClinVar e uma projeção de frequências é uma contribuição de software; não demonstra método estatístico novo ou descoberta biológica.

## Limitações técnicas com impacto científico

- Indels precisam ter representação idêntica; não há normalização com FASTA.
- VCF ClinVar não cobre todas as classes de variantes e alelos de interesse.
- A camada de anotação não resolve herança, fase, segregação, ancestralidade local ou causalidade.
- A poda de LD usa genótipos públicos agrupados de 2014; não elimina necessariamente LD de longo alcance ou específico de população.
- O painel K=21 e o subconjunto de marcadores podem produzir componentes correlacionados ou difíceis de distinguir.
- O app depende da referência e da orientação declaradas no arquivo de entrada; um cabeçalho convertido não prova ambas.

## Critério para uma alegação científica futura

Definir previamente a pergunta, o método comparador, as amostras independentes, as métricas e as análises de sensibilidade. Reportar erros e casos de falha, além dos resultados favoráveis. Até esse estudo existir, utilizar o termo **projeção exploratória**.
