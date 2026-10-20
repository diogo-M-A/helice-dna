# Referências externas

Os índices ClinVar, os ZIPs públicos de genótipos e o painel populacional NHGRI processado não fazem parte do repositório. Eles são obtidos pelos scripts a partir dos provedores originais e mantidos apenas na instalação local. Consulte os termos e as condições dos respectivos provedores antes de reutilizar ou redistribuir esses dados.

| Recurso | Uso no programa | Referência |
|---|---|---|
| ClinVar / NCBI | Classificações agregadas de variantes, condições e revisão | https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/ |
| Ensembl REST / dbSNP | Coordenadas e alelos dos 20 marcadores iniciais | https://rest.ensembl.org/ |
| NHGRI / CRGGH | Frequências A2 de 21 componentes e genótipos públicos para poda de LD | https://www.genome.gov/about-nhgri/Center-for-Research-on-Genomics-and-Global-Health/Research |
| Baker, Rotimi & Shriner (2017) | Origem científica do painel de projeção | https://doi.org/10.1038/s41598-017-01837-7 |
| PLINK | Especificação BED usada pelo leitor implementado neste projeto; nenhum binário PLINK é distribuído | https://www.cog-genomics.org/plink/1.9/formats#bed |
| CPIC / ClinPGx | Links e contexto de marcadores farmacogenéticos; não são implementadas recomendações completas | https://cpicpgx.org/guidelines/ |

O snapshot `data/markers.json` contém coordenadas, identificadores e alelos de referência públicos, com URLs e data de obtenção. Não contém genótipos pessoais.

Os exemplos de entrada e as fixtures de testes são sintéticos. A captura da interface foi produzida com a sessão vazia.

Não foi atribuída uma licença de distribuição ao código neste repositório privado. Torná-lo público com licença explícita e revisar a redistribuição de dados de terceiros são decisões posteriores do titular do projeto.
