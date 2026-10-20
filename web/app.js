'use strict';
const $ = (id) => document.getElementById(id);
const esc = (v) => String(v ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const num = (v) => Number(v || 0).toLocaleString('pt-BR');
const link = (url, label) => `<a href="${esc(/^https:\/\//.test(url || '') ? url : '#')}" target="_blank" rel="noopener noreferrer">${esc(label)} ↗</a>`;
const meta = {
  overview:['Visão geral','Seu DNA, em perspectiva.','Explore o que seus dados podem revelar — e conheça os limites de cada análise.'],
  medical:['Saúde e doenças','Genética médica, com contexto.','Variantes observadas e evidências públicas para orientar uma investigação responsável.'],
  pharma:['Farmacogenética','A genética por trás dos medicamentos.','Marcadores de metabolismo e transporte, com fontes e limitações explícitas.'],
  nutrition:['Nutrição','Uma peça do seu metabolismo.','Tendências genéticas não substituem sintomas, exames ou avaliação nutricional.'],
  performance:['Performance','Biologia não é destino.','Explore traços musculares sem transformar um SNP em uma prescrição de treino.'],
  ancestry:['Ancestralidade','Sua história precisa de referências.','Comparações populacionais dependem de um painel compatível e bem caracterizado.'],
  explorer:['Explorar variantes','Olhe de perto seus dados.','Busque um marcador ou uma região e consulte a chamada que realmente veio no arquivo.'],
  methods:['Método e fontes','Evidência que você pode conferir.','Saiba como as análises são feitas, o que cobrem e onde cada afirmação se apoia.']
};
const labels = {observed:'Alelo observado',not_observed:'Alelo não observado',unavailable:'Não avaliado',pathogenic:'Patogênica / provável',conflicting:'Classificações conflitantes',uncertain:'Significado incerto',benign:'Benigna / provável',other:'Outra associação'};
let token = '', report = null, status = null, view = 'overview', lastJob = null, polling = false, clinicalPage = 0;
let clinicalFilter = {query:'',group:'priority',stars:'0'};
let markerFilter = 'all';
const chip = (value,text) => `<span class="chip ${esc(value)}">${esc(text || labels[value] || value)}</span>`;
const button = (text,viewName) => `<button class="link-button" data-go="${viewName}">${text} →</button>`;
const empty = (title,description,action='') => `<div class="empty"><div class="empty-symbol">⫶</div><h3>${esc(title)}</h3><p>${esc(description)}</p>${action}</div>`;
const note = (text,type='') => `<div class="note ${type}">${text}</div>`;
const progressBar = (percentage) => `<div class="progress-track"><div class="progress-fill" style="width:${Math.max(0,Math.min(100,Number(percentage)||0))}%"></div></div>`;

async function api(path,options={}) {
  const response = await fetch(path,{...options,headers:{'X-DNA-Token':token,...options.headers}});
  let data;
  try { data = await response.json(); } catch { throw new Error('O servidor não respondeu como esperado.'); }
  if(!response.ok) throw new Error(data.error || 'Não foi possível concluir a operação.');
  return data;
}
function toast(text) { $('toast').textContent=text; $('toast').hidden=false; clearTimeout(toast.timer); toast.timer=setTimeout(()=>$('toast').hidden=true,6500); }
function goto(name) {
  if(!meta[name]) name='overview';
  view=name; markerFilter='all';
  history.replaceState(null,'',`#${name}`);
  render(); window.scrollTo({top:0,behavior:'instant'});
}
function render() {
  document.querySelectorAll('.nav').forEach(b=>b.classList.toggle('active',b.dataset.view===view));
  $('breadcrumb').textContent=meta[view][0]; $('page-title').textContent=meta[view][1]; $('page-description').textContent=meta[view][2];
  $('export-html').hidden=!report;
  $('content').innerHTML=({overview,medical,pharma:()=>markerPage('pharma'),nutrition:()=>markerPage('nutrition'),performance:()=>markerPage('performance'),ancestry,explorer,methods}[view])();
  wireContent();
}

function helix() {
  let rungs='';
  for(let y=0;y<=290;y+=16) {
    const angle=y/36, x1=135+70*Math.sin(angle), x2=135-70*Math.sin(angle);
    rungs+=`<line x1="${x1}" y1="${y}" x2="${x2}" y2="${y}" stroke="#7fa68e" stroke-width="2" opacity=".6"/><circle cx="${x1}" cy="${y}" r="3.5" fill="#547e63"/><circle cx="${x2}" cy="${y}" r="3.5" fill="#a2bda5"/>`;
  }
  const path=(sign)=>Array.from({length:150},(_,i)=>`${i?'L':'M'}${135+sign*70*Math.sin(i*2/36)},${i*2}`).join(' ');
  return `<div class="helix" aria-hidden="true"><svg viewBox="0 0 270 300"><path d="${path(1)}" stroke="#5c8469" fill="none" stroke-width="2.5"/><path d="${path(-1)}" stroke="#9bb69e" fill="none" stroke-width="2.5"/>${rungs}</svg></div>`;
}

function areaCards() {
  const areas=[['medical','✚','Saúde e doenças','Investigue variantes associadas a condições hereditárias, com apoio do ClinVar.'],['pharma','◒','Farmacogenética','Conheça marcadores que podem influenciar o metabolismo de medicamentos.'],['nutrition','♧','Nutrição e performance','Explore lactase, metabolismo do folato e traços musculares.'],['ancestry','◎','Ancestralidade','Compare seus dados com populações de um painel de referência compatível.']];
  return `<div class="grid four">${areas.map(([name,icon,title,desc])=>{
    const markers=report?.panel.filter(p=>p.category===name)||[];
    const coverage=markers.filter(p=>p.status!=='unavailable').length;
    const caption=!report?'Aguardando arquivo':name==='ancestry'?(report.ancestry.available?'Referência carregada':'Referência necessária'):name==='medical'?(report.clinical.available?'ClinVar disponível':'Painel inicial disponível'):`${coverage}/${markers.length} marcadores avaliáveis`;
    return `<article class="card area-card"><div class="icon-box">${icon}</div><h3>${title}</h3><p>${desc}</p><div class="card-bottom"><span>${caption}</span><button class="link-button" data-go="${name}" aria-label="Abrir ${title}">↗</button></div></article>`;
  }).join('')}</div>`;
}

function overview() {
  if(!report) return `<section class="hero"><div class="hero-copy"><div class="eyebrow">UM NOVO OLHAR PARA SEUS DADOS</div><h2>Conhecer seu DNA começa<br>por entender o que ele diz.</h2><p>Reúna marcadores, fontes científicas e limitações em um único lugar. Tudo processado no seu computador.</p><button class="link-button text-link" data-import>Importar meu primeiro arquivo →</button></div>${helix()}<div class="hero-code">A · T · C · G</div></section>
    <div class="section-top"><h2>Quatro caminhos para explorar</h2><span class="section-label">SEUS DADOS. SEU CONTEXTO.</span></div>${areaCards()}
    <div class="bottom-grid"><section class="card"><h3>Do arquivo à evidência</h3><div class="step"><span class="step-index">01</span><div><strong>Importe o que você já tem</strong><p>TXT / CSV da Genera e 23andMe, ou VCF de variantes de sequenciamento.</p></div></div><div class="step"><span class="step-index">02</span><div><strong>Confira referência e qualidade</strong><p>Assembly, orientação, chamadas ausentes e duplicatas entram na análise.</p></div></div><div class="step"><span class="step-index">03</span><div><strong>Explore com as fontes ao lado</strong><p>Veja o que foi observado, o que não foi testado e o que precisa de confirmação.</p></div></div></section><section class="card privacy-card"><div class="privacy-art">⌂</div><h3>Seus dados continuam seus.</h3><p>Os genótipos são processados na memória do servidor local. Não usamos serviços de interpretação remota, telemetria ou contas.</p><p>Somente as referências públicas são baixadas. Relatórios são exportados quando você solicita.</p>${button('Entender o método','methods')}</section></div>`;
  const s=report.summary, available=report.panel.filter(p=>p.status!=='unavailable').length;
  const clinicalCount=report.clinical.summary.pathogenic || 0;
  return `<div class="file-strip"><div class="file-main"><div class="file-icon">≋</div><div><div class="file-name">${esc(s.name)}</div><div class="file-meta">${esc(s.format)} · ${esc(s.build)} · Fita positiva${s.sample?' · Amostra '+esc(s.sample):''}</div></div></div>${chip('not_observed','ARQUIVO PROCESSADO')}</div>
    <div class="grid four"><article class="metric"><div class="metric-label">Posições únicas importadas</div><div class="metric-value">${num(s.loci)}</div><div class="metric-sub">Não representa cobertura de todo o genoma</div></article><article class="metric"><div class="metric-label">Chamadas utilizáveis no arquivo</div><div class="metric-value accent">${num(s.call_rate)}<small>%</small></div><div class="metric-sub">${num(s.no_call)} posições sem chamada utilizável</div></article><article class="metric"><div class="metric-label">Painel inicial de marcadores</div><div class="metric-value">${available}<small> / ${report.panel.length}</small></div><div class="metric-sub">Posições avaliáveis; sem score de saúde</div></article><article class="metric"><div class="metric-label">Correspondências P / LP no ClinVar</div><div class="metric-value">${report.clinical.available?num(clinicalCount):'—'}</div><div class="metric-sub">${report.clinical.available?'Classificação da base; exige revisão clínica':'Base ClinVar ainda não instalada'}</div></article></div>
    ${note('<strong>O que este arquivo permite concluir?</strong> Uma posição sem alelo de interesse não exclui doença. Chips testam locais selecionados; mesmo VCFs não abrangem todos os tipos de alteração genética. Não existe um teste de DNA que preveja todas as doenças.')}
    <div class="section-top"><h2>Explore suas análises</h2><span class="section-label">EVIDÊNCIA ANTES DE INTERPRETAÇÃO</span></div>${areaCards()}
    <div class="bottom-grid"><section class="card"><h3>Cobertura do painel inicial</h3><p class="small muted">Contagem dos marcadores deste aplicativo. Não é cobertura clínica dos genes.</p>${[['medical','Saúde'],['pharma','Medicamentos'],['nutrition','Nutrição'],['performance','Performance']].map(([cat,title])=>{const all=report.panel.filter(p=>p.category===cat),n=all.filter(p=>p.status!=='unavailable').length;return `<div class="coverage-line"><span class="coverage-name">${title}</span>${progressBar(100*n/all.length)}<span class="coverage-value">${n} / ${all.length}</span></div>`;}).join('')}</section><section class="card"><h3>Qualidade e proveniência</h3><p class="small muted">${num(s.stats.duplicates)} registros duplicados · ${num(s.stats.conflicts)} conflitos · ${esc(s.build)}</p>${s.warnings.slice(0,2).map(w=>`<p class="small muted">${esc(w)}</p>`).join('')}${button('Ver auditoria do arquivo','explorer')}</section></div>`;
}

function markerCard(m) {
  return `<article class="result-card ${m.status}"><div class="result-meta"><span class="gene">${esc(m.gene)}</span>${chip(m.status)}</div><h3>${esc(m.title)}</h3><div class="result-genotype"><span class="genotype">${esc(m.genotype || '—')}</span><span>${esc(m.rsid)} · ${esc(m.chrom)}:${num(m.pos)}<br>${m.copies===null?'Sem interpretação':`${m.copies} cópia(s) do alelo ${esc(m.effect_allele)}`}</span></div><p>${esc(m.interpretation)}</p><div class="evidence"><span>${esc(m.evidence)}</span>${link(m.source,'Fonte')}</div></article>`;
}
function markerPage(category) {
  if(!report) return empty('Primeiro, importe seu DNA.','Os resultados aparecem somente após a leitura de um arquivo real.','<button class="button primary" data-import>Importar arquivo</button>');
  let intro='';
  if(category==='pharma') intro=note('<strong>Marcadores, não uma prescrição.</strong> Este painel não determina diplótipos completos, CNVs de CYP2D6 ou tipagem HLA. Alelos não testados, fase, interações e função renal/hepática podem mudar a interpretação. Não inicie, suspenda ou ajuste medicamentos com este relatório.');
  if(category==='nutrition') intro=note('<strong>Lactose, folato e vitamina D.</strong> Os resultados são associações específicas e não medem absorção ou concentração de vitaminas. Doença celíaca exige avaliação própria; SNPs isolados não determinam HLA-DQ2/DQ8 nem diagnosticam “sensibilidade ao glúten”.');
  if(category==='performance') intro=note('<strong>Nenhum SNP define seu potencial.</strong> ACTN3 é um traço funcional estudado, mas não prevê de forma útil o treino ideal ou o sucesso esportivo. Não são geradas recomendações de intensidade, modalidade ou suplementação a partir desse marcador.');
  const all=report.panel.filter(m=>m.category===category), visible=all.filter(m=>markerFilter==='all'||m.status===markerFilter);
  return `${intro}<div class="section-top"><div><h2>Painel exploratório</h2><p>${all.filter(m=>m.status!=='unavailable').length} de ${all.length} marcadores avaliáveis · ${esc(report.summary.build)} · alelos na fita positiva</p></div></div><div class="pill-tabs">${[['all','Todos'],['observed','Alelo observado'],['not_observed','Não observado'],['unavailable','Não avaliado']].map(([v,t])=>`<button data-marker-filter="${v}" class="${markerFilter===v?'active':''}">${t}</button>`).join('')}</div><div class="report-cards">${visible.map(markerCard).join('')}</div>${!visible.length?empty('Nenhum marcador neste filtro.','Altere o filtro para ver os demais registros.'):''}${category==='nutrition'?`<div class="section-top"><h2>Glúten e doença celíaca</h2></div><article class="card"><h3>HLA não determinado por este painel</h3><p class="small muted">A presença de tipos HLA compatíveis indica suscetibilidade, mas não confirma doença celíaca. Não inferimos HLA por SNPs substitutos. Sensibilidade não celíaca ao glúten não possui um diagnóstico genético estabelecido por este aplicativo.</p>${link('https://www.niddk.nih.gov/health-information/digestive-diseases/celiac-disease/diagnosis','Diagnóstico de doença celíaca • NIDDK')}</article>`:''}`;
}

function clinicalMatches() {
  if(!report?.clinical.available) return [];
  const query=clinicalFilter.query.toLowerCase().trim();
  return report.clinical.findings.filter(f=>(clinicalFilter.group==='all'||(clinicalFilter.group==='priority'?['pathogenic','conflicting'].includes(f.group):f.group===clinicalFilter.group)) && f.stars>=Number(clinicalFilter.stars) && (!query||`${f.gene} ${f.conditions} ${f.rsid} ${f.variation_id} ${f.classification}`.toLowerCase().includes(query)));
}
function clinicalTable() {
  const filtered=clinicalMatches(),pages=Math.max(1,Math.ceil(filtered.length/25));
  clinicalPage=Math.min(clinicalPage,pages-1);
  const shown=filtered.slice(clinicalPage*25,(clinicalPage+1)*25);
  return `<div class="table-wrap"><table><thead><tr><th>Gene / variante</th><th>Condição registrada no ClinVar</th><th>Genótipo</th><th>Classificação</th><th>Revisão</th><th></th></tr></thead><tbody>${shown.map(f=>`<tr><td><span class="table-gene">${esc(f.gene.split('|').map(g=>g.split(':')[0]).join(', '))}</span><br><span class="mono">${esc(f.rsid || 'VCV '+f.variation_id)}</span></td><td><span class="table-condition">${esc(f.conditions.slice(0,160))}${f.conditions.length>160?'…':''}</span></td><td class="mono">${esc(f.genotype)}<br><small>${f.copies} × ${esc(f.alt)}</small></td><td>${chip(f.group)}</td><td title="${esc(f.review)}">${'★'.repeat(f.stars)}${'☆'.repeat(4-f.stars)}<br><small>${f.stars}/4</small></td><td><button class="link-button" data-detail="${esc(f.variation_id)}">Detalhes ↗</button></td></tr>`).join('')||'<tr><td colspan="6">Nenhuma correspondência neste filtro. Isso não exclui doença ou variantes não testadas.</td></tr>'}</tbody></table></div><div class="pagination"><span>${num(filtered.length)} correspondências no filtro · página ${clinicalPage+1} de ${pages}</span><div><button class="button" id="previous-page" ${clinicalPage===0?'disabled':''}>←</button> <button class="button" id="next-page" ${clinicalPage>=pages-1?'disabled':''}>→</button></div></div>`;
}
function medical() {
  if(!report) return markerPage('medical');
  const c=report.clinical;
  const clinicalSection=c.available?`<div class="grid four">${[['pathogenic','Patogênica / provável'],['conflicting','Conflitante'],['uncertain','Significado incerto'],['benign','Benigna / provável']].map(([key,title])=>`<article class="metric"><div class="metric-label">${title}</div><div class="metric-value">${num(c.summary[key])}</div><div class="metric-sub">Correspondências, não doenças</div></article>`).join('')}</div><div class="filters"><label class="wide">Buscar gene, condição ou rsID<input id="clinical-query" value="${esc(clinicalFilter.query)}" placeholder="Ex.: BRCA1, BRCA2, LDLR, câncer…"></label><label>Classificação<select id="clinical-group">${[['priority','Prioritárias: P/LP e conflitos'],['all','Todas as classificações'],...Object.keys(labels).filter(k=>['pathogenic','conflicting','uncertain','benign','other'].includes(k)).map(k=>[k,labels[k]])].map(([v,t])=>`<option value="${v}" ${v===clinicalFilter.group?'selected':''}>${t}</option>`).join('')}</select></label><label>Revisão mínima<select id="clinical-stars">${[0,1,2,3,4].map(n=>`<option value="${n}" ${String(n)===clinicalFilter.stars?'selected':''}>${n}+ estrela(s)</option>`).join('')}</select></label></div><div id="clinical-table">${clinicalTable()}</div><p class="small muted">ClinVar ${esc(c.metadata.file_date)} · ${num(c.metadata.records)} variantes públicas indexadas · ${num(c.matched_loci)} posições com alelos correspondentes.${c.omitted?' '+num(c.omitted)+' associações omitidas do detalhamento por limite de memória.':''}</p>${note(esc(c.message),'teal')}`:empty('Amplie a análise com o ClinVar.','O motor está pronto para cruzar suas variantes com todos os genes da base, incluindo BRCA1 e BRCA2. Instale o VCF público correspondente à referência do seu arquivo.','<button class="button primary" data-reference>Abrir bases de referência</button>');
  return `${note('<strong>Rastreamento exploratório de variantes.</strong> Uma classificação patogênica no ClinVar não significa que você tenha a doença ou que vá desenvolvê-la. Variantes raras em chips podem ser falsos positivos; estrelas do ClinVar não validam a chamada do seu arquivo. Herança, zigosidade, fase, penetrância, qualidade e condição específica precisam de revisão. VUS e conflitos não estabelecem predisposição.')}<div class="section-top"><div><h2>Busca ampla · ClinVar</h2><p>Correspondência exata de referência, posição e alelos. Nenhum diagnóstico automatizado.</p></div></div>${clinicalSection}
    <article class="card"><h3>BRCA1 / BRCA2 e desenvolvimento</h3><p class="small muted">Mesmo sem correspondências, estes genes não podem ser considerados “negativos”. O chip da Genera cobre posições selecionadas. Sequenciamento clínico, análise de deleções/duplicações, expansões e outros testes podem ser necessários conforme a suspeita. O ClinVar cobre variantes catalogadas, não toda possibilidade de doença.</p>${button('Consultar uma variante no arquivo','explorer')}</article>
    <div class="section-top"><div><h2>Marcadores iniciais de saúde</h2><p>Leitura pontual, independente da instalação do ClinVar.</p></div></div><div class="report-cards">${report.panel.filter(m=>m.category==='medical').map(markerCard).join('')}</div>`;
}

function ancestry() {
  const a=report?.ancestry;
  const pct=(v)=>Number(v).toLocaleString('pt-BR',{minimumFractionDigits:1,maximumFractionDigits:1});
  const rows=(values)=>values.map(p=>`<tr><td>${esc(p.label||p.population)}</td><td><strong>${pct(p.percent)}%</strong></td><td>${p.lower===undefined?'Não calculada':`${pct(p.lower)}–${pct(p.upper)}%`}</td></tr>`).join('');
  const regional=a?.regional_proportions;
  const result=a?.available?`${note(esc(a.limitation))}<article class="card"><div class="eyebrow">PROJEÇÃO EXPLORATÓRIA · FREQUÊNCIAS DE REFERÊNCIA FIXAS</div><h2>${regional?'Agrupamentos dos componentes do modelo':'Componentes populacionais estimados'}</h2><p class="small muted">${num(a.markers)} SNPs utilizáveis · ${a.chromosomes} cromossomos · ${esc(a.assembly)}</p>${(regional||a.proportions).map(p=>`<div class="ancestry-row"><strong>${esc(p.label||p.population)}</strong>${progressBar(p.percent)}<span>${pct(p.percent)}%</span></div>`).join('')}${regional?`<p class="small muted">Agrupamentos são somas dos componentes abaixo, não categorias étnicas. A separação entre Europa, norte da África e Ásia ocidental depende do modelo.</p><div class="table-wrap"><table><thead><tr><th>Agrupamento</th><th>Estimativa</th><th>Faixa de reamostragem</th></tr></thead><tbody>${rows(regional)}</tbody></table></div>`:''}<p class="small muted" style="margin-top:18px">Referência: ${esc(a.source)}</p><p class="small muted">${esc(a.method)} ${a.converged?'Convergência numérica atingida.':'Convergência não atingida; os percentuais não devem ser considerados estáveis.'}</p>${a.uncertainty?`<p class="small muted">${a.uncertainty.replicates} reamostragens em ${a.uncertainty.blocks} blocos de 5 Mb. ${esc(a.uncertainty.method)} ${esc(a.uncertainty.limitation)}</p>`:''}</article>${regional?`<div class="section-top"><h2>Os 21 componentes originais</h2></div><div class="table-wrap"><table><thead><tr><th>Componente do painel</th><th>Estimativa</th><th>Faixa de reamostragem</th></tr></thead><tbody>${rows(a.proportions)}</tbody></table></div>`:''}`:empty('Percentuais ainda não estimados.',a?.message||'Importe seu DNA e prepare uma referência populacional compatível para iniciar a projeção.');
  return `${result}<div class="bottom-grid"><article class="card"><h3>${status?.ancestry_reference?.available?'Painel NHGRI instalado':'Preparar referência NHGRI'}</h3><p class="small muted">Prepare a referência pública GRCh37 com 21 componentes usando tools/nhgri_ancestry.py --prepare-only. Quando instalada, ela é aplicada ao importar um arquivo compatível. Uma análise salva só é reutilizada quando versão e hashes do arquivo e do painel coincidem.</p><p class="small muted">Para outros arquivos, a interface calcula a projeção; as 200 reamostragens e o relatório detalhado são gerados pelo utilitário tools/nhgri_ancestry.py. NumPy é necessário.</p><label class="button">↑ Importar outro painel JSON<input type="file" id="ancestry-file" accept=".json" hidden ${report?'':'disabled'}></label><p class="small muted" style="margin-top:14px">Painéis personalizados: 2–32 componentes, pelo menos 500 SNPs informativos em 5 cromossomos, fonte e controle de LD.</p></article><article class="card privacy-card"><h3>Componentes não são identidades.</h3><p>O resultado depende das populações incluídas no modelo. Não determina nacionalidade, pertencimento étnico, tribos nem proporções de povos históricos.</p><p>Componentes pequenos podem representar ruído, compartilhamento antigo ou limitações da referência. Faixas de reamostragem não corrigem vieses do painel ou erros sistemáticos do arquivo.</p>${link('https://www.genome.gov/about-nhgri/Center-for-Research-on-Genomics-and-Global-Health/Research','Painel público NHGRI / CRGGH')}</article></div>`;
}

function explorer() {
  if(!report) return markerPage('medical');
  const s=report.summary,largest=Math.max(...Object.values(s.chromosomes));
  return `<article class="card"><h3>Buscar nos seus genótipos</h3><p class="small muted">Informe um rsID ou uma coordenada na referência ${esc(s.build)}. Até 100 resultados por busca.</p><form id="search-form" class="search-row"><input id="variant-query" aria-label="rsID ou coordenada" placeholder="rs4988235 ou 2:136608646 ou 17:43044000-43045000" required><button class="button primary">Buscar variante →</button></form><div id="search-results" aria-live="polite"></div></article>
    <div class="section-top"><h2>Auditoria da importação</h2></div><article class="card"><dl class="key-value"><dt>Formato</dt><dd>${esc(s.format)}</dd><dt>Referência</dt><dd>${esc(s.build)}</dd><dt>Fita</dt><dd>Genômica positiva</dd><dt>Posições utilizáveis</dt><dd>${num(s.called)} / ${num(s.loci)} (${num(s.call_rate)}%)</dd><dt>Heterozigosidade autossômica</dt><dd>${num(s.heterozygosity)}% entre SNPs diploides utilizáveis; não é um diagnóstico de qualidade.</dd><dt>Duplicatas</dt><dd>${num(s.stats.duplicates)} registros; ${num(s.stats.conflicts)} conflitos excluídos</dd><dt>Qualidade incompleta no VCF</dt><dd>${num(s.stats.quality_incomplete)} registros</dd><dt>Contigs não suportados</dt><dd>${num(s.stats.unsupported_contigs)} registros ignorados</dd><dt>SHA-256 do arquivo</dt><dd class="mono">${esc(s.sha256)}</dd>${Object.entries(s.excluded).map(([k,v])=>`<dt>${esc(k)}</dt><dd>${num(v)} posições</dd>`).join('')}</dl></article>
    <div class="section-top"><h2>Distribuição de posições por cromossomo</h2><span class="section-label">NÃO É COBERTURA DO GENOMA</span></div><div class="chr-grid">${Object.entries(s.chromosomes).map(([ch,n])=>`<div class="chr-bar"><strong>CHR ${esc(ch)}</strong>${progressBar(100*n/largest)}<span>${num(n)} posições</span></div>`).join('')}</div>
    ${s.warnings.map(w=>note(esc(w),'teal')).join('')}`;
}

function methods() {
  const sources=[
    ['https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/','ClinVar · NCBI','Classificações de variantes, condições e status de revisão. VCF GRCh37 / GRCh38 obtido e indexado localmente.'],
    ['https://rest.ensembl.org/','Ensembl REST / dbSNP','Coordenadas e alelos do painel inicial, verificados em GRCh37 e GRCh38. Snapshot local em data/markers.json.'],
    ['https://cpicpgx.org/guidelines/','CPIC / ClinPGx','Diretrizes farmacogenéticas. Este aplicativo mostra marcadores associados, sem implementar prescrições ou diplótipos completos.'],
    ['https://medlineplus.gov/genetics/understanding/dtcgenetictesting/dtcresults/','MedlinePlus Genetics','Alcance e limites da interpretação de testes genéticos de consumo.'],
    ['https://pubmed.ncbi.nlm.nih.gov/29565420/','Tandy-Connor et al. · confirmação de dados de consumo','Estudo sobre falsos positivos em dados brutos e necessidade de confirmação clínica. Não se aplica uma taxa fixa de erro ao seu arquivo.'],
    ['https://www.niddk.nih.gov/health-information/digestive-diseases/celiac-disease/diagnosis','NIDDK · doença celíaca','HLA, avaliação diagnóstica e limitações da interpretação genética isolada.'],
    ['https://www.cdc.gov/folic-acid/data-research/mthfr/index.html','CDC · MTHFR e folato','Variantes comuns não implicam necessidade automática de formas específicas de folato.'],
    ['https://pmc.ncbi.nlm.nih.gov/articles/PMC4680136/','Consenso sobre genética e esporte','Limites da previsão de desempenho e seleção de treinamento a partir de testes comerciais.'],
    ['https://www.genome.gov/about-nhgri/Center-for-Research-on-Genomics-and-Global-Health/Research','NHGRI / CRGGH · painel ancestral','Painel de projeção com 21 componentes, frequências do alelo A2, coordenadas GRCh37 e genótipos públicos para controle empírico de LD.'],
    ['https://doi.org/10.1038/s41598-017-01837-7','Baker, Rotimi & Shriner (2017)','Estudo e frequências de componentes do painel global. A projeção implementada aqui não substitui a validação externa.']
  ];
  return `<div class="grid two"><article class="card"><h3>O que o motor faz</h3><ol class="notice-list"><li>Valida assembly, orientação e estrutura do arquivo.</li><li>Exclui chamadas ausentes, conflitantes ou com qualidade insuficiente quando informada.</li><li>Compara coordenadas e alelos com o painel inicial e com o VCF público do ClinVar.</li><li>Expõe classificação, condição, revisão, fonte e limitações para cada achado.</li></ol></article><article class="card"><h3>Formatos escolhidos</h3><p class="small muted"><strong>TXT / CSV / TSV:</strong> rsid, chromosome, position, genotype; admite allele1 + allele2 com cabeçalho. Alelos A/C/G/T na fita positiva.</p><p class="small muted"><strong>VCF / VCF.GZ:</strong> variantes já chamadas, GT haploide ou diploide, amostra selecionada, GRCh37/38. DP ≥ 10 e GQ ≥ 20 quando presentes.</p><p class="small muted"><strong>FASTQ / BAM / CRAM:</strong> precisam de alinhamento, controle de qualidade e chamada de variantes em pipeline externo antes de gerar um VCF.</p></article></div>
    <div class="section-top"><h2>Limites que fazem parte do resultado</h2></div><article class="card"><ul class="notice-list">${(report?.limitations||['Ferramenta exploratória, sem validação clínica. Não prevê todas as doenças nem diagnostica condições.','Ausência em chip ou VCF não é resultado negativo. Nenhuma imputação é aplicada.','Não avalia CNVs, expansões, rearranjos, HLA, fase entre registros ou diplótipos completos.','Indels no VCF precisam corresponder exatamente à representação ClinVar; não fazemos normalização com FASTA.','Não calcula risco absoluto, doses ou escores poligênicos clínicos.']).map(l=>`<li>${esc(l)}</li>`).join('')}</ul><p class="small muted">As estrelas do ClinVar representam nível de revisão da classificação, não a qualidade do seu genótipo. Classificações podem mudar; a data da base aparece no relatório.</p></article>
    <div class="section-top"><h2>Fontes e documentação</h2><span class="section-label">CONSULTADAS EM 19 OUT 2026</span></div><div class="source-list">${sources.map(([url,title,description])=>`<article class="source-item">${link(url,title)}<p>${esc(description)}</p></article>`).join('')}</div>
    <div class="section-top"><h2>Privacidade e exportação</h2></div><article class="card"><p class="small muted">O servidor aceita conexões apenas em 127.0.0.1. Genótipos ficam na memória; uploads usam um arquivo temporário local que é apagado após a leitura. Uma interrupção abrupta pode deixar esse temporário no diretório temporário do sistema. O arquivo original permanece onde você o guardou.</p><p class="small muted">Limpar a sessão ou encerrar o processo libera os dados da aplicação. Isso não apaga relatórios já exportados, backups, arquivo de paginação ou o arquivo original. Bases públicas são mantidas na pasta data. Abrir uma fonte externa acessa aquele site normalmente.</p><div class="heading-actions"><button class="button" id="export-json" ${report?'':'disabled'}>↓ Relatório JSON</button><button class="button" id="export-csv" ${report?'':'disabled'}>↓ Variantes ClinVar CSV</button></div></article>`;
}

function showDetail(id) {
  const f=report.clinical.findings.find(f=>f.variation_id===id); if(!f) return;
  $('detail-content').innerHTML=`<div class="eyebrow">CORRESPONDÊNCIA EXATA · CLINVAR</div><h2>${esc(f.gene.split('|').map(g=>g.split(':')[0]).join(', ') || 'Variante registrada')}</h2>${chip(f.group)}<dl class="key-value"><dt>Identificador</dt><dd>${esc(f.rsid || '—')} · VCV ${esc(f.variation_id)}</dd><dt>Posição</dt><dd>${esc(report.summary.build)} ${esc(f.chrom)}:${num(f.pos)} ${esc(f.ref)} → ${esc(f.alt)}</dd><dt>Chamada</dt><dd>${esc(f.genotype)} · ${f.copies} cópia(s) do ALT</dd><dt>Classificação original</dt><dd>${esc(f.classification)}</dd><dt>Revisão</dt><dd>${f.stars}/4 · ${esc(f.review.replaceAll('_',' '))}</dd><dt>Condições na base</dt><dd>${esc(f.conditions)}</dd><dt>Conflitos</dt><dd>${esc(f.conflict || 'Não informados nesse campo')}</dd><dt>Qualidade</dt><dd>${esc(f.quality)}</dd></dl>${note('A associação registrada na base não determina se a condição se aplica a você. Modo de herança, fase, penetrância e confirmação laboratorial não são resolvidos por este cruzamento.')}${link(f.source,'Abrir registro original no ClinVar')}`;
  $('detail-dialog').showModal();
}

function wireTable() {
  $('previous-page')?.addEventListener('click',()=>{clinicalPage--;updateClinicalTable();});
  $('next-page')?.addEventListener('click',()=>{clinicalPage++;updateClinicalTable();});
  document.querySelectorAll('[data-detail]').forEach(b=>b.addEventListener('click',()=>showDetail(b.dataset.detail)));
}
function updateClinicalTable() { $('clinical-table').innerHTML=clinicalTable(); wireTable(); }
function wireContent() {
  document.querySelectorAll('[data-go]').forEach(b=>b.addEventListener('click',()=>goto(b.dataset.go)));
  document.querySelectorAll('[data-import]').forEach(b=>b.addEventListener('click',()=>$('import-dialog').showModal()));
  document.querySelectorAll('[data-reference]').forEach(b=>b.addEventListener('click',openReferences));
  document.querySelectorAll('[data-marker-filter]').forEach(b=>b.addEventListener('click',()=>{markerFilter=b.dataset.markerFilter;render();}));
  $('clinical-query')?.addEventListener('input',e=>{clinicalFilter.query=e.target.value;clinicalPage=0;updateClinicalTable();});
  $('clinical-group')?.addEventListener('change',e=>{clinicalFilter.group=e.target.value;clinicalPage=0;updateClinicalTable();});
  $('clinical-stars')?.addEventListener('change',e=>{clinicalFilter.stars=e.target.value;clinicalPage=0;updateClinicalTable();});
  wireTable();
  $('search-form')?.addEventListener('submit',async e=>{
    e.preventDefault(); const target=$('search-results'); target.textContent='Buscando…';
    try {
      const {results}=await api('/api/search?q='+encodeURIComponent($('variant-query').value.trim()));
      target.innerHTML=results.length?`<div class="table-wrap"><table><thead><tr><th>Marcador</th><th>Posição</th><th>Genótipo</th><th>REF / ALT</th><th>Qualidade / exclusão</th></tr></thead><tbody>${results.map(r=>`<tr><td>${esc(r.rsids.join(', ')||'—')}</td><td>${esc(r.chrom)}:${num(r.pos)}</td><td class="mono">${esc(r.genotype)}</td><td class="mono">${esc(r.ref||'—')} / ${esc(r.alts.join(', ')||'—')}</td><td>${esc(r.reason||'Chamada utilizável; não confirma validade clínica')}<br><small>DP ${esc(r.dp??'—')} · GQ ${esc(r.gq??'—')}</small></td></tr>`).join('')}</tbody></table></div>`:note('Posição ou rsID não encontrado. Ausência de registro não significa ausência de uma variante.');
    } catch(err) {target.innerHTML=`<p class="error">${esc(err.message)}</p>`;}
  });
  $('ancestry-file')?.addEventListener('change',async e=>{
    if(!e.target.files[0]) return;
    try {await upload(e.target.files[0],{kind:'ancestry'});toast('Painel ancestral recebido; calculando o modelo.');await poll();} catch(err){toast(err.message);}
  });
  $('export-json')?.addEventListener('click',()=>saveBlob(JSON.stringify(report,null,2),'application/json','helice-relatorio.json'));
  $('export-csv')?.addEventListener('click',exportCSV);
}

function saveBlob(content,type,filename) {
  const blob=new Blob([content],{type}),url=URL.createObjectURL(blob),a=document.createElement('a');
  a.href=url; a.download=filename; a.click(); setTimeout(()=>URL.revokeObjectURL(url),3000);
}
function exportCSV() {
  const keys=['gene','rsid','chrom','pos','ref','alt','genotype','copies','classification','conditions','review','stars','source'];
  const cell=v=>'"'+String(v??'').replace(/^(\s*[=+@-])/,"'$1").replaceAll('"','""')+'"';
  const rows=(report?.clinical.findings||[]).map(f=>keys.map(k=>cell(f[k])).join(';'));
  saveBlob('\ufeff'+[keys.join(';'),...rows].join('\r\n'),'text/csv;charset=utf-8','helice-clinvar.csv');
}
function exportHTML() {
  if(!report) return;
  const s=report.summary,c=report.clinical,a=report.ancestry;
  const ancestryHTML=ancestryExportHTML(a);
  const findings=c.findings||[];
  const html=`<!doctype html><html lang="pt-BR"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><meta http-equiv="Content-Security-Policy" content="default-src 'none'; style-src 'unsafe-inline'; base-uri 'none'"><title>Hélice • Relatório exploratório</title><style>body{font:14px/1.7 system-ui,sans-serif;color:#18343e;margin:40px auto;max-width:1100px;padding:0 22px}h1,h2,h3{line-height:1.3}h1{font-size:34px}h2{margin-top:36px}table{border-collapse:collapse;width:100%;font-size:11px}th,td{border:1px solid #dce3df;padding:9px;text-align:left;vertical-align:top}th{background:#eef3eb}.card{border:1px solid #dce3df;padding:16px;margin:12px 0;break-inside:avoid}.note{background:#fbf4e5;padding:18px}a{color:#18766b}small{color:#687c80}code{overflow-wrap:anywhere}@media print{body{margin:0;max-width:none}h2{break-after:avoid}tr{break-inside:avoid}}</style><h1>hélice. / Relatório exploratório</h1><p>${esc(s.name)} · ${esc(s.build)} · ${esc(s.format)}</p><p>${num(s.loci)} posições únicas · ${num(s.called)} chamadas utilizáveis · ${num(s.call_rate)}% do arquivo</p><p><small>SHA-256: <code>${esc(s.sha256)}</code></small></p><div class="note">${report.limitations.map(l=>`<p>${esc(l)}</p>`).join('')}</div><h2>Auditoria</h2>${s.warnings.map(w=>`<p>${esc(w)}</p>`).join('')}<p>Duplicatas: ${num(s.stats.duplicates)}; conflitos: ${num(s.stats.conflicts)}.</p><h2>Marcadores iniciais</h2>${report.panel.map(m=>`<article class="card"><h3>${esc(m.title)} · ${esc(m.gene)}</h3><p>${esc(m.rsid)} · ${esc(m.chrom)}:${num(m.pos)} · ${esc(m.genotype||'Não avaliado')} · ${esc(labels[m.status])} · alelo de interesse ${esc(m.effect_allele)}</p><p>${esc(m.interpretation)}</p><p><small>${esc(m.evidence)} · ${link(m.source,'Fonte')} · ${link(m.coordinate_source,'Coordenadas')}</small></p></article>`).join('')}<h2>ClinVar</h2>${c.available?`<p>Edição ${esc(c.metadata.file_date)} · ${num(c.metadata.records)} variantes na base · ${num(c.total_findings)} correspondências. ${esc(c.message)}</p><p>Fonte: ${link(c.metadata.source,'ClinVar')} · SHA-256 da referência: <code>${esc(c.metadata.sha256)}</code></p><table><thead><tr><th>Gene / variante</th><th>Chamada</th><th>Classificação / revisão</th><th>Condições registradas</th></tr></thead><tbody>${findings.map(f=>`<tr><td>${esc(f.gene)}<br>${esc(f.chrom)}:${num(f.pos)} ${esc(f.ref)}→${esc(f.alt)}<br>${link(f.source,f.rsid||f.variation_id)}</td><td>${esc(f.genotype)}<br>${f.copies} × ALT</td><td>${esc(f.classification)}<br>${f.stars}/4 · ${esc(f.review)}</td><td>${esc(f.conditions)}</td></tr>`).join('')}</tbody></table>${c.omitted?`<p>${num(c.omitted)} associações omitidas do detalhamento.</p>`:''}`:'<p>Base não instalada. Busca ampla não realizada.</p>'}<h2>Ancestralidade</h2>${ancestryHTML}<hr><p><small>Hélice ${esc(report.version)} · Documento privado exportado localmente. Este relatório contém informações genéticas sensíveis.</small></p></html>`;
  saveBlob(html,'text/html;charset=utf-8','helice-relatorio.html');toast('Relatório HTML exportado. Abra o arquivo para ler ou imprimir em PDF.');
}

function ancestryExportHTML(a) {
  if(!a.available)return '<p>Ancestralidade não estimada: '+esc(a.message||'painel não disponível')+'</p>';
  const n=v=>Number(v).toLocaleString('pt-BR',{minimumFractionDigits:1,maximumFractionDigits:1});
  const table=values=>`<table><thead><tr><th>Componente / agrupamento</th><th>Estimativa</th><th>Faixa de reamostragem</th></tr></thead><tbody>${values.map(p=>`<tr><td>${esc(p.label||p.population)}</td><td>${n(p.percent)}%</td><td>${p.lower===undefined?'Não calculada':`${n(p.lower)}–${n(p.upper)}%`}</td></tr>`).join('')}</tbody></table>`;
  return `<p>${esc(a.method)} · ${num(a.markers)} SNPs · ${a.chromosomes} cromossomos · Fonte: ${esc(a.source)}</p>${a.regional_proportions?'<h3>Agrupamentos do modelo</h3>'+table(a.regional_proportions):''}<h3>Componentes de referência</h3>${table(a.proportions)}${a.uncertainty?`<p>${a.uncertainty.replicates} reamostragens em ${a.uncertainty.blocks} blocos de 5 Mb. ${esc(a.uncertainty.method)} ${esc(a.uncertainty.limitation)}</p>`:''}<p>${esc(a.limitation)}</p>`;
}

async function upload(file,params) {
  const qs=new URLSearchParams({...params,name:file.name});
  const result=await api('/api/upload?'+qs,{method:'POST',headers:{'Content-Type':'application/octet-stream'},body:file});
  lastJob=null;
  return result;
}
function renderDatabases() {
  $('database-list').innerHTML=(status?.databases||[]).map(db=>`<div class="database-row"><strong>${esc(db.build)}</strong> ${db.available?chip('not_observed','Instalada'):chip('unavailable','Não instalada')}<small>${db.available?`${num(db.records)} variantes · edição ${esc(db.file_date)} · indexada ${esc(db.imported)}`:'Baixe a base correspondente ao seu arquivo.'}</small></div>`).join('');
}
function openReferences() {
  if(report) $('reference-build').value=report.summary.build;
  renderDatabases(); $('reference-error').textContent=''; $('reference-dialog').showModal();
}
async function poll() {
  if(polling) return;
  polling=true;
  try {
    status=await api('/api/status');
    const j=status.job,key=JSON.stringify(j);
    $('job-banner').hidden=!j.running&&!j.error;
    $('job-banner').className='job-banner'+(j.running?' busy':'')+(j.error?' failed':'');
    $('job-banner').textContent=j.error||`${j.message}${j.count?' · '+num(j.count):''}`;
    $('import-button').disabled=j.running; $('download-reference').disabled=j.running;
    if(!j.running&&status.loaded&&(!report||lastJob!==key)) {
      report=await api('/api/report');render();renderDatabases();
    }
    if(!status.loaded&&report&&!j.running) {report=null;render();}
    lastJob=key;
  } catch(err) { $('job-banner').hidden=false;$('job-banner').className='job-banner failed';$('job-banner').textContent='Conexão com o servidor local interrompida. Abra novamente o aplicativo para continuar.'; }
  finally {polling=false;}
}

document.querySelectorAll('[data-view]').forEach(b=>b.addEventListener('click',()=>goto(b.dataset.view)));
$('import-button').addEventListener('click',()=>$('import-dialog').showModal());
$('reference-button').addEventListener('click',openReferences);
$('export-html').addEventListener('click',exportHTML);
$('dna-file').addEventListener('change',()=>{$('selected-filename').textContent=$('dna-file').files[0]?.name||'Escolha ou arraste seu arquivo';});
$('dropzone').addEventListener('dragover',e=>{e.preventDefault();$('dropzone').classList.add('dragging');});
$('dropzone').addEventListener('dragleave',()=>$('dropzone').classList.remove('dragging'));
$('dropzone').addEventListener('drop',e=>{e.preventDefault();$('dropzone').classList.remove('dragging');if(e.dataTransfer.files.length){$('dna-file').files=e.dataTransfer.files;$('dna-file').dispatchEvent(new Event('change'));}});
$('upload-form').addEventListener('submit',async e=>{
  e.preventDefault(); const file=$('dna-file').files[0];if(!file)return;
  $('upload-error').textContent='';$('upload-submit').disabled=true;$('upload-submit').textContent='Transferindo para o servidor local…';
  try {await upload(file,{kind:'dna',build:$('build-select').value,strand:$('strand-select').value,sample:$('sample-input').value.trim()});$('import-dialog').close();lastJob=null;await poll();}
  catch(err){$('upload-error').textContent=err.message;}
  finally{$('upload-submit').disabled=false;$('upload-submit').textContent='Importar e analisar →';}
});
$('download-reference').addEventListener('click',async()=>{
  try {await api('/api/download-clinvar?build='+$('reference-build').value,{method:'POST'});$('reference-dialog').close();lastJob=null;await poll();}
  catch(err){$('reference-error').textContent=err.message;}
});
$('reference-file').addEventListener('change',async e=>{
  if(!e.target.files[0])return;
  try{await upload(e.target.files[0],{kind:'clinvar',build:$('reference-build').value});$('reference-dialog').close();lastJob=null;await poll();}
  catch(err){$('reference-error').textContent=err.message;}
});
$('clear-button').addEventListener('click',async()=>{
  try{await api('/api/clear',{method:'POST'});report=null;clinicalFilter={query:'',group:'priority',stars:'0'};goto('overview');toast('Genótipos removidos da memória da sessão. Arquivos originais e relatórios exportados foram preservados.');}
  catch(err){toast(err.message);}
});
window.addEventListener('hashchange',()=>{view=meta[location.hash.slice(1)]?location.hash.slice(1):'overview';render();});
(async()=>{
  view=meta[location.hash.slice(1)]?location.hash.slice(1):'overview';render();
  try {const session=await api('/api/session');token=session.token;await poll();setInterval(poll,1500);} catch(err){toast('Não foi possível conectar ao servidor local.');}
})();
