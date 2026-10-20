"""Baixa o painel público NHGRI, verifica LD e analisa DNA apenas localmente.

Fontes: Baker, Rotimi & Shriner 2017; genótipos públicos Shriner et al. 2014.
Não transmite marcadores, genótipos, nomes ou resultados pessoais.
"""
import os
# Matrizes pequenas: evitar sobrecusto de dezenas de threads BLAS.
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse
import csv
import hashlib
import html
import io
import json
import pathlib
import sys
import time
import urllib.request
import zipfile
from datetime import datetime,timezone
from collections import Counter
import numpy as np
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]))
from genetics import ROOT,parse_file
from ancestry_projection import fit_mixture,bootstrap_blocks

DATA=ROOT/'data'/'ancestry'
PAGE='https://www.genome.gov/about-nhgri/Center-for-Research-on-Genomics-and-Global-Health/Research'
PAPER='https://doi.org/10.1038/s41598-017-01837-7'
FREQUENCIES_URL='https://www.genome.gov/sites/default/files/media/files/2019-03/reference_21_P.txt.zip'
GENOTYPES_URL='https://www.genome.gov/sites/default/files/media/files/2019-03/Ancestry_2014.zip'
FREQUENCIES_SHA='7aef883b9e08e77049998240167e5f053ae4273c403ffa0a4e2a11e91147c54e'
GENOTYPES_SHA='7b7ff043ee365c9d986347c163243e07d52737c12c107e27c6ebbe682c49c78e'
PANEL_PATH=DATA/'nhgri_k21_grch37.json'
LABELS={
    'Northern Asian':'Norte da Ásia','West-Central African':'África centro-ocidental',
    'Omotic':'Omótico (leste africano)','Northern European':'Norte da Europa',
    'Japanese':'Japonês','Eastern African':'Leste da África','Kalash':'Kalash',
    'Northern African':'Norte da África','Sino-Tibetan':'Sino-tibetano','Circumpolar':'Circumpolar',
    'Amerindian':'Indígena das Américas','Southern Asian':'Sul da Ásia','Western Asian':'Ásia ocidental',
    'Central African':'África central','Southeastern Asian':'Sudeste da Ásia','South Indian':'Sul da Índia',
    'Southern European':'Sul da Europa','Southern African':'Sul da África','Oceanian':'Oceania',
    'Western African':'Oeste da África','Arabian':'Arábico',
}
GROUPS={
    'Europa':['Northern European','Southern European'],
    'África subsaariana':['West-Central African','Omotic','Eastern African','Central African','Southern African','Western African'],
    'Indígena das Américas':['Amerindian'],
    'Norte da África / Ásia ocidental':['Northern African','Western Asian','Arabian'],
    'Sul da Ásia':['Southern Asian','South Indian','Kalash'],
    'Leste / norte da Ásia e circumpolar':['Northern Asian','Japanese','Sino-Tibetan','Circumpolar','Southeastern Asian'],
    'Oceania':['Oceanian'],
}


def say(message,count=0):
    print(f'{message}'+(f' • {count}' if count else ''),flush=True)


def fetch_public(url,expected):
    DATA.mkdir(exist_ok=True)
    target=DATA/url.rsplit('/',1)[1]
    if not target.exists():
        partial=target.with_suffix('.part')
        with urllib.request.urlopen(url,timeout=120) as response,partial.open('wb') as stream:
            while chunk:=response.read(1048576): stream.write(chunk)
        os.replace(partial,target)
    actual=hashlib.sha256(target.read_bytes()).hexdigest()
    if actual!=expected:
        raise ValueError(f'Snapshot público mudou ou está corrompido: {target.name}. Revise a versão antes de atualizar o SHA esperado.')
    return target


def decode_bed(raw,n_samples,n_markers):
    """Dosagem do alelo A2 no formato PLINK SNP-major oficial."""
    if raw[:3]!=bytes.fromhex('6c1b01'): raise ValueError('BED não está em formato PLINK SNP-major.')
    width=(n_samples+3)//4
    if len(raw)!=3+width*n_markers: raise ValueError('BED/BIM/FAM têm dimensões incompatíveis.')
    packed=np.frombuffer(raw,dtype=np.uint8,offset=3).reshape(n_markers,width)
    shifts=np.arange(4,dtype=np.uint8)*2
    codes=((packed[:,:,None]>>shifts)&3).reshape(n_markers,width*4)[:,:n_samples]
    return np.asarray([0,np.nan,1,2],dtype=np.float32)[codes]


def ld_prune(rows,reference,reference_index,threshold=.1,window=1000000,exclude_palindromic=True):
    """Poda empírica gulosa com r² de dosagens nos indivíduos públicos.

Ausências do painel público são centradas na média apenas para calcular LD;
nenhum genótipo do usuário é imputado. Não equivale a independência perfeita.
"""
    retained=[]
    active=[]
    counts=Counter()
    max_retained_r2=0.0
    for r in sorted(rows,key=lambda r:(int(r['Chromosome']),int(r['Position']))):
        ch,pos=r['Chromosome'],int(r['Position'])
        pair={r['A1'],r['A2']}
        if len(pair)!=2 or not pair<=set('ACGT'):
            counts['non_snv']+=1;continue
        if exclude_palindromic and pair in ({'A','T'},{'C','G'}):
            counts['palindromic']+=1;continue
        index=reference_index.get((ch,pos))
        if index is None:
            counts['missing_in_ld_reference']+=1;continue
        vector=reference[index].astype(float)
        available=np.isfinite(vector)
        if available.mean()<.95:
            counts['reference_callrate_below_95']+=1;continue
        mean=float(vector[available].mean())
        if min(mean/2,1-mean/2)<.01:
            counts['reference_maf_below_1']+=1;continue
        vector=np.where(available,vector-mean,0)
        norm=float(np.linalg.norm(vector))
        if not norm:
            counts['monomorphic_reference']+=1;continue
        vector/=norm
        active=[a for a in active if a[0]==ch and pos-a[1]<=window]
        r2=max((float(a[2]@vector)**2 for a in active),default=0)
        if r2>threshold:
            counts['ld_removed']+=1;continue
        max_retained_r2=max(max_retained_r2,r2)
        active.append((ch,pos,vector))
        retained.append(r)
    return retained,dict(counts),max_retained_r2


def prepare_panel(threshold=.1,window=1000000):
    freq_file=fetch_public(FREQUENCIES_URL,FREQUENCIES_SHA)
    bed_file=fetch_public(GENOTYPES_URL,GENOTYPES_SHA)
    with zipfile.ZipFile(freq_file) as archive:
        reader=csv.DictReader(io.TextIOWrapper(archive.open('reference_21_P.txt'),encoding='utf-8'),delimiter='\t')
        populations=reader.fieldnames[5:]
        rows=list(reader)
    if len(populations)!=21 or set(populations)!=set(LABELS):
        raise ValueError('Cabeçalho de componentes não corresponde ao snapshot conhecido.')
    with zipfile.ZipFile(bed_file) as archive:
        bim=[line.split() for line in archive.read('Ancestry_2014.bim').decode().splitlines() if line.strip()]
        fam=[line.split() for line in archive.read('Ancestry_2014.fam').decode().splitlines() if line.strip()]
        reference=decode_bed(archive.read('Ancestry_2014.bed'),len(fam),len(bim))
    reference_index={(r[0],int(r[3])):i for i,r in enumerate(bim)}
    for r in rows:
        index=reference_index.get((r['Chromosome'],int(r['Position'])))
        if index is not None and {r['A1'],r['A2']}!={bim[index][4],bim[index][5]}:
            raise ValueError('Alelos dos dois snapshots públicos divergem na mesma posição.')
    say('Verificando LD nos genótipos públicos',len(fam))
    retained,excluded,max_r2=ld_prune(rows,reference,reference_index,threshold,window)
    markers=[]
    for r in retained:
        markers.append(dict(chrom=r['Chromosome'],pos=int(r['Position']),rsid=r['Marker'],
                            effect_allele=r['A2'],other_allele=r['A1'],
                            frequencies={p:float(r[p]) for p in populations}))
    panel={'assembly':'GRCh37','strand':'+','source':f'NHGRI/CRGGH; Baker, Rotimi & Shriner (2017); projeção K=21; {PAPER}',
           'ld_pruned':True,'populations':populations,'population_labels':LABELS,'markers':markers,
           'provenance':{'download_page':PAGE,'frequencies_url':FREQUENCIES_URL,'frequencies_sha256':FREQUENCIES_SHA,
                         'ld_reference_url':GENOTYPES_URL,'ld_reference_sha256':GENOTYPES_SHA,
                         'paper':PAPER,'retrieved':datetime.now(timezone.utc).date().isoformat(),'allele_frequency_orientation':'A2, conforme NHGRI; coordenadas compatíveis com 1000G fase 3 / GRCh37',
                         'original_markers':len(rows),'ld_reference_individuals':len(fam),'ld_reference_markers':len(bim)},
           'quality_control':{'method':'Poda empírica gulosa de dosagens no painel público de 2014, em ordem genômica; SNPs palindrômicos excluídos.',
                              'window_bp':window,'r2_threshold':threshold,'maximum_retained_pair_r2_in_window':max_r2,
                              'reference_maf_min':.01,'reference_call_rate_min':.95,'excluded':excluded,
                              'retained':len(markers),'limitation':'r² em painel agrupado; pode restar LD específico de populações e LD de longo alcance. Painel de 2014 não contém todos os indivíduos do estudo de 2017.'}}
    return panel


def match_genome(genome,panel):
    if genome.build!=panel['assembly']: raise ValueError('Arquivo pessoal e painel usam assemblies diferentes.')
    selected=[]
    counts=Counter()
    for r in panel['markers']:
        call=genome.calls.get((r['chrom'],r['pos']))
        if not call: counts['absent']+=1;continue
        if not call.alleles:counts['no_call']+=1;continue
        if len(call.alleles)!=2:counts['non_diploid']+=1;continue
        if r.get('rsid') and call.rsids and r['rsid'] not in call.rsids:counts['rsid_mismatch']+=1;continue
        if any(a not in (r['effect_allele'],r['other_allele']) for a in call.alleles):counts['allele_mismatch']+=1;continue
        selected.append((r,call.alleles.count(r['effect_allele'])))
    if len(selected)<500 or len({r['chrom'] for r,_ in selected})<5:
        raise ValueError('Sobreposição insuficiente após controle de qualidade.')
    return selected,dict(counts)


def analyze(genome,panel,replicates=200):
    selected,excluded=match_genome(genome,panel)
    pops=panel['populations']
    f=np.asarray([[r['frequencies'][p] for p in pops] for r,_ in selected])
    g=np.asarray([dose for _,dose in selected])
    if np.linalg.matrix_rank(f[:,:-1]-f[:,-1,None],tol=1e-5)<len(pops)-1:
        raise ValueError('Componentes não identificáveis no subconjunto de marcadores disponível.')
    say('Projetando genótipos nos componentes de referência',len(selected))
    fit=fit_mixture(f,g)
    if not fit['converged']:raise ValueError(f"Projeção não convergiu; gap KKT={fit['kkt_gap']}")
    q=fit['q']
    blocks=[(r['chrom'],r['pos']//5000000) for r,_ in selected]
    bootstrap=bootstrap_blocks(f,g,blocks,replicates=replicates,initial=q,progress=say)
    if bootstrap['nonconverged']:
        raise ValueError(f"{bootstrap['nonconverged']} reamostragens não convergiram; revisão numérica necessária.")
    proportions=[]
    for i,p in enumerate(pops):
        proportions.append({'population':p,'label':LABELS[p],'percent':round(float(q[i])*100,4),
                            'lower':round(float(bootstrap['lower'][i])*100,4),'upper':round(float(bootstrap['upper'][i])*100,4),
                            'small_component':bool(q[i]<.01),'boundary_in_resamples':bool(bootstrap['lower'][i]<1e-6)})
    regional=[]
    for name,members in GROUPS.items():
        ix=[pops.index(p) for p in members]
        values=bootstrap['samples'][:,ix].sum(axis=1)
        regional.append({'population':name,'percent':round(float(q[ix].sum())*100,4),
                         'lower':round(float(np.percentile(values,2.5))*100,4),'upper':round(float(np.percentile(values,97.5))*100,4),
                         'members':members})
    # Sensibilidade a regiões cromossômicas: ajustar retirando cada cromossomo.
    chromosome_results=[]
    for ch in sorted({r['chrom'] for r,_ in selected},key=int):
        keep=np.asarray([r['chrom']!=ch for r,_ in selected])
        sensitivity=fit_mixture(f[keep],g[keep],initial=q,tolerance=1e-7)
        if not sensitivity['converged']:raise ValueError('Ajuste sem um cromossomo não convergiu.')
        chromosome_results.append({'excluded_chromosome':ch,'percent':[float(v)*100 for v in sensitivity['q']]})
    # Simulação interna: controle de implementação; não é validação externa.
    rng=np.random.default_rng(20261019)
    synthetic_q=np.zeros(len(pops))
    for p,v in [('Northern European',.35),('Southern European',.30),('Amerindian',.20),('West-Central African',.15)]:synthetic_q[pops.index(p)]=v
    synthetic_g=rng.binomial(2,f@synthetic_q)
    recovered=fit_mixture(f,synthetic_g)
    if not recovered['converged']:raise ValueError('Controle sintético não convergiu.')
    validation={'type':'Simulação interna com o mesmo painel e SNPs disponíveis; não é validação populacional independente.',
                'seed':20261019,'true_percent':dict(zip(pops,(synthetic_q*100).tolist())),
                'estimated_percent':dict(zip(pops,(recovered['q']*100).tolist()))}
    return {'available':True,'source':panel['source'],'assembly':genome.build,'markers':len(selected),
            'chromosomes':len({r['chrom'] for r,_ in selected}),'converged':True,'iterations':fit['iterations'],
            'kkt_gap':fit['kkt_gap'],'log_likelihood_per_marker':fit['log_likelihood_per_marker'],
            'proportions':sorted(proportions,key=lambda p:-p['percent']),
            'regional_proportions':sorted(regional,key=lambda p:-p['percent']),
            'method':'Projeção de máxima verossimilhança binomial em 21 componentes NHGRI, frequências fixas do alelo A2 e soma dos pesos igual a 1.',
            'uncertainty':{'method':'Percentis 2,5–97,5 de reamostragem de blocos físicos de 5 Mb, condicionados às frequências de referência fixas.',
                           'replicates':replicates,'seed':bootstrap['seed'],'blocks':bootstrap['blocks'],'nonconverged':0,
                           'limitation':'Faixas de reamostragem não abrangem incerteza das frequências, vieses do painel, erros sistemáticos de genotipagem, assembly/orientação ou escolha de componentes.'},
            'chromosome_sensitivity':chromosome_results,'synthetic_validation':validation,'matching_excluded':excluded,
            'reference_provenance':panel['provenance'],'reference_quality_control':panel['quality_control'],
            'reference_panel_markers':len(panel['markers']),
            'limitation':'Estimativa exploratória condicionada ao modelo K=21 de 2017, sem calibração para este chip ou validação externa individual. Componentes menores que 1% ou com reamostragens na fronteira zero são instáveis. Não determina etnia, nacionalidade, tribo ou povos históricos; Norte/Sul da Europa não identifica países ou ancestrais individuais.'}


def write_html(result,summary,path):
    e=html.escape
    regional=result['regional_proportions']
    def n(x):return f'{x:.1f}'.replace('.',',')
    def rows(values):
        return ''.join(f'<tr><td>{e(r.get("label",r["population"]))}</td><td><strong>{n(r["percent"])}%</strong></td><td>{n(r["lower"])}–{n(r["upper"])}%</td></tr>' for r in values)
    bars=''.join(f'<div class="barrow"><span>{e(r["population"])}</span><div class="track"><div style="width:{max(0,min(100,r["percent"]))}%"></div></div><strong>{n(r["percent"])}%</strong></div>' for r in regional)
    document=f'''<!doctype html><html lang="pt-BR"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta http-equiv="Content-Security-Policy" content="default-src 'none'; style-src 'unsafe-inline'; base-uri 'none'"><title>Hélice · ancestralidade Genera</title><style>body{{font:15px/1.7 'Segoe UI',sans-serif;color:#18343e;background:#f5f6f2;margin:0}}main{{max-width:1050px;margin:40px auto;padding:0 25px}}h1{{font-size:35px;letter-spacing:-1px;margin-bottom:5px}}h2{{margin-top:35px;font-size:23px}}.eyebrow{{font-size:11px;letter-spacing:2px;color:#18766b}}.card{{background:white;border:1px solid #dce3df;border-radius:10px;padding:26px;margin:20px 0}}.muted,small{{color:#687c80}}.note{{background:#fcf7ec;border:1px solid #e2d4b9;border-radius:8px;padding:20px}}.barrow{{display:flex;align-items:center;gap:16px;margin:18px 0}}.barrow>span{{width:270px;font-size:13px}}.barrow>strong{{width:70px;text-align:right}}.track{{background:#edf1e9;height:13px;border-radius:6px;flex:1;overflow:hidden}}.track>div{{height:100%;background:#18766b}}table{{width:100%;border-collapse:collapse;font-size:13px}}th,td{{padding:12px;border-bottom:1px solid #e1e8df;text-align:left}}th{{background:#edf2eb}}code{{overflow-wrap:anywhere;font-size:11px}}a{{color:#18766b}}li{{margin:7px 0}}@media(max-width:650px){{main{{padding:0 16px}}.barrow{{flex-wrap:wrap}}.barrow>span{{width:100%}}.card{{padding:18px}}}}@media print{{body{{background:white}}main{{margin:0}}.card{{break-inside:avoid}}tr{{break-inside:avoid}}}}</style><main><div class="eyebrow">HÉLICE · ANÁLISE LOCAL · 19 OUT 2026</div><h1>Ancestralidade · arquivo da Genera</h1><p class="muted">Projeção exploratória nos componentes populacionais de referência NHGRI K=21.</p><div class="note">{e(result['limitation'])}</div><section class="card"><h2>Resumo por agrupamentos do modelo</h2><p class="muted">Estes grupos somam componentes do painel; não são categorias étnicas ou percentuais de nacionalidades.</p>{bars}<table><thead><tr><th>Agrupamento</th><th>Estimativa</th><th>Faixa de reamostragem</th></tr></thead><tbody>{rows(regional)}</tbody></table></section><section class="card"><h2>Os 21 componentes originais</h2><p class="muted">Rótulos traduzidos; nomes originais e definição dos agrupamentos estão no JSON. Valores pequenos e próximos de zero não devem ser tomados como ascendência confirmada.</p><table><thead><tr><th>Componente do painel</th><th>Estimativa</th><th>Faixa de reamostragem</th></tr></thead><tbody>{rows(result['proportions'])}</tbody></table></section><section class="card"><h2>Dados e controle de qualidade</h2><ul><li>{result['markers']:,} SNPs utilizáveis em {result['chromosomes']} cromossomos autossômicos, entre {result['reference_panel_markers']:,} marcadores da referência após filtros.</li><li>Referência GRCh37; correspondência exata de posição, rsID e alelos; nenhuma imputação, complementação de fita ou troca de alelos por suposição.</li><li>SNPs A/T e C/G excluídos por ambiguidade de fita. Poda empírica r² ≤ {result['reference_quality_control']['r2_threshold']}, janela de {result['reference_quality_control']['window_bp']//1000} kb, nos {result['reference_provenance']['ld_reference_individuals']} indivíduos públicos de 2014.</li><li>Convergência numérica atingida; gap KKT {result['kkt_gap']:.3g}; {result['iterations']} iterações.</li><li>{result['uncertainty']['replicates']} reamostragens em {result['uncertainty']['blocks']} blocos de 5 Mb, semente fixa; todas convergiram.</li><li>Verificação de sensibilidade retirando cada cromossomo, registrada no JSON.</li></ul><p><strong>Arquivo analisado:</strong> {e(summary['name'])}<br><strong>SHA-256:</strong> <code>{e(summary['sha256'])}</code></p><p>{e(result['uncertainty']['method'])} {e(result['uncertainty']['limitation'])}</p></section><section class="card"><h2>Referências e interpretação</h2><p><a href="{e(PAGE)}">NHGRI / CRGGH · painel público de projeção</a> · <a href="{e(PAPER)}">Baker, Rotimi &amp; Shriner (2017)</a></p><p>O painel original reúne frequências de 21 componentes e 19.075 SNPs. Foi preparado no contexto de um estudo com amostras globais. O modelo e o conjunto de referências influenciam fortemente a resolução; resultados não são diretamente intercambiáveis com o relatório comercial da Genera.</p><p>{e(result['reference_quality_control']['limitation'])}</p><p>Uma simulação interna testa a implementação, mas não avalia acurácia em pessoas reais. Não foram atribuídos países, tribos ou povos antigos a partir de pequenos componentes.</p></section><p class="muted"><small>Relatório privado, gerado neste computador. Nenhum genótipo pessoal foi enviado para obter as referências públicas.</small></p></main></html>'''
    path.write_text(document,encoding='utf-8')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=pathlib.Path)
    parser.add_argument('--prepare-only',action='store_true',help='Baixa e prepara somente a referência pública, sem DNA pessoal.')
    parser.add_argument('--output',type=pathlib.Path,default=ROOT/'reports'/'ancestralidade_genera.json')
    parser.add_argument('--bootstrap',type=int,default=200)
    parser.add_argument('--rebuild-panel',action='store_true')
    args=parser.parse_args()
    if not args.prepare_only and args.input is None:parser.error('Informe --input ou use --prepare-only para preparar somente a referência pública.')
    if not args.prepare_only and args.bootstrap<100:raise ValueError('Use pelo menos 100 reamostragens para o relatório final.')
    if not PANEL_PATH.exists() or args.rebuild_panel:
        panel=prepare_panel()
        PANEL_PATH.write_text(json.dumps(panel,ensure_ascii=False),encoding='utf-8')
    else:panel=json.loads(PANEL_PATH.read_text(encoding='utf-8'))
    say('Painel público preparado',len(panel['markers']))
    if args.prepare_only:
        raise SystemExit(0)
    genome=parse_file(args.input)
    result=analyze(genome,panel,replicates=args.bootstrap)
    result['panel_sha256']=hashlib.sha256(PANEL_PATH.read_bytes()).hexdigest()
    payload={'version':'1.1.0','summary':genome.summary(),'ancestry':result}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,ensure_ascii=False,indent=2),encoding='utf-8')
    write_html(result,payload['summary'],args.output.with_suffix('.html'))
    with args.output.with_suffix('.csv').open('w',encoding='utf-8-sig',newline='') as stream:
        writer=csv.writer(stream,delimiter=';')
        writer.writerow(['component','label','percent','resample_2.5_percent','resample_97.5_percent'])
        for row in result['proportions']:writer.writerow([row['population'],row['label'],row['percent'],row['lower'],row['upper']])
    say('Análise concluída e relatórios privados salvos')
    print(json.dumps({'markers':result['markers'],'regional_proportions':result['regional_proportions']},ensure_ascii=True))
