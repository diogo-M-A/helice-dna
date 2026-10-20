"""Mistura exploratória por máxima verossimilhança com frequências fornecidas.

Não entrega etnia, parentesco, haplogrupos nem proporções de povos históricos.
"""
import json
import math
import hashlib
from genetics import chromosome, ROOT


def installed_result(genome,progress=lambda *args:None):
    """Referência pública instalada; cache pessoal só se o SHA do arquivo coincidir."""
    path=ROOT/'data'/'ancestry'/'nhgri_k21_grch37.json'
    if not path.exists():
        return {'available':False,'message':'Nenhum painel de referência instalado. Percentuais não foram estimados.'}
    panel=json.loads(path.read_text(encoding='utf-8'))
    if panel.get('assembly')!=genome.build:
        return {'available':False,'message':'O painel NHGRI instalado é GRCh37 e não corresponde à referência deste arquivo.'}
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    cached=ROOT/'reports'/'ancestralidade_genera.json'
    if cached.exists():
        payload=json.loads(cached.read_text(encoding='utf-8'))
        if payload.get('version')=='1.1.0' and payload.get('summary',{}).get('sha256')==genome.sha256 and payload.get('ancestry',{}).get('panel_sha256')==digest:
            result=payload['ancestry']
            result['cache_verified_by_input_sha256']=True
            return result
    progress('Projetando no painel populacional NHGRI',0)
    try:
        result=estimate(genome,panel)
    except (ValueError,ImportError) as exc:
        return {'available':False,'message':str(exc)}
    result['reference_provenance']=panel.get('provenance',{})
    result['reference_quality_control']=panel.get('quality_control',{})
    result['panel_sha256']=digest
    for row in result['proportions']:
        row['label']=panel.get('population_labels',{}).get(row['population'],row['population'])
    return result


def estimate(genome, panel):
    if panel.get('assembly') != genome.build or panel.get('strand') != '+':
        raise ValueError('Painel ancestral deve ter o mesmo assembly e usar strand +.')
    if not panel.get('source') or not panel.get('ld_pruned'):
        raise ValueError('Painel deve declarar source e ld_pruned=true. Use referências independentes com proveniência conhecida.')
    populations=panel.get('populations',[])
    if not 2<=len(populations)<=32 or len(set(populations)) != len(populations) or not all(isinstance(p,str) and p for p in populations):
        raise ValueError('São necessários 2 a 32 componentes distintos.')
    rows=[]
    chroms=set()
    seen=set()
    for row in panel.get('markers',[]):
        ch=chromosome(row['chrom'])
        pos=int(row['pos'])
        key=(ch,pos)
        if key in seen: raise ValueError('Painel ancestral tem posições duplicadas.')
        seen.add(key)
        if ch not in [str(n) for n in range(1,23)]: continue
        effect=row['effect_allele']
        other=row['other_allele']
        if effect not in 'ACGT' or other not in 'ACGT' or len(effect)!=1 or len(other)!=1 or effect==other:
            raise ValueError('Alelos inválidos no painel ancestral.')
        freq=[float(row['frequencies'][p]) for p in populations]
        if any(not math.isfinite(p) or not 0<=p<=1 for p in freq): raise ValueError('Frequências devem estar entre 0 e 1.')
        call=genome.calls.get(key)
        if not call or not call.alleles or len(call.alleles)!=2 or any(a not in (effect,other) for a in call.alleles): continue
        if row.get('rsid') and call.rsids and row['rsid'] not in call.rsids: continue
        if max(freq)-min(freq)<0.05: continue
        rows.append((call.alleles.count(effect),freq))
        chroms.add(ch)
    if len(rows)<500 or len(chroms)<5:
        raise ValueError(f'Cobertura insuficiente: {len(rows)} SNPs informativos em {len(chroms)} cromossomos; mínimo 500 e 5.')
    # Identificabilidade básica: referências idênticas não podem ser separadas.
    for a in range(len(populations)):
        for b in range(a):
            if sum((f[a]-f[b])**2 for _,f in rows)/len(rows)<0.0001:
                raise ValueError('Populações de referência indistinguíveis nos marcadores disponíveis.')
    k=len(populations)
    # As diferenças de frequências precisam ter posto completo sob soma(q)=1.
    # Referências que são combinações exatas de outras não são identificáveis.
    dimension=k-1
    gram=[[sum((f[a]-f[-1])*(f[b]-f[-1]) for _,f in rows)/len(rows) for b in range(dimension)] for a in range(dimension)]
    for pivot in range(dimension):
        largest=max(range(pivot,dimension),key=lambda i:abs(gram[i][pivot]))
        gram[pivot],gram[largest]=gram[largest],gram[pivot]
        if abs(gram[pivot][pivot])<1e-6:
            raise ValueError('Painel não identifica uma mistura única: populações linearmente dependentes ou pouco diferenciadas.')
        for row in range(pivot+1,dimension):
            factor=gram[row][pivot]/gram[pivot][pivot]
            for col in range(pivot,dimension): gram[row][col]-=factor*gram[pivot][col]
    if k>12:
        try:
            from ancestry_projection import fit_mixture
        except ImportError:
            raise ValueError('Painéis com mais de 12 componentes requerem NumPy. Instale a dependência opcional numpy>=1.26 no Python usado para abrir o app.')
        fit=fit_mixture([f for _,f in rows],[g for g,_ in rows])
        return dict(available=True,source=panel['source'],assembly=genome.build,markers=len(rows),chromosomes=len(chroms),
                    converged=fit['converged'],iterations=fit['iterations'],kkt_gap=fit['kkt_gap'],
                    proportions=[{'population':p,'percent':round(100*float(v),2)} for p,v in zip(populations,fit['q'])],
                    method='Projeção de máxima verossimilhança binomial; frequências de componentes fixas; otimização restrita ao simplex.',
                    limitation='Estimativa exploratória condicionada ao painel. Não determina etnia, nacionalidade ou povos históricos; componentes pequenos podem ser instáveis. Sem intervalos de reamostragem nesta execução da interface.')
    q=[1/k]*k
    converged=False
    for iteration in range(2000):
        expected=[0.0]*k
        for dosage,freq in rows:
            p=min(1-1e-9,max(1e-9,sum(qi*fi for qi,fi in zip(q,freq))))
            for j in range(k):
                expected[j]+=dosage*q[j]*freq[j]/p+(2-dosage)*q[j]*(1-freq[j])/(1-p)
        updated=[v/(2*len(rows)) for v in expected]
        if max(abs(a-b) for a,b in zip(q,updated))<1e-7:
            q=updated
            converged=True
            break
        q=updated
    return dict(available=True,source=panel['source'],assembly=genome.build,markers=len(rows),chromosomes=len(chroms),
                converged=converged,iterations=iteration+1,
                proportions=[{'population':p,'percent':round(100*v,2)} for p,v in zip(populations,q)],
                method='Máxima verossimilhança binomial; frequências populacionais fixas; SNPs independentes declarados pelo fornecedor.',
                limitation='Estimativa exploratória sem validação clínica ou intervalos de confiança. Depende da seleção de populações e do controle de LD do painel. Não determina etnia, nacionalidade ou proporção de povos históricos.')
