"""ClinVar VCF local: correspondência exata de assembly + posição + REF + ALT.

Não interpreta causalidade, modo de herança, penetrância, pares em trans ou diagnósticos.
"""
import gzip
import hashlib
import json
import os
import sqlite3
import time
import urllib.request
from contextlib import closing
from pathlib import Path
from urllib.parse import unquote
from genetics import ROOT, chromosome, BUILDS, detect_build

DATA = ROOT / 'data'


def database_path(build):
    if build not in BUILDS: raise ValueError('Assembly inválido.')
    return DATA / f'clinvar-{build}.sqlite'


def info(build):
    path = database_path(build)
    if not path.exists(): return {'available': False, 'build': build}
    with closing(sqlite3.connect(path)) as db:
        metadata = dict(db.execute('SELECT key,value FROM metadata'))
    return {'available': True, 'build': build, **metadata}


def decode_info(raw):
    return {k: unquote(v) for field in raw.split(';') if '=' in field for k, v in [field.split('=', 1)]}


def import_clinvar(path, build, progress=lambda *args: None, source='Arquivo ClinVar local'):
    DATA.mkdir(exist_ok=True)
    target = database_path(build)
    tmp = target.with_suffix('.building.sqlite')
    tmp.unlink(missing_ok=True)
    count = skipped = 0
    date = 'não informada'
    sha = hashlib.sha256()
    with open(path,'rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024),b''): sha.update(chunk)
    with open(path,'rb') as stream: compressed = stream.read(2) == b'\x1f\x8b'
    opener = gzip.open if compressed else open
    db = sqlite3.connect(tmp)
    try:
        db.executescript('PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF; CREATE TABLE variants(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,variation TEXT,details TEXT); CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT);')
        batch = []
        header_seen = False
        assembly_seen = False
        with opener(path,'rt',encoding='utf-8-sig') as stream:
            for line in stream:
                if line.startswith('##'):
                    if line.startswith('##fileDate='): date=line.strip().split('=',1)[1]
                    # NCBI reference headers have GRCh37 / GRCh38 in the FASTA name.
                    if line.startswith(('##reference=', '##assembly=')):
                        declared = detect_build(line)
                        if declared:
                            assembly_seen = True
                            if declared != build: raise ValueError('ClinVar e genótipos usam assemblies diferentes.')
                    continue
                if line.startswith('#CHROM\t'):
                    header_seen = True
                    continue
                if line.startswith('#') or not line.strip(): continue
                if not header_seen: raise ValueError('Base precisa ser um VCF oficial do ClinVar com #CHROM.')
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) < 8: raise ValueError('Registro ClinVar truncado.')
                ch,pos,variation,ref,alt,_,_,raw = fields[:8]
                ch=chromosome(ch)
                annotation=decode_info(raw)
                if ch not in [str(n) for n in range(1,23)]+['X','Y','MT'] or ',' in alt or not ref or not alt or any(c not in 'ACGT' for c in ref+alt):
                    skipped+=1
                    continue
                # Multi-ALT annotations and symbolic alleles are deliberately excluded.
                if not annotation.get('CLNSIG') or not variation.isdigit():
                    skipped+=1
                    continue
                details={
                    'classification':annotation.get('CLNSIG',''),
                    'review':annotation.get('CLNREVSTAT',''),
                    'conditions':annotation.get('CLNDN','').replace('_',' '),
                    'gene':annotation.get('GENEINFO',''),
                    'conflict':annotation.get('CLNSIGCONF',''),
                    'origin':annotation.get('ORIGIN',''),
                    'rsid': 'rs'+annotation['RS'] if annotation.get('RS','').isdigit() else '',
                }
                batch.append((ch,int(pos),ref,alt,variation,json.dumps(details,ensure_ascii=False)))
                count+=1
                if len(batch)>=10000:
                    db.executemany('INSERT INTO variants VALUES (?,?,?,?,?,?)',batch)
                    batch=[]
                if count%100000 == 0: progress('Indexando ClinVar',count)
        if not assembly_seen: raise ValueError('Base sem declaração de referência GRCh37/GRCh38 reconhecível.')
        if not count: raise ValueError('Nenhuma variante ClinVar compatível.')
        db.executemany('INSERT INTO variants VALUES (?,?,?,?,?,?)',batch)
        progress('Criando índice local',count)
        db.execute('CREATE INDEX locus ON variants(chrom,pos)')
        metadata={'source':source,'file_date':date,'imported':time.strftime('%Y-%m-%d'), 'records':str(count),'skipped':str(skipped),'sha256':sha.hexdigest(),'assembly':build}
        db.executemany('INSERT INTO metadata VALUES (?,?)',metadata.items())
        db.commit()
        db.close()
        os.replace(tmp,target)
        return info(build)
    except BaseException:
        db.close()
        tmp.unlink(missing_ok=True)
        raise


def download(build, progress=lambda *args: None):
    database_path(build)  # valida antes de construir URL
    DATA.mkdir(exist_ok=True)
    url=f'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_{build}/clinvar.vcf.gz'
    target=DATA/f'clinvar-{build}.vcf.gz'
    partial=target.with_suffix('.part')
    try:
        req=urllib.request.Request(url,headers={'User-Agent':'HeliceLocal/1.0 (public ClinVar download)'})
        with urllib.request.urlopen(req,timeout=120) as response, partial.open('wb') as stream:
            total=int(response.headers.get('Content-Length',0))
            written=0
            while chunk:=response.read(1024*1024):
                stream.write(chunk)
                written+=len(chunk)
                progress(f'Baixando referência pública • {written//1048576} MB' + (f' / {total//1048576} MB' if total else ''),written)
        os.replace(partial,target)
        return import_clinvar(target,build,progress,source=url)
    finally:
        partial.unlink(missing_ok=True)


def severity(details):
    text=details['classification'].lower()
    if 'conflict' in text or details.get('conflict'):
        return 'conflicting'
    if text in ('pathogenic','likely_pathogenic','pathogenic/likely_pathogenic'):
        return 'pathogenic'
    if 'uncertain' in text:
        return 'uncertain'
    if text in ('benign','likely_benign','benign/likely_benign'):
        return 'benign'
    return 'other'


def stars(review):
    if 'practice_guideline' in review: return 4
    if 'reviewed_by_expert_panel' in review: return 3
    if 'multiple_submitters' in review and 'no_conflicts' in review: return 2
    if 'criteria_provided' in review and 'no_assertion' not in review: return 1
    return 0


def analyze(genome, progress=lambda *args: None):
    metadata=info(genome.build)
    if not metadata['available']:
        return dict(available=False, summary={}, findings=[], total_findings=0,
                    message='Instale a referência pública ClinVar para pesquisar variantes em todos os genes da base.')
    counts={k:0 for k in ('pathogenic','conflicting','uncertain','benign','other')}
    findings=[]
    matched_loci=set()
    public_at_locus=0
    unsupported_indels=0
    limit=20000
    # Percorre o índice público uma vez: evita centenas de milhares de consultas SQL.
    with closing(sqlite3.connect(database_path(genome.build))) as db:
        for i,(ch,pos,ref,alt,variation,raw) in enumerate(db.execute('SELECT chrom,pos,ref,alt,variation,details FROM variants'),1):
            if i%200000==0: progress('Cruzando variantes com ClinVar',i)
            call=genome.calls.get((ch,pos))
            if not call or not call.alleles: continue
            public_at_locus+=1
            if genome.format=='VCF':
                if call.ref != ref or alt not in call.alts or alt not in call.alleles:
                    continue
            else:
                if len(ref)!=1 or len(alt)!=1:
                    unsupported_indels+=1
                    continue
                if alt not in call.alleles or any(a not in (ref,alt) for a in call.alleles): continue
                if call.rsids:
                    rsid=json.loads(raw).get('rsid')
                    if rsid and rsid not in call.rsids: continue
            details=json.loads(raw)
            group=severity(details)
            counts[group]+=1
            matched_loci.add((ch,pos))
            finding=dict(**details,group=group,chrom=ch,pos=pos,ref=ref,alt=alt,
                         genotype=call.display(),copies=call.alleles.count(alt),variation_id=variation,
                         stars=stars(details['review']),source=f'https://www.ncbi.nlm.nih.gov/clinvar/variation/{variation}/',
                         quality='DP/GQ presentes' if call.dp is not None and call.gq is not None else 'Qualidade de chamada não verificável integralmente')
            # Conservar achados prioritários mesmo depois do limite de associações comuns.
            if group in ('pathogenic','conflicting') or len(findings)<limit:
                findings.append(finding)
    order={'pathogenic':0,'conflicting':1,'uncertain':2,'other':3,'benign':4}
    findings.sort(key=lambda f:(order[f['group']],-f['stars'],f['gene']))
    return dict(available=True,metadata=metadata,summary=counts,findings=findings,
                total_findings=sum(counts.values()),matched_loci=len(matched_loci),public_at_locus=public_at_locus,
                omitted=max(0,sum(counts.values())-len(findings)),unsupported_chip_indels=unsupported_indels,
                message='Correspondências de variantes observadas, não diagnósticos. Classificação agregada da base; condições e herança exigem revisão individual. Ausência de achados não exclui doença. Indels exigem a mesma representação exata; nenhuma normalização com FASTA foi aplicada.')
