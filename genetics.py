"""Importação conservadora de genótipos e VCF, sem dependências externas."""
from __future__ import annotations
import csv
import gzip
import hashlib
import io
import json
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from catalog import MARKERS

ROOT = Path(__file__).resolve().parent
BUILDS = ('GRCh37', 'GRCh38')
CHROMS = {str(n) for n in range(1, 23)} | {'X', 'Y', 'MT'}


def chromosome(value):
    value = str(value).strip().upper().removeprefix('CHR')
    return {'23': 'X', '24': 'Y', '25': 'XY', '26': 'MT', 'M': 'MT'}.get(value, value)


def detect_build(text):
    found = set()
    for pattern, build in [(r'grch37|hg19|build\s*37|assembly\s*37', 'GRCh37'),
                           (r'grch38|hg38|build\s*38|assembly\s*38', 'GRCh38')]:
        if re.search(pattern, text, re.I):
            found.add(build)
    if len(found) > 1:
        raise ValueError('Cabeçalho contém referências genômicas conflitantes.')
    return next(iter(found), None)


@dataclass(slots=True)
class Call:
    chrom: str
    pos: int
    rsids: tuple[str, ...]
    alleles: tuple[str, ...] | None
    ref: str = ''
    alts: tuple[str, ...] = ()
    reason: str = ''
    dp: float | None = None
    gq: float | None = None
    phased: bool = False

    def display(self):
        return '|'.join(self.alleles) if self.alleles and self.phased else '/'.join(self.alleles or ()) or 'Sem chamada'


@dataclass
class Genome:
    name: str
    build: str
    format: str
    sha256: str
    calls: dict = field(default_factory=dict)
    by_rsid: dict = field(default_factory=dict)
    stats: Counter = field(default_factory=Counter)
    warnings: list = field(default_factory=list)
    samples: list = field(default_factory=list)
    sample: str = ''
    strand: str = '+'

    def add(self, call):
        key = (call.chrom, call.pos)
        self.stats['records'] += 1
        old = self.calls.get(key)
        if old:
            self.stats['duplicates'] += 1
            old_gt=tuple(sorted(old.alleles)) if old.alleles else None
            new_gt=tuple(sorted(call.alleles)) if call.alleles else None
            if old_gt != new_gt or (old.ref and call.ref and old.ref != call.ref):
                old.alleles = None
                old.reason = 'Duplicatas conflitantes'
                self.stats['conflicts'] += 1
            return
        self.calls[key] = call
        for rsid in call.rsids:
            if rsid in self.by_rsid and self.by_rsid[rsid] != key:
                self.by_rsid[rsid] = None
            else:
                self.by_rsid[rsid] = key

    def summary(self):
        quality = Counter(c.reason for c in self.calls.values() if c.reason)
        called = sum(c.alleles is not None for c in self.calls.values())
        diploid = [c for c in self.calls.values() if c.alleles and len(c.alleles) == 2 and c.chrom.isdigit() and all(len(a) == 1 for a in c.alleles)]
        return dict(name=self.name, build=self.build, format=self.format, sha256=self.sha256,
                    strand=self.strand, sample=self.sample, samples=self.samples,
                    loci=len(self.calls), called=called, no_call=len(self.calls)-called,
                    call_rate=round(100*called/max(1,len(self.calls)),2),
                    heterozygosity=round(100*sum(c.alleles[0] != c.alleles[1] for c in diploid)/max(1,len(diploid)),2),
                    chromosomes=dict(Counter(c.chrom for c in self.calls.values())),
                    stats=dict(self.stats), excluded=dict(quality), warnings=self.warnings)


def parse_file(path, name=None, build='auto', strand='auto', sample=None, progress=lambda *args: None):
    path = Path(path)
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    with path.open('rb') as stream:
        compressed = stream.read(2) == b'\x1f\x8b'
    opener = gzip.open if compressed else open
    with opener(path, 'rt', encoding='utf-8-sig', errors='strict') as stream:
        headers = []
        first = None
        for line in stream:
            if len(line) > 10_000_000:
                raise ValueError('Linha excede o limite de tamanho.')
            if line.startswith('#') or not line.strip():
                headers.append(line)
                if line.startswith('#CHROM\t'):
                    break
                if len(headers) > 100000:
                    raise ValueError('Cabeçalho excessivamente longo.')
            else:
                first = line
                break
        header = ''.join(headers)
        vcf = any(h.startswith('##fileformat=VCF') for h in headers)
        declared = detect_build(header)
        if build != 'auto' and build not in BUILDS:
            raise ValueError('Referência inválida. Escolha GRCh37 ou GRCh38.')
        if declared and build != 'auto' and declared != build:
            raise ValueError('A referência escolhida contradiz o cabeçalho. Corrija o arquivo de origem.')
        assembly = declared if build == 'auto' else build
        if not assembly:
            raise ValueError('Referência genômica desconhecida. Selecione GRCh37 ou GRCh38 conforme a documentação do laboratório.')
        positive = bool(re.search(r'plus\s+strand|positive\s+strand|strand\s*[:=]\s*\+|fita\s+positiva', header, re.I))
        negative = bool(re.search(r'minus\s+strand|negative\s+strand|strand\s*[:=]\s*-|fita\s+negativa', header, re.I))
        if not vcf and negative:
            raise ValueError('Cabeçalho declara fita negativa. Converta para a fita positiva antes de importar.')
        if not vcf and strand not in ('auto', '+'):
            raise ValueError('Converta para a fita positiva antes de importar. Não fazemos complementação automática.')
        if not vcf and strand == 'auto' and not positive:
            raise ValueError('Fita desconhecida. Confirme na documentação que os alelos usam a fita genômica positiva e selecione essa opção.')
        genome = Genome(name or path.name, assembly, 'VCF' if vcf else 'Genotipagem TXT/CSV', digest.hexdigest())
        if not declared:
            genome.warnings.append('Referência informada manualmente; não foi confirmada pelo arquivo.')
        if not vcf:
            genome.warnings.extend(['Genotipagem não equivale a sequenciamento completo. Posições ausentes não foram avaliadas.', 'Cabeçalho convertido não comprova a origem, a referência nem a orientação do arquivo original. Verifique com o laboratório.'])
            if not positive:
                genome.warnings.append('Fita positiva informada manualmente pelo usuário.')
        if vcf:
            cols = next((h.rstrip('\r\n').split('\t') for h in headers if h.startswith('#CHROM\t')), None)
            if not cols or len(cols) < 10:
                raise ValueError('VCF precisa de cabeçalho #CHROM, FORMAT e pelo menos uma amostra com GT.')
            genome.samples = cols[9:]
            if sample and sample not in genome.samples:
                raise ValueError('Amostra solicitada não existe no VCF.')
            if not sample and len(genome.samples) > 1:
                raise ValueError('VCF contém múltiplas amostras. Informe o nome da amostra no campo correspondente.')
            genome.sample = sample or genome.samples[0]
            sample_col = 9 + genome.samples.index(genome.sample)
            for number, line in enumerate(stream, 1):
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) <= sample_col:
                    raise ValueError(f'VCF truncado na linha de dados {number}.')
                ch = chromosome(fields[0])
                if ch not in CHROMS:
                    genome.stats['unsupported_contigs'] += 1
                    continue
                try:
                    pos = int(fields[1])
                    if pos < 1: raise ValueError()
                except ValueError:
                    raise ValueError(f'Posição inválida no VCF, linha de dados {number}.')
                ref = fields[3].upper()
                alts = tuple(fields[4].upper().split(','))
                values = dict(zip(fields[8].split(':'), fields[sample_col].split(':')))
                gt = values.get('GT', '.')
                reason = ''
                dp = gq = None
                try:
                    dp = float(values['DP']) if values.get('DP', '.') != '.' else None
                    gq = float(values['GQ']) if values.get('GQ', '.') != '.' else None
                    if dp is not None and (not __import__('math').isfinite(dp) or dp < 0): raise ValueError()
                    if gq is not None and (not __import__('math').isfinite(gq) or gq < 0): raise ValueError()
                except ValueError:
                    reason = 'Qualidade inválida'
                if fields[6] not in ('PASS', '.'):
                    reason = 'FILTER não aprovado'
                elif (dp is not None and dp < 10) or (gq is not None and gq < 20):
                    reason = 'DP < 10 ou GQ < 20'
                if dp is None or gq is None or fields[6] == '.':
                    genome.stats['quality_incomplete'] += 1
                alleles = None
                if '.' in gt or not gt:
                    reason = reason or 'GT ausente ou parcial'
                else:
                    try:
                        indices = [int(x) for x in re.split(r'[/|]', gt)]
                        choices = (ref,) + alts
                        if len(indices) not in (1, 2) or any(i < 0 or i >= len(choices) for i in indices):
                            raise ValueError()
                        alleles = tuple(choices[i] for i in indices)
                        if any(not re.fullmatch('[ACGT]+', a) for a in alleles):
                            reason = reason or 'Alelo simbólico ou estrutural não suportado'
                        if ch.isdigit() and len(alleles)!=2:
                            reason = reason or 'Ploidia autossômica não diploide'
                    except ValueError:
                        reason = reason or 'GT inválido'
                if reason:
                    alleles = None
                genome.add(Call(ch, pos, tuple(r for r in fields[2].split(';') if r.startswith('rs')), alleles, ref, alts, reason, dp, gq, '|' in gt))
                if number % 100000 == 0: progress('Lendo variantes', number)
            genome.warnings.extend(['VCF: DP mínimo 10 e GQ mínimo 20 quando disponíveis; FILTER diferente de PASS ou ponto é excluído. Qualidade ausente é explicitamente contabilizada.', 'Ausência de variante em VCF não implica genótipo de referência. Blocos gVCF, CNVs, expansões, variantes estruturais e fase entre registros não são interpretados.'])
        else:
            if first is None:
                raise ValueError('Arquivo sem registros de genótipos.')
            delimiter = '\t' if '\t' in first else ';' if ';' in first else ',' if ',' in first else None
            def split(line):
                return next(csv.reader([line], delimiter=delimiter)) if delimiter else line.split()
            first_fields = [x.strip().strip('"') for x in split(first)]
            is_header = any(x.lower() in ('rsid', 'snp', 'marker', 'chromosome', 'cromossomo') for x in first_fields)
            mapping = None
            if is_header:
                aliases = {'rsid': ('rsid','snp','marker','id'), 'chrom': ('chromosome','chrom','chr','cromossomo'), 'pos': ('position','pos','posicao','posição'), 'gt': ('genotype','genotipo','genótipo','result'), 'a1': ('allele1','allele_1','alelo1'), 'a2': ('allele2','allele_2','alelo2')}
                lower = [x.lower() for x in first_fields]
                mapping = {k: next((lower.index(a) for a in names if a in lower), None) for k,names in aliases.items()}
                if any(mapping[k] is None for k in ('rsid','chrom','pos')) or (mapping['gt'] is None and (mapping['a1'] is None or mapping['a2'] is None)):
                    raise ValueError('CSV deve conter rsid, chromosome, position e genotype (ou allele1 e allele2).')
            import itertools
            lines = stream if is_header else itertools.chain([first], stream)
            for number, line in enumerate(lines, 1):
                if not line.strip() or line.startswith('#'): continue
                fields = [x.strip().strip('"') for x in split(line)]
                try:
                    if mapping:
                        rsid, ch, pos = (fields[mapping[k]] for k in ('rsid','chrom','pos'))
                        raw = fields[mapping['gt']] if mapping['gt'] is not None else fields[mapping['a1']] + fields[mapping['a2']]
                    else:
                        rsid, ch, pos, raw = fields[:4]
                    ch = chromosome(ch)
                    pos = int(pos)
                    if pos < 1: raise ValueError()
                except (IndexError, ValueError):
                    raise ValueError(f'Linha de genotipagem inválida: {number}. Importação cancelada para evitar relatório parcial.')
                if ch not in CHROMS:
                    genome.stats['unsupported_contigs'] += 1
                    continue
                gt = raw.upper().replace('/', '').replace('|','').replace(' ', '')
                valid = bool(re.fullmatch('[ACGT]{1,2}', gt))
                reason = '' if valid else 'Sem chamada' if gt in ('--','00','0','NC','NN','N','NA','.','..','') else 'Indel ou alelo não suportado no chip'
                if valid and ch.isdigit() and len(gt)!=2:
                    valid=False
                    reason='Ploidia autossômica não diploide'
                genome.add(Call(ch, pos, (rsid,) if re.fullmatch(r'rs\d+',rsid) else (), tuple(gt) if valid else None, reason=reason))
                if number % 100000 == 0: progress('Lendo marcadores', number)
    if not genome.calls:
        raise ValueError('Nenhuma posição humana suportada encontrada.')
    return genome


def panel_results(genome):
    path = ROOT / 'data' / 'markers.json'
    if not path.exists():
        raise ValueError('Snapshot de coordenadas ausente: execute tools/fetch_markers.py.')
    coordinates = {r['rsid']: r for r in json.loads(path.read_text(encoding='utf-8'))['records']}
    results = []
    for marker in MARKERS:
        m = dict(marker)
        coord = coordinates[m['rsid']]['assemblies'][genome.build]
        call = genome.calls.get((coord['chrom'], coord['pos']))
        m.update(chrom=coord['chrom'], pos=coord['pos'], genotype=None, copies=None, status='unavailable')
        if call and call.alleles and all(len(a) == 1 and a in coord['alleles'] for a in call.alleles):
            # rsID não basta: coordenada e alelos precisam concordar.
            if call.rsids and m['rsid'] not in call.rsids:
                m['interpretation'] = 'Identificador do registro diverge do catálogo. Posição não interpretada.'
            elif call.chrom.isdigit() and len(call.alleles) != 2:
                m['interpretation'] = 'Chamada autossômica não diploide. Posição não interpretada.'
            else:
                count = call.alleles.count(m['effect_allele'])
                m.update(status='observed' if count else 'not_observed', genotype=call.display(), copies=count,
                         interpretation=m['present'] if count else m['absent'])
        else:
            m['interpretation'] = ('Posição excluída: ' + (call.reason or 'alelos incompatíveis com o catálogo') if call else 'Marcador ausente do arquivo. Não é um resultado negativo.')
        m['coordinate_source'] = coord['source']
        for k in ('present','absent'): m.pop(k)
        results.append(m)
    return results


def search_calls(genome, query, limit=100):
    query = query.strip()
    if re.fullmatch(r'rs\d+', query):
        key = genome.by_rsid.get(query)
        calls = [genome.calls[key]] if key else []
    elif re.fullmatch(r'(?:chr)?[\dXYMT]+:\d+(?:-\d+)?', query, re.I):
        ch, region = query.split(':')
        ch = chromosome(ch)
        start, _, end = region.partition('-')
        start, end = int(start), int(end or start)
        if end < start or end-start > 1000000:
            raise ValueError('Use intervalo crescente de até 1 milhão de bases.')
        calls = [c for (chrom,pos), c in genome.calls.items() if chrom == ch and start <= pos <= end][:limit]
    else:
        raise ValueError('Busque rsID ou coordenada, por exemplo rs4988235 ou 2:136608646.')
    return [dict(chrom=c.chrom, pos=c.pos, rsids=c.rsids, genotype=c.display(), ref=c.ref, alts=c.alts, reason=c.reason, dp=c.dp, gq=c.gq, phased=c.phased) for c in calls]
