"""Prepara frequências de referências a partir de VCF público e população por amostra.

Seleção por distância física NÃO substitui LD pruning empírico. Entrada deve ser um
VCF de SNPs autossômicos previamente filtrados e podados por LD (ex.: PLINK).
"""
import argparse
import csv
import gzip
import json
import pathlib
import sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]))
from genetics import chromosome, detect_build


def prepare(vcf,samples,assembly,source,min_alleles=40):
    with open(samples,encoding='utf-8-sig') as stream:
        reader=csv.DictReader(stream,delimiter='\t')
        mapping={row['sample']:row['population'] for row in reader}
    pops=sorted(set(mapping.values()))
    if not 2<=len(pops)<=12: raise ValueError('Escolha 2 a 12 populações de referência.')
    markers=[]
    columns=None
    declared=None
    opener=gzip.open if str(vcf).endswith('.gz') else open
    with opener(vcf,'rt',encoding='utf-8-sig') as stream:
        for line in stream:
            if line.startswith('##'):
                if line.startswith(('##reference=','##assembly=')):
                    detected=detect_build(line)
                    if detected and detected!=assembly: raise ValueError('Assembly do VCF incompatível com o solicitado.')
                    declared=detected or declared
                continue
            if line.startswith('#CHROM'):
                names=line.rstrip().split('\t')[9:]
                columns={p:[i+9 for i,name in enumerate(names) if mapping.get(name)==p] for p in pops}
                if any(not values for values in columns.values()): raise ValueError('Alguma população não tem amostras no VCF.')
                continue
            if line.startswith('#') or not line.strip(): continue
            if columns is None: raise ValueError('VCF sem cabeçalho de amostras.')
            fields=line.rstrip().split('\t')
            ch=chromosome(fields[0])
            if ch not in [str(n) for n in range(1,23)] or len(fields[3])!=1 or len(fields[4])!=1 or fields[3] not in 'ACGT' or fields[4] not in 'ACGT' or fields[6]!='PASS': continue
            fmt=fields[8].split(':')
            if 'GT' not in fmt: continue
            gi=fmt.index('GT')
            frequencies={}
            for pop,indices in columns.items():
                alt_count=total=0
                for i in indices:
                    values=fields[i].split(':')
                    if len(values)<=gi: continue
                    gt=values[gi].replace('|','/').split('/')
                    if len(gt)!=2 or any(a not in ('0','1') for a in gt): continue
                    alt_count+=gt.count('1')
                    total+=2
                if total<min_alleles: break
                frequencies[pop]=alt_count/total
            if len(frequencies)!=len(pops) or max(frequencies.values())-min(frequencies.values())<0.05: continue
            markers.append({'chrom':ch,'pos':int(fields[1]),'other_allele':fields[3],'effect_allele':fields[4],'frequencies':frequencies})
    if not declared: raise ValueError('VCF deve declarar GRCh37 ou GRCh38 em ##reference ou ##assembly.')
    return {'assembly':assembly,'strand':'+','source':source,'ld_pruned':True,'populations':pops,'markers':markers}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--vcf',required=True,type=pathlib.Path)
    p.add_argument('--samples',required=True,type=pathlib.Path,help='TSV com cabeçalho sample e population')
    p.add_argument('--assembly',required=True,choices=['GRCh37','GRCh38'])
    p.add_argument('--source',required=True)
    p.add_argument('--output',required=True,type=pathlib.Path)
    p.add_argument('--confirm-ld-pruned',action='store_true',required=True,help='Confirma que o VCF de referência já passou por poda de LD empírica')
    p.add_argument('--min-alleles',type=int,default=40)
    args=p.parse_args()
    panel=prepare(args.vcf,args.samples,args.assembly,args.source,args.min_alleles)
    args.output.write_text(json.dumps(panel,ensure_ascii=False),encoding='utf-8')
    print(f"Painel preparado: {len(panel['markers'])} SNPs. Isso não é uma estimativa de ancestralidade pessoal.")
