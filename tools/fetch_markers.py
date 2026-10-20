"""Baixa somente referências públicas, nunca genótipos pessoais."""
import concurrent.futures
import json
import pathlib
import urllib.request
from datetime import date
import sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from catalog import MARKERS

ROOT = pathlib.Path(__file__).resolve().parents[1]


def fetch(item):
    rsid = item['rsid']
    result = {'rsid': rsid, 'assemblies': {}}
    for build, host in [('GRCh37', 'grch37.rest.ensembl.org'), ('GRCh38', 'rest.ensembl.org')]:
        url = f'https://{host}/variation/human/{rsid}?content-type=application/json'
        with urllib.request.urlopen(url, timeout=60) as response:
            data = json.load(response)
        matches = [m for m in data['mappings'] if m['assembly_name'] == build and m['coord_system'] == 'chromosome' and m['strand'] == 1 and m['start'] == m['end'] and m['seq_region_name'] in [str(n) for n in range(1,23)] + ['X','Y']]
        if len(matches) != 1:
            raise ValueError(f'Mapeamento não único para {rsid} {build}')
        m = matches[0]
        if item['effect_allele'] not in m['allele_string'].split('/'):
            raise ValueError(f'Alelo incompatível: {rsid} {build}')
        result['assemblies'][build] = {'chrom': m['seq_region_name'], 'pos': m['start'], 'alleles': m['allele_string'].split('/'), 'source': url}
    return result


if __name__ == '__main__':
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        records = list(pool.map(fetch, MARKERS))
    path = ROOT / 'data' / 'markers.json'
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps({'retrieved': str(date.today()), 'provider': 'Ensembl REST / dbSNP', 'records': records}, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'{len(records)} marcadores validados em GRCh37 e GRCh38.')
