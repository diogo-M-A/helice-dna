"""Relatório JSON local por linha de comando; não faz uploads."""
import argparse
import json
import pathlib
import sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]))
from genetics import parse_file
from server import make_report

if __name__=='__main__':
    parser=argparse.ArgumentParser(description='Analisa DNA localmente e exporta JSON')
    parser.add_argument('input',type=pathlib.Path)
    parser.add_argument('--output',type=pathlib.Path,required=True)
    parser.add_argument('--build',default='auto',choices=['auto','GRCh37','GRCh38'])
    parser.add_argument('--strand',default='auto',choices=['auto','+'])
    parser.add_argument('--sample')
    args=parser.parse_args()
    genome=parse_file(args.input,build=args.build,strand=args.strand,sample=args.sample)
    report=make_report(genome)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(f"Relatório local salvo; {len(genome.calls)} posições; {sum(m['status']!='unavailable' for m in report['panel'])} marcadores do painel avaliáveis.")
