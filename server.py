"""Hélice: aplicação local de exploração genética. Python >= 3.11."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import secrets
import tempfile
import threading
import traceback
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
from urllib.parse import urlsplit, parse_qs
import webbrowser
import clinical
from genetics import ROOT, parse_file, panel_results, search_calls
from ancestry import estimate,installed_result

TOKEN=secrets.token_urlsafe(32)
LOCK=threading.RLock()
STATE={'genome':None,'report':None,'job':{'running':False,'message':'Pronto para importar','count':0,'error':None}}
MAX_UPLOAD=512*1024*1024


def progress(message,count=0):
    with LOCK:
        STATE['job'].update(message=message,count=count)


def make_report(genome):
    progress('Preparando painel de marcadores')
    report={'version':'1.1.0','summary':genome.summary(),'panel':panel_results(genome),
            'clinical':clinical.analyze(genome,progress),
            'ancestry':installed_result(genome,progress),
            'limitations':[
                'Ferramenta exploratória, sem validação clínica. Achados precisam de confirmação por teste clínico e interpretação profissional.',
                'Variantes raras em chips podem ser falsos positivos. Uma classificação bem revisada no ClinVar não valida a chamada do seu arquivo.',
                'Não calcula risco absoluto nem exclui doenças. Um marcador não observado se refere apenas àquela posição.',
                'BRCA1/BRCA2, câncer hereditário e distúrbios do desenvolvimento exigem avaliação de variantes e tipos de alteração que podem estar ausentes do chip.',
                'Não determina diplótipos completos, HLA, CYP2D6/CNVs, ancestralidade histórica ou resposta individual a dieta e exercício.',
                'Sem imputação, alinhamento de FASTQ/BAM/CRAM, normalização de indels com FASTA ou cálculo de escores poligênicos clínicos.'
            ]}
    return report


def start_job(function):
    with LOCK:
        if STATE['job']['running']:
            raise ValueError('Uma operação está em andamento. Aguarde a conclusão.')
        STATE['job']={'running':True,'message':'Iniciando','count':0,'error':None}
    def work():
        try:
            function()
            with LOCK: STATE['job'].update(running=False,message='Concluído',error=None)
        except Exception as exc:
            # Genótipos e caminhos pessoais não são registrados em logs.
            with LOCK: STATE['job'].update(running=False,message='Operação não concluída',error=str(exc))
    threading.Thread(target=work,daemon=True).start()


def load_genome(path,name=None,build='auto',strand='auto',sample=None,temporary=False):
    try:
        genome=parse_file(path,name,build,strand,sample,progress)
        report=make_report(genome)
        with LOCK:
            STATE['genome']=genome
            STATE['report']=report
    finally:
        if temporary: Path(path).unlink(missing_ok=True)


class Handler(BaseHTTPRequestHandler):
    server_version='HeliceLocal/1.0'

    def log_message(self,*args): pass

    def guard(self,api=False):
        host=self.headers.get('Host','')
        allowed={f'127.0.0.1:{self.server.server_port}',f'localhost:{self.server.server_port}'}
        if host not in allowed: raise PermissionError('Host não permitido.')
        origin=self.headers.get('Origin')
        if origin and origin not in {f'http://{h}' for h in allowed}:
            raise PermissionError('Origem não permitida.')
        if api and not secrets.compare_digest(self.headers.get('X-DNA-Token',''),TOKEN):
            raise PermissionError('Sessão local inválida. Recarregue a página.')

    def send(self,status,payload,kind='application/json; charset=utf-8',extra=None):
        data=json.dumps(payload,ensure_ascii=False).encode() if kind.startswith('application/json') else payload if isinstance(payload,bytes) else payload.encode('utf-8')
        self.send_response(status)
        self.send_header('Content-Type',kind)
        self.send_header('Content-Length',str(len(data)))
        self.send_header('Cache-Control','no-store')
        self.send_header('X-Content-Type-Options','nosniff')
        self.send_header('Referrer-Policy','no-referrer')
        self.send_header('Content-Security-Policy',"default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; connect-src 'self'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'; form-action 'self'")
        if extra:
            for k,v in extra.items(): self.send_header(k,v)
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        try:
            path=urlsplit(self.path).path
            self.guard(api=path.startswith('/api/') and path!='/api/session')
            params=parse_qs(urlsplit(self.path).query)
            if path=='/api/session':
                self.send(200,{'token':TOKEN,'version':'1.1.0'})
            elif path=='/api/status':
                with LOCK:
                    self.send(200,{'job':dict(STATE['job']),'loaded':STATE['report'] is not None,'databases':[clinical.info(b) for b in ('GRCh37','GRCh38')],
                                   'ancestry_reference':{'available':(ROOT/'data'/'ancestry'/'nhgri_k21_grch37.json').exists(),'assembly':'GRCh37','components':21}})
            elif path=='/api/report':
                with LOCK: report=STATE['report']
                if not report: raise ValueError('Importe um arquivo primeiro.')
                self.send(200,report)
            elif path=='/api/search':
                with LOCK: genome=STATE['genome']
                if not genome: raise ValueError('Importe um arquivo primeiro.')
                self.send(200,{'results':search_calls(genome,params.get('q',[''])[0])})
            elif path in ('/','/index.html','/app.js','/style.css'):
                file=ROOT/'web'/('index.html' if path in ('/','/index.html') else path[1:])
                kind='text/html; charset=utf-8' if file.suffix=='.html' else 'text/javascript; charset=utf-8' if file.suffix=='.js' else 'text/css; charset=utf-8'
                self.send(200,file.read_bytes(),kind)
            elif path=='/favicon.ico': self.send(204,b'','image/x-icon')
            else: self.send(404,{'error':'Não encontrado.'})
        except PermissionError as exc: self.send(403,{'error':str(exc)})
        except (ValueError,KeyError) as exc: self.send(400,{'error':str(exc)})
        except (BrokenPipeError,ConnectionResetError): pass
        except Exception: self.send(500,{'error':'Erro interno ao ler os dados.'})

    def do_POST(self):
        tmp=None
        try:
            self.guard(api=True)
            parsed=urlsplit(self.path)
            params={k:v[0] for k,v in parse_qs(parsed.query).items()}
            with LOCK:
                if STATE['job']['running']: raise ValueError('Uma operação está em andamento.')
            if parsed.path=='/api/upload':
                length=int(self.headers.get('Content-Length','0'))
                if not 0<length<=MAX_UPLOAD: raise ValueError('Arquivo vazio ou maior que 512 MB. Para VCF maior use a opção --input no terminal.')
                name=Path(params.get('name','arquivo.txt')).name
                # Pasta temporária do sistema; apagada assim que a operação termina.
                fd,tmp=tempfile.mkstemp(prefix='helice-')
                with os.fdopen(fd,'wb') as stream:
                    remaining=length
                    self.connection.settimeout(120)
                    while remaining:
                        chunk=self.rfile.read(min(1024*1024,remaining))
                        if not chunk: raise ValueError('Upload incompleto.')
                        stream.write(chunk)
                        remaining-=len(chunk)
                kind=params.get('kind','dna')
                retained=tmp
                if kind=='dna':
                    start_job(lambda:load_genome(retained,name,params.get('build','auto'),params.get('strand','auto'),params.get('sample') or None,True))
                elif kind=='clinvar':
                    build=params.get('build','GRCh37')
                    def import_reference():
                        try:
                            clinical.import_clinvar(retained,build,progress,source=name)
                            with LOCK: genome=STATE['genome']
                            if genome:
                                report=make_report(genome)
                                with LOCK: STATE['report']=report
                        finally: Path(retained).unlink(missing_ok=True)
                    start_job(import_reference)
                elif kind=='ancestry':
                    with LOCK: genome=STATE['genome']
                    if not genome: raise ValueError('Importe seu DNA antes do painel.')
                    def import_panel():
                        try:
                            if Path(retained).stat().st_size>50*1024*1024: raise ValueError('Painel ancestral excede 50 MB.')
                            panel=json.loads(Path(retained).read_text(encoding='utf-8-sig'))
                            result=estimate(genome,panel)
                            with LOCK: STATE['report']['ancestry']=result
                        finally: Path(retained).unlink(missing_ok=True)
                    start_job(import_panel)
                else: raise ValueError('Tipo de importação inválido.')
                tmp=None
                self.send(202,{'ok':True})
            elif parsed.path=='/api/download-clinvar':
                build=params.get('build','GRCh37')
                clinical.database_path(build)
                def install():
                    clinical.download(build,progress)
                    with LOCK: genome=STATE['genome']
                    if genome:
                        report=make_report(genome)
                        with LOCK: STATE['report']=report
                start_job(install)
                self.send(202,{'ok':True})
            elif parsed.path=='/api/clear':
                with LOCK:
                    STATE['genome']=None
                    STATE['report']=None
                self.send(200,{'ok':True})
            else: self.send(404,{'error':'Não encontrado.'})
        except PermissionError as exc: self.send(403,{'error':str(exc)})
        except (ValueError,KeyError) as exc: self.send(400,{'error':str(exc)})
        except (BrokenPipeError,ConnectionResetError): pass
        except Exception: self.send(500,{'error':'Falha ao receber o arquivo. Confira o formato e tente novamente.'})
        finally:
            if tmp: Path(tmp).unlink(missing_ok=True)


def main():
    parser=argparse.ArgumentParser(description='Hélice — análise genética local')
    parser.add_argument('--port',type=int,default=8765)
    parser.add_argument('--input',type=Path)
    parser.add_argument('--build',default='auto',choices=['auto','GRCh37','GRCh38'])
    parser.add_argument('--strand',default='auto',choices=['auto','+'])
    parser.add_argument('--sample')
    parser.add_argument('--no-browser',action='store_true')
    args=parser.parse_args()
    server=ThreadingHTTPServer(('127.0.0.1',args.port),Handler)
    server.daemon_threads=True
    if args.input:
        start_job(lambda:load_genome(args.input,build=args.build,strand=args.strand,sample=args.sample))
    url=f'http://127.0.0.1:{server.server_port}'
    print(f'Hélice local disponível em {url}',flush=True)
    print('Processamento local. Ctrl+C encerra e libera os genótipos da sessão.',flush=True)
    if not args.no_browser: webbrowser.open(url)
    try: server.serve_forever()
    except KeyboardInterrupt: pass
    finally:
        server.server_close()
        with LOCK: STATE.update(genome=None,report=None)


if __name__=='__main__': main()
