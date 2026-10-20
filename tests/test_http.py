import json
import pathlib
import sys
import tempfile
import threading
import time
import unittest
import urllib.error
import urllib.request
from functools import partial
from unittest.mock import patch
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]))
import server
from genetics import ROOT
from unittest.mock import MagicMock


class HTTPTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.http=server.ThreadingHTTPServer(('127.0.0.1',0),server.Handler)
        cls.http.daemon_threads=True
        cls.thread=threading.Thread(target=cls.http.serve_forever,daemon=True)
        cls.thread.start()
        cls.base=f'http://127.0.0.1:{cls.http.server_port}'
    @classmethod
    def tearDownClass(cls):
        cls.http.shutdown();cls.http.server_close();cls.thread.join(timeout=3)
    def request(self,path,method='GET',data=None,headers=None):
        req=urllib.request.Request(self.base+path,method=method,data=data,headers=headers or {})
        with urllib.request.urlopen(req,timeout=20) as response:
            return response.status,response.headers,response.read()
    def token(self): return json.loads(self.request('/api/session')[2])['token']
    def test_static_security_headers(self):
        code,headers,data=self.request('/')
        self.assertEqual(code,200)
        self.assertIn('Hélice'.encode('utf-8'),data)
        self.assertIn("connect-src 'self'",headers['Content-Security-Policy'])
        self.assertEqual(headers['Cache-Control'],'no-store')
    def test_token_required(self):
        with self.assertRaises(urllib.error.HTTPError) as error: self.request('/api/report')
        self.assertEqual(error.exception.code,403)
    def test_foreign_origin_blocked(self):
        with self.assertRaises(urllib.error.HTTPError) as error:
            self.request('/api/status',headers={'X-DNA-Token':self.token(),'Origin':'https://example.com'})
        self.assertEqual(error.exception.code,403)
    def test_foreign_host_blocked(self):
        with self.assertRaises(urllib.error.HTTPError) as error:
            self.request('/api/session',headers={'Host':'example.com'})
        self.assertEqual(error.exception.code,403)
    def test_status_reports_missing_ancestry_reference(self):
        root=MagicMock()
        target=root.__truediv__.return_value
        target.__truediv__.return_value=target
        target.exists.return_value=False
        with patch.object(server,'ROOT',root):
            state=json.loads(self.request('/api/status',headers={'X-DNA-Token':self.token()})[2])
        self.assertFalse(state['ancestry_reference']['available'])
    def test_upload_analyze_search_failed_import_and_clear(self):
        workspace=(ROOT/'.runtime').resolve();workspace.mkdir(exist_ok=True)
        headers={'X-DNA-Token':self.token(),'Content-Type':'application/octet-stream'}
        clinical_result={'available':False,'summary':{},'findings':[],'total_findings':0,'message':'Synthetic test'}
        original_mkstemp=tempfile.mkstemp
        with patch.object(server.clinical,'analyze',return_value=clinical_result), patch.object(server.tempfile,'mkstemp',side_effect=partial(original_mkstemp,dir=workspace)):
            payload=(ROOT/'examples/exemplo_sintetico.vcf').read_bytes()
            self.assertEqual(self.request('/api/upload?name=synthetic.vcf',method='POST',data=payload,headers=headers)[0],202)
            for _ in range(200):
                state=json.loads(self.request('/api/status',headers=headers)[2])
                if not state['job']['running']: break
                time.sleep(.02)
            self.assertIsNone(state['job']['error'])
            result=json.loads(self.request('/api/report',headers=headers)[2])
            self.assertEqual(result['summary']['name'],'synthetic.vcf')
            self.assertEqual(result['summary']['called'],3)
            found=json.loads(self.request('/api/search?q=rs4988235',headers=headers)[2])
            self.assertEqual(found['results'][0]['genotype'],'G/A')
            self.request('/api/upload?name=invalid.txt',method='POST',data=b'invalid',headers=headers)
            for _ in range(200):
                state=json.loads(self.request('/api/status',headers=headers)[2])
                if not state['job']['running']: break
                time.sleep(.02)
            self.assertIsNotNone(state['job']['error'])
            self.assertEqual(json.loads(self.request('/api/report',headers=headers)[2])['summary']['name'],'synthetic.vcf')
            self.request('/api/clear',method='POST',headers=headers)
            self.assertFalse(json.loads(self.request('/api/status',headers=headers)[2])['loaded'])
            self.assertFalse(list(workspace.glob('helice-*')))


if __name__=='__main__': unittest.main()
