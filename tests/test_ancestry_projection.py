import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import pathlib
import sys
import unittest
import tempfile
import hashlib
import json
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]))
from ancestry_projection import fit_mixture,project_simplex,bootstrap_blocks
from tools.nhgri_ancestry import decode_bed,ld_prune
import ancestry
from genetics import Genome,ROOT


class ProjectionTests(unittest.TestCase):
    def setUp(self):
        self.rng=np.random.default_rng(58)
        self.f=self.rng.uniform(.05,.95,(1500,21))
        self.q=np.zeros(21);self.q[[1,7,15]]=[.2,.5,.3]
        self.g=2*(self.f@self.q)
    def test_recovers_known_expected_dosage_in_21_components(self):
        result=fit_mixture(self.f,self.g)
        self.assertTrue(result['converged'])
        np.testing.assert_allclose(result['q'],self.q,atol=1e-5)
        self.assertLess(result['kkt_gap'],1e-8)
    def test_allele_inversion_invariant(self):
        a=fit_mixture(self.f,self.g)
        b=fit_mixture(1-self.f,2-self.g)
        np.testing.assert_allclose(a['q'],b['q'],atol=1e-6)
    def test_boundary_solution_and_simplex(self):
        result=fit_mixture(self.f,2*self.f[:,3])
        self.assertTrue(result['converged'])
        self.assertAlmostEqual(result['q'][3],1,places=5)
        self.assertTrue((result['q']>=0).all())
        self.assertAlmostEqual(result['q'].sum(),1,places=12)
    def test_weighted_fit_matches_explicit_duplicates(self):
        f=self.f[:100,:3];g=self.g[:100]
        weights=np.tile([1,2],50)
        a=fit_mixture(f,g,weights=weights)
        b=fit_mixture(np.repeat(f,weights,axis=0),np.repeat(g,weights))
        np.testing.assert_allclose(a['q'],b['q'],atol=1e-6)
    def test_invalid_frequency_rejected(self):
        f=self.f.copy();f[0,0]=float('nan')
        with self.assertRaises(ValueError):fit_mixture(f,self.g)
        with self.assertRaises(ValueError):fit_mixture(self.f,self.g,weights=np.zeros(len(self.g)))
    def test_bootstrap_is_reproducible(self):
        f=self.f[:80,:3];g=self.g[:80];blocks=[str(i//5) for i in range(80)]
        a=bootstrap_blocks(f,g,blocks,replicates=4,seed=4)
        b=bootstrap_blocks(f,g,blocks,replicates=4,seed=4)
        np.testing.assert_array_equal(a['samples'],b['samples'])
    def test_simplex_projection(self):
        q=project_simplex(np.asarray([-2,2,.5]))
        self.assertTrue((q>=0).all());self.assertAlmostEqual(q.sum(),1)


class PublicReferenceTests(unittest.TestCase):
    def test_plink_bed_a2_dosage_and_missing(self):
        values=decode_bed(bytes.fromhex('6c1b01e4'),4,1)[0]
        self.assertEqual(values[0],0)
        self.assertTrue(np.isnan(values[1]))
        self.assertEqual(values[2],1)
        self.assertEqual(values[3],2)
    def test_bed_bad_magic_and_dimensions(self):
        with self.assertRaises(ValueError):decode_bed(b'bad',4,1)
        with self.assertRaises(ValueError):decode_bed(bytes.fromhex('6c1b01'),4,1)
    def test_pruning_removes_correlated_nearby_snp(self):
        rows=[{'Chromosome':'1','Position':str(pos),'A1':'A','A2':'G'} for pos in (10,20,2000000)]
        g=np.asarray([[0,0,1,1,2,2]]*3,dtype=float)
        index={('1',int(row['Position'])):i for i,row in enumerate(rows)}
        kept,excluded,maxr2=ld_prune(rows,g,index,threshold=.1,window=1000)
        self.assertEqual(len(kept),2)
        self.assertEqual(excluded['ld_removed'],1)
        self.assertEqual(maxr2,0)
    def test_palindromic_reference_not_used(self):
        rows=[{'Chromosome':'1','Position':'10','A1':'A','A2':'T'}]
        kept,excluded,_=ld_prune(rows,np.asarray([[0,1,2]],dtype=float),{('1',10):0})
        self.assertFalse(kept);self.assertEqual(excluded['palindromic'],1)


class PersonalCacheTests(unittest.TestCase):
    def test_requires_identical_input_panel_and_version(self):
        workspace=(ROOT/'.runtime').resolve();workspace.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=workspace,prefix='cache-test-') as temp:
            root=pathlib.Path(temp).resolve()
            self.assertTrue(root.is_relative_to(workspace))
            (root/'data/ancestry').mkdir(parents=True)
            (root/'reports').mkdir()
            panel={'assembly':'GRCh37','strand':'+','source':'SYNTHETIC','ld_pruned':True,'populations':['A','B'],'markers':[]}
            path=root/'data/ancestry/nhgri_k21_grch37.json'
            path.write_text(json.dumps(panel),encoding='utf-8')
            result={'available':True,'panel_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'proportions':[{'population':'A','percent':60}]}
            payload={'version':'1.1.0','summary':{'sha256':'SYNTHETIC-HASH'},'ancestry':result}
            (root/'reports/ancestralidade_genera.json').write_text(json.dumps(payload),encoding='utf-8')
            genome=Genome('synthetic','GRCh37','chip','SYNTHETIC-HASH')
            with patch.object(ancestry,'ROOT',root):
                self.assertTrue(ancestry.installed_result(genome)['cache_verified_by_input_sha256'])
                genome.sha256='DIFFERENT'
                self.assertFalse(ancestry.installed_result(genome)['available'])
                genome.sha256='SYNTHETIC-HASH'
                path.write_text(json.dumps(panel)+' ',encoding='utf-8')
                self.assertFalse(ancestry.installed_result(genome)['available'])


if __name__=='__main__':unittest.main()
