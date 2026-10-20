import gzip
import json
import pathlib
import random
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]))
from genetics import Call, Genome, parse_file, panel_results, search_calls, ROOT
import clinical
from ancestry import estimate

HEADER='# reference=GRCh37\n# strand=+\n'
VCF='##fileformat=VCFv4.2\n##reference=GRCh37\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTEST\n'


def workspace_temp():
    root=(ROOT/'.runtime').resolve()
    root.mkdir(exist_ok=True)
    result=tempfile.TemporaryDirectory(dir=root,prefix='tests-')
    assert pathlib.Path(result.name).resolve().is_relative_to(root)
    return result


class ParsingTests(unittest.TestCase):
    def setUp(self):
        self.temp=workspace_temp()
        self.path=pathlib.Path(self.temp.name)/'input.txt'
    def tearDown(self): self.temp.cleanup()
    def parse(self,text,**kwargs):
        self.path.write_text(text,encoding='utf-8')
        return parse_file(self.path,**kwargs)
    def test_chip_example(self):
        g=parse_file(ROOT/'examples/exemplo_sintetico.txt')
        p={m['rsid']:m for m in panel_results(g)}
        self.assertEqual(g.build,'GRCh37')
        self.assertEqual(p['rs4988235']['copies'],1)
        self.assertEqual(p['rs334']['status'],'not_observed')
        self.assertEqual(p['rs3918290']['status'],'unavailable')
    def test_csv_header_and_separate_alleles(self):
        g=self.parse(HEADER+'rsid;chromosome;position;allele1;allele2\nrs1;chr1;10;A;G\n')
        self.assertEqual(g.calls[('1',10)].alleles,('A','G'))
    def test_unknown_build_rejected(self):
        with self.assertRaisesRegex(ValueError,'Referência genômica desconhecida'):
            self.parse('# strand=+\nrs1 1 10 AG\n')
    def test_manual_build_and_strand(self):
        g=self.parse('rs1 1 10 AG\n',build='GRCh37',strand='+')
        self.assertEqual(len(g.warnings),4)
    def test_build_contradiction(self):
        with self.assertRaisesRegex(ValueError,'contradiz'):
            self.parse(HEADER+'rs1 1 10 AG\n',build='GRCh38')
    def test_mixed_builds_rejected(self):
        with self.assertRaisesRegex(ValueError,'conflitantes'):
            self.parse(HEADER+'# GRCh38\nrs1 1 10 AG\n')
    def test_unknown_strand_rejected(self):
        with self.assertRaisesRegex(ValueError,'Fita desconhecida'):
            self.parse('# reference=GRCh37\nrs1 1 10 AG\n')
    def test_negative_strand_rejected(self):
        with self.assertRaisesRegex(ValueError,'fita positiva'):
            self.parse(HEADER+'rs1 1 10 AG\n',strand='-')
    def test_nocalls_indels_and_duplicates(self):
        g=self.parse(HEADER+'rs1 1 10 AG\nrs1 1 10 AA\nrs2 1 11 --\nrs3 1 12 DI\n')
        self.assertEqual(g.summary()['called'],0)
        self.assertEqual(g.stats['conflicts'],1)
        self.assertEqual(search_calls(g,'rs1')[0]['reason'],'Duplicatas conflitantes')
    def test_unphased_allele_order_is_not_a_conflict(self):
        g=self.parse(HEADER+'rs1 1 10 AG\nrs1 1 10 GA\n')
        self.assertEqual(g.summary()['called'],1)
        self.assertEqual(g.stats['conflicts'],0)
    def test_malformed_cancels_import(self):
        with self.assertRaisesRegex(ValueError,'Importação cancelada'):
            self.parse(HEADER+'rs1 1 10 AG\nbad line\n')
    def test_coordinate_rsid_disagreement(self):
        g=self.parse(HEADER+'rs999 2 136608646 AA\n')
        m=next(m for m in panel_results(g) if m['rsid']=='rs4988235')
        self.assertEqual(m['status'],'unavailable')
    def test_rsid_at_wrong_coordinate_never_matches(self):
        g=self.parse(HEADER+'rs4988235 2 100 AA\n')
        m=next(m for m in panel_results(g) if m['rsid']=='rs4988235')
        self.assertEqual(m['status'],'unavailable')
    def test_vcf_multiallelic_phased_and_no_calls(self):
        g=self.parse(VCF+'1\t10\trs1\tA\tC,G\t99\tPASS\t.\tGT:DP:GQ\t1|2:30:99\n1\t11\t.\tA\tG\t99\tPASS\t.\tGT:DP:GQ\t0/.:30:99\n')
        self.assertEqual(g.calls[('1',10)].alleles,('C','G'))
        self.assertTrue(g.calls[('1',10)].phased)
        self.assertIsNone(g.calls[('1',11)].alleles)
    def test_vcf_low_quality_filter_and_haploid(self):
        lines=[('10','PASS','0/1:9:99'),('11','PASS','0/1:30:19'),('12','LowQual','0/1:30:99'),('13','PASS','1:30:99'),('14','PASS','3/1:30:99')]
        g=self.parse(VCF+''.join(f'1\t{pos}\t.\tA\tG\t99\t{flt}\t.\tGT:DP:GQ\t{gt}\n' for pos,flt,gt in lines))
        self.assertEqual(g.summary()['called'],0)
    def test_multisample_needs_selection(self):
        header=VCF.rstrip()+'\tOTHER\n'
        row='1\t10\t.\tA\tG\t99\tPASS\t.\tGT\t0/0\t0/1\n'
        with self.assertRaisesRegex(ValueError,'múltiplas amostras'):
            self.parse(header+row)
        self.assertEqual(self.parse(header+row,sample='OTHER').calls[('1',10)].alleles,('A','G'))
    def test_gzip(self):
        with gzip.open(self.path,'wt',encoding='utf-8') as stream: stream.write(HEADER+'rs1 1 10 AG\n')
        self.assertEqual(parse_file(self.path).summary()['called'],1)
    def test_search_interval_validation(self):
        g=self.parse(HEADER+'rs1 1 10 AG\n')
        self.assertEqual(len(search_calls(g,'1:10-20')),1)
        with self.assertRaises(ValueError): search_calls(g,'1:20-10')
        with self.assertRaises(ValueError): search_calls(g,'arbitrary')
    def test_pseudoautosomal_25_not_mitochondrial(self):
        g=self.parse(HEADER+'rs1 25 10 AG\nrs2 26 10 A\n')
        self.assertEqual(g.summary()['loci'],1)
        self.assertIn(('MT',10),g.calls)


class ClinicalTests(unittest.TestCase):
    def setUp(self):
        self.temp=workspace_temp()
        self.patch=patch.object(clinical,'DATA',pathlib.Path(self.temp.name))
        self.patch.start()
        self.path=pathlib.Path(self.temp.name)/'public.vcf'
        self.path.write_text('##fileformat=VCFv4.2\n##reference=GRCh37\n##fileDate=2026-10-19\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
            '1\t10\t100\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;CLNDN=Example_condition;GENEINFO=TEST:1;RS=1\n'
            '1\t11\t101\tC\tT\t.\t.\tCLNSIG=Conflicting_classifications_of_pathogenicity;CLNREVSTAT=criteria_provided,_conflicting_classifications;CLNDN=Other;GENEINFO=TEST:1\n'
            '1\t12\t102\tA\tAT\t.\t.\tCLNSIG=Likely_pathogenic;CLNREVSTAT=criteria_provided,_single_submitter;CLNDN=Indel;GENEINFO=TEST:1\n',encoding='utf-8')
        clinical.import_clinvar(self.path,'GRCh37')
    def tearDown(self):
        self.patch.stop();self.temp.cleanup()
    def genome(self): return Genome('synthetic','GRCh37','Genotipagem TXT/CSV','hash')
    def test_chip_exact_alleles_not_rsid_only(self):
        g=self.genome()
        g.add(Call('1',10,('rs1',),('A','G')))
        g.add(Call('1',11,(),('C','T')))
        result=clinical.analyze(g)
        self.assertEqual(result['summary']['pathogenic'],1)
        self.assertEqual(result['summary']['conflicting'],1)
        self.assertEqual(result['findings'][0]['stars'],3)
    def test_reference_homozygous_is_not_hit(self):
        g=self.genome();g.add(Call('1',10,('rs1',),('A','A')))
        self.assertEqual(clinical.analyze(g)['total_findings'],0)
    def test_incompatible_chip_alleles_rejected(self):
        g=self.genome();g.add(Call('1',10,('rs1',),('C','G')))
        self.assertEqual(clinical.analyze(g)['total_findings'],0)
    def test_nocall_is_not_hit(self):
        g=self.genome();g.add(Call('1',10,('rs1',),None,reason='sem chamada'))
        self.assertEqual(clinical.analyze(g)['total_findings'],0)
    def test_vcf_ref_mismatch(self):
        g=self.genome();g.format='VCF';g.add(Call('1',10,('rs1',),('C','G'),'C',('G',)))
        self.assertEqual(clinical.analyze(g)['total_findings'],0)
    def test_vcf_exact_indel(self):
        g=self.genome();g.format='VCF';g.add(Call('1',12,(),('A','AT'),'A',('AT',)))
        self.assertEqual(clinical.analyze(g)['summary']['pathogenic'],1)
    def test_other_assembly_never_matches(self):
        g=self.genome();g.build='GRCh38';g.add(Call('1',10,('rs1',),('A','G')))
        self.assertFalse(clinical.analyze(g)['available'])
    def test_clinvar_wrong_assembly(self):
        with self.assertRaisesRegex(ValueError,'assemblies diferentes'):
            clinical.import_clinvar(self.path,'GRCh38')
    def test_failed_update_preserves_database(self):
        self.path.write_text('invalid',encoding='utf-8')
        with self.assertRaises(ValueError): clinical.import_clinvar(self.path,'GRCh37')
        self.assertEqual(clinical.info('GRCh37')['records'],'3')
    def test_classification_conflict_not_pathogenic(self):
        self.assertEqual(clinical.severity({'classification':'Pathogenic','conflict':'Benign(1)'}),'conflicting')


class AncestryTests(unittest.TestCase):
    def test_synthetic_known_mixture_and_quality_requirements(self):
        rng=random.Random(41)
        g=Genome('synthetic','GRCh37','chip','hash')
        panel={'assembly':'GRCh37','strand':'+','source':'SYNTHETIC TEST ONLY','ld_pruned':True,'populations':['A','B'],'markers':[]}
        for i in range(1200):
            ch=str(i%10+1);pos=10+i*1000
            f1,f2=(.9,.1) if i%3 else (.2,.8)
            probability=.7*f1+.3*f2
            alleles=tuple('G' if rng.random()<probability else 'A' for _ in range(2))
            g.add(Call(ch,pos,(),alleles))
            panel['markers'].append({'chrom':ch,'pos':pos,'other_allele':'A','effect_allele':'G','frequencies':{'A':f1,'B':f2}})
        result=estimate(g,panel)
        self.assertAlmostEqual(result['proportions'][0]['percent'],70,delta=6)
        self.assertTrue(result['converged'])
        panel['assembly']='GRCh38'
        with self.assertRaisesRegex(ValueError,'mesmo assembly'): estimate(g,panel)
    def test_low_coverage_no_fabricated_proportions(self):
        g=Genome('synthetic','GRCh37','chip','hash')
        p={'assembly':'GRCh37','strand':'+','source':'TEST','ld_pruned':True,'populations':['A','B'],'markers':[]}
        with self.assertRaisesRegex(ValueError,'Cobertura insuficiente'): estimate(g,p)
    def test_linearly_dependent_populations_rejected(self):
        g=Genome('synthetic','GRCh37','chip','hash')
        panel={'assembly':'GRCh37','strand':'+','source':'TEST ONLY','ld_pruned':True,'populations':['A','B','C'],'markers':[]}
        for i in range(600):
            ch=str(i%5+1);pos=i+1000
            g.add(Call(ch,pos,(),('A','G')))
            f=.9 if i%2 else .8
            panel['markers'].append({'chrom':ch,'pos':pos,'other_allele':'A','effect_allele':'G','frequencies':{'A':f,'B':1-f,'C':.5}})
        with self.assertRaisesRegex(ValueError,'mistura única'): estimate(g,panel)


if __name__=='__main__': unittest.main()
