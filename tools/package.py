"""Empacota apenas código e exemplos sintéticos. Nunca inclui DNA pessoal ou relatórios."""
import pathlib
import zipfile

ROOT=pathlib.Path(__file__).resolve().parents[1]
FILES=['server.py','genetics.py','clinical.py','ancestry.py','ancestry_projection.py','catalog.py','requirements.txt','README.md','THIRD_PARTY_NOTICES.md','iniciar.cmd','data/markers.json']
FILES.extend(str(p.relative_to(ROOT)).replace('\\','/') for p in (ROOT/'docs').glob('*.md'))
if (ROOT/'docs/images/interface.jpg').exists():
    FILES.append('docs/images/interface.jpg')
if (ROOT/'data/ancestry/nhgri_k21_grch37.json').exists():
    FILES.append('data/ancestry/nhgri_k21_grch37.json')
for folder,patterns in [('web',['*.html','*.css','*.js']),('tools',['*.py']),('tests',['*.py','*.cjs']),('examples',['*.txt','*.vcf'])]:
    for pattern in patterns:
        FILES.extend(str(p.relative_to(ROOT)).replace('\\','/') for p in (ROOT/folder).glob(pattern))

if __name__=='__main__':
    output=ROOT/'helice-portatil.zip'
    with zipfile.ZipFile(output,'w',compression=zipfile.ZIP_DEFLATED) as archive:
        for name in sorted(set(FILES)):
            archive.write(ROOT/name,'helice/'+name)
    print(f'Pacote criado com {len(set(FILES))} arquivos de código, documentação e exemplos sintéticos. Bases ClinVar grandes não incluídas.')
