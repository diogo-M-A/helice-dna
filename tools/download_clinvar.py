import pathlib
import sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]))
import clinical

if __name__ == '__main__':
    build=sys.argv[1] if len(sys.argv)>1 else 'GRCh37'
    last=[0]
    def progress(message,count):
        import time
        if time.monotonic()-last[0]>5:
            print(message,flush=True)
            last[0]=time.monotonic()
    result=clinical.download(build,progress)
    print(f"ClinVar {build}: {result['records']} variantes; edição {result['file_date']}")
