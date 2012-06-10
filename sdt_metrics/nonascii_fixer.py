# run with python 3k

with open('_sdt_metrics.bak', encoding='latin-1') as fid:
    with open('_sdt_metrics.py','wb') as fod:
        s = fid.read()
        
        fod.write(s.encode("ascii", "ignore"))
        
