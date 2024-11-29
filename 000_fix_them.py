
fnames = [
'sps_match_injectiontrajectory.py',
'YASP_injectionOrbitReference.txt',
'YASP_injectionTrajectory.txt',
'sps_cpymad.py',
]

for fname in fnames:
    with open(fname, 'r') as fid:
        content = fid.read()

    content = content.replace('\n\n', '\n')

    with open(fname, 'w') as fid:
        fid.write(content)
