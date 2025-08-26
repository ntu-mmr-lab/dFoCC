# ## CCP4 & Phenix wrappers
# 
# As CCP4 and Phenix commands can not be natively called in Python, we map them into Python functions.

import subprocess

def ccp4(program_name: str):
    def program(script='', **kwargs):
        run_ccp4 = [program_name]
        for key, val in kwargs.items():
            run_ccp4.extend([key, val])
        # `ccp4_program keyin1 ./key1 keyin2 ./key2 <<< $script`
        complete_process = subprocess.run(run_ccp4, input=script.encode(), capture_output=True)
        return complete_process.stdout.decode()
    return program

mapmask = ccp4('mapmask')
sfall = ccp4('sfall')
fft = ccp4('fft')
sftools = ccp4('sftools')
mapdump = ccp4('mapdump')
overlapmap = ccp4('overlapmap')

def phenix(program_name: str, capture_output=True):
    def program(*args, **kwargs):
        run_phenix = [f'phenix.{program_name}']
        run_phenix.extend(args)
        run_phenix.extend([
            f'{key}={val}' for key, val in kwargs.items()
        ])
        # `map_box key_in1=./key1 key_in2=./key2`
        complete_process = subprocess.run(run_phenix, capture_output=capture_output)
        if not capture_output: return ''
        return complete_process.stdout.decode() + complete_process.stderr.decode()
    return program

map_box = phenix('map_box')
fobs_minus_fobs_map = phenix('fobs_minus_fobs_map')
pdb_atom_selection = phenix('pdb_atom_selection')
phenix_refine = phenix('refine', capture_output=False)
