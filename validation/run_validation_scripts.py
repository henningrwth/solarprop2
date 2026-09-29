#! /usr/bin/env python3

"""
Run all tests of validation models consecutively.
"""

import argparse
import os
import shlex
import subprocess

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument('-r', '--reuse', action='store_true', help='Re-use Green matrix files from previous runs.')
args = argparser.parse_args()

logfilename = 'run_validation_scripts.out'

test_scripts = ['test_yamada_fig4.py',
                'test_jokipii_kopriva_fig2.py',
                'test_potgieter_moraal_fig6.py',
                'test_burger_potgieter_fig5.py',
                'test_burger_2012_fig2.py',
                'test_strauss_2012_fig15.py',
                'test_strauss_2012_fig11.py',
                'test_strauss_2012_fig5.py',
                'test_strauss_2012_variations.py',
                'test_helmod.py',
                ]

with open(logfilename, 'w', encoding='utf-8') as logfile:

    for test_script in test_scripts:

        header = f'======== {test_script} ========'
        print(header)
        logfile.write(header)
        logfile.write('\n')

        # get script with full path and add command-line arguments
        script_exe = os.path.join(os.getcwd(), test_script)
        cmd = f'{script_exe} -n 10000 -b'
        if args.reuse:
            cmd += ' -r'

        res = subprocess.run(shlex.split(cmd), check=False, capture_output=True, text=True, cwd=os.getcwd())
        logfile.write(res.stdout)
        if res.stderr:
            logfile.write('======== ERRORS: ========\n')
            logfile.write(res.stderr)
            print('======== ERRORS: ========')
            print(res.stderr)
        if res.returncode:
            logfile.write(f'Exited with return code {res.returncode}')
        res.check_returncode()
        logfile.write('\n')


print(f'Logfile: {logfilename}')
