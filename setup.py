from setuptools import setup, find_packages

setup(name='MCI',
        version='0.1',
        packages=['mci'],
#        packages=find_packages(),
        url='https://github.com/gicsaw/MCI',
        license='MIT LICENSE',
        author='Seung Hwan Hong',
        author_email='gicsaw0@gmail.com',
        description='',
        scripts=['bin/draw_binding_residue.py',
                 'bin/draw_binding_residue_list.py',
                 'bin/find_binding_residue.py',
                 'bin/find_pf_receptor.py',
                 'bin/find_pocket.py',
                 'bin/mci_draw.py',
                 'bin/mci_draw_receptor.py',
                 'bin/mci_score.py',
                 'bin/mci_score_vhts.py'
                ]
)


