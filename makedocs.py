'''
This script is used to update ThothPy documentation automatically.
Requires pdoc, install it with `pip install pdoc`.
It also requires ThotPy itself; installation instructions can be found in the README.md file.
Run this script as `python3 makedocs.py`.
'''

import thotpy as th

readme = './README.md'
temp_readme = './_README_temp.md'
version_path = './thotpy/core.py'

fix_dict ={
    '[core](https://pablogila.github.io/ThotPy/thotpy/core.html)'         : '`thotpy.core`',
    '[call](https://pablogila.github.io/ThotPy/thotpy/call.html)'         : '`thotpy.call`',
    '[file](https://pablogila.github.io/ThotPy/thotpy/file.html)'         : '`thotpy.file`',
    '[find](https://pablogila.github.io/ThotPy/thotpy/find.html)'         : '`thotpy.find`',
    '[extract](https://pablogila.github.io/ThotPy/thotpy/extract.html)'   : '`thotpy.extract`',
    '[text](https://pablogila.github.io/ThotPy/thotpy/text.html)'         : '`thotpy.text`',
    '[qe](https://pablogila.github.io/ThotPy/thotpy/qe.html)'             : '`thotpy.qe`',
    '[phonopy](https://pablogila.github.io/ThotPy/thotpy/phonopy.html)'   : '`thotpy.phonopy`',
} 

version = th.find.lines(r"version\s*=", version_path, -1, 0, False, True)[0]
version = th.extract.string(version, 'version', None, True)

print(f'Updating README to {version}...')
th.text.replace_line(f'# ThotPy {version}', '# ThotPy v', readme, 1)

print('Updating docs with Pdoc...')
th.file.from_template(readme, temp_readme, None, fix_dict)
th.call.bash(f"pdoc ./thotpy/ -o ./docs --mermaid --math --footer-text='ThotPy {version} documentation'")
th.file.remove(temp_readme)
print('Done!')

