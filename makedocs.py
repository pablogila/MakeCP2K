'''
This script is used to update Thoth documentation automatically.
Requires pdoc, install it with `pip install pdoc`.
It also requires Thoth itself; installation instructions can be found in the README.md file.
Run this script as `python3 pdoc.py`.
'''

import thoth as th

readme = './README.md'
temp_readme = './_README_temp.md'
version_path = './thoth/autoload.py'

fix_dict ={
    '[autoload](https://pablogila.github.io/Thoth/thoth/autoload.html)' : '`thoth.autoload`',
    '[alias](https://pablogila.github.io/Thoth/thoth/alias.html)'       : '`thoth.alias`',
    '[call](https://pablogila.github.io/Thoth/thoth/call.html)'         : '`thoth.call`',
    '[file](https://pablogila.github.io/Thoth/thoth/file.html)'         : '`thoth.file`',
    '[find](https://pablogila.github.io/Thoth/thoth/find.html)'         : '`thoth.find`',
    '[extract](https://pablogila.github.io/Thoth/thoth/extract.html)'   : '`thoth.extract`',
    '[text](https://pablogila.github.io/Thoth/thoth/text.html)'         : '`thoth.text`',
    '[qe](https://pablogila.github.io/Thoth/thoth/call.html)'           : '`thoth.qe`',
    '[phonopy](https://pablogila.github.io/Thoth/thoth/phonopy.html)'   : '`thoth.phonopy`',
} 

version = th.find.lines(r"version\s*=", version_path, -1, 0, False, True)[0]
version = th.extract.string(version, 'version', None, True)

print(f'Updating README to {version}...')
th.text.replace_line(f'# Thoth {version}', '# Thoth v', readme, 1)

print('Updating docs with Pdoc...')
th.file.from_template(readme, temp_readme, None, fix_dict)
th.call.bash(f"pdoc ./thoth/ -o ./docs --mermaid --math --footer-text='Thoth {version} documentation'")
th.file.remove(temp_readme)
print('Done!')

