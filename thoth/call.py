'''
# Description
Functions to simplify calling bash scripts and related.

# Index
- `shell()`
- `git()`
- `here()`

---
'''


import subprocess
import datetime
import sys
import os
import io
from .common import *


def shell(command, cwd=None, verbose=True):
    '''
    Run a shell `command`, inside an optional `cwd` directory.
    If empty, the current working directory will be used.
    Returns the result of the command used.
    '''
    if verbose:
        print(f'$ {command}')
    result = subprocess.run(
        command,
        cwd=cwd,
        shell=True, 
        text=True, 
        capture_output=True
    )
    if verbose and result.returncode == 0 and result.stdout:
        print(result.stdout)
    elif result.returncode != 0:
        error_message = (
            f"thoth.call.shell: Command failed with exit code {result.returncode}.\n"
            f"{result.stderr.strip()}"
        )
        raise RuntimeError(error_message)


def git(path=None, verbose=True) -> None:
    '''Update'''
    if path:
        os.chdir(path)
    date = datetime.datetime.now().strftime("%Y.%m.%dT%H:%M")
    shell("git fetch", path, verbose)
    shell("git add .", path, verbose)
    shell(f'git commit -m "Automatic push on {date} with Thoth {version}"', path, verbose)
    shell("git push", path, verbose)
    print("Git updated!", path, verbose)
    return None


def here(run_here=True):
    '''
    Returns the directory where the current script lies.
    By default, it also runs the rest of the script from said directory;
    This is really useful to run scripts from the VSCode terminal, etc.
    You might want to override this behaviour if you just want to know the path of the current script;
    to do so, set `run_here=False`.
    '''
    caller = os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0])))
    if run_here:
        os.chdir(caller)
    return caller

