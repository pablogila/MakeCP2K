'''
# Description
Functions to simplify calling bash scripts and related.

# Index
- `bash()`
- `git()`
- `here()`

---
'''


import subprocess
import datetime
import sys
import os
from .common import *


def bash(command,
         cwd=None,
         verbose=True,
         return_anyway=False
         ):
    '''
    Run a bash shell `command`, inside an optional `cwd` directory.
    If empty, the current working directory will be used.
    Prints the running command and outputs by default, override this with `verbose=False`.
    Returns the result of the command used, except for when
    errors are raised automatically; set `return_anyway=True` to override this.
    '''
    if verbose:
        print(f'$ {command}')
    result = subprocess.run(command, cwd=cwd, shell=True, text=True, capture_output=True)
    if verbose and result.returncode == 0 and result.stdout:
        print(result.stdout)
    elif result.returncode != 0:
        error_message = (
            f"thoth.call.bash: Command failed with exit code {result.returncode}.\n"
            f"{result.stderr.strip()}"
        )
        if not return_anyway:
            raise RuntimeError(error_message)
    return result


def git(path=None,
        verbose=True,
        message=None,
        tag=None
        ) -> None:
    '''Automatically update a Git repository.'''
    if path:
        os.chdir(path)
    bash("git fetch", path, verbose)
    bash("git add .", path, verbose)
    if not message:
        date = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M")
        message = f'Automatic push on {date} with Thoth {version}'
    bash(f'git commit -m "{message}"', path, verbose)
    if tag:
        bash(f'git tag -a {tag} HEAD -m {message}', path, verbose)
        bash("git push origin main --tags", path, verbose)
    else:
        bash("git push", path, verbose)
    print("Git updated!")
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

