'''
# Description
Functions to handle bash calls and related operations on Linux systems.

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
from .core import *


def bash(
        command:str,
        cwd=None,
        verbose:bool=True,
        return_anyway:bool=False
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


def git(
        path=None,
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


def here(folder=None) -> str:
    '''
    Runs the rest of the script inside the specified `folder`.
    If none is provided, it runs from the same directory where the current script lies.
    This is really useful to run scripts from the VSCode terminal, etc.
    Returns the path of the used `folder`, or the path of the script if folder is not provided.\n
    Note that this changes not only the working directory of your script,
    but also of other scripts that import and run your script.
    '''
    if folder:
        caller = os.path.abspath(folder)
    else:
        caller = os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0])))
    os.chdir(caller)
    return caller

