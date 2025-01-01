'''
# Description
Functions to move files around.

# Index
- `get()`
- `get_list()`
- `copy()`
- `move()`
- `remove()`
- `rename()`
- `rename_on_subfolders()`
- `copy_to_subfolders()`
- `from_template()`

---
'''


import os
import shutil


def get(
        filename:str,
        filters=None,
        return_anyway:bool=False
        ) -> str:
    '''
    Check if the given `filename` exists in the currrent working directory
    or in the full path, and returns its full path as a string.\n
    Raises an error if the file is not found, unless `return_anyway=True`,
    in which case it returns None. This can be used to personalize errors.\n
    If the provided string is a directory, it checks the files inside it.
    if there is only one file inside, it returns said file;
    if there are more files, it tries to filter them with the `filters` keyword(s) to return a single file.
    If this fails, try using more strict filers to return a single file.
    '''
    if os.path.isfile(filename):
        return os.path.abspath(filename)
    elif os.path.isdir(filename):
        files = get_list(filename, filters, abspath=True)
    elif return_anyway:
        return None
    else:
        raise FileNotFoundError('Nothing found at ' + str(filename))
    # Return a single file
    if len(files) == 1:
        return files[0]
    elif return_anyway:
        return None
    elif len(files) == 0:
        raise FileNotFoundError('The following directory is empty (maybe due to the filters):\n' + filename)
    else:
        raise FileExistsError(f'More than one file found, please apply a more strict filter. Found:\n{files}')


def get_list(
        folder:str,
        filters=None,
        abspath:bool=True
    ) -> list:
    '''
    Takes a `folder`, filters the content with the `filters` keyword(s) if provided, and returns a list with the matches.
    The full paths are returned by default; to get only the base names, set `abspath=False`.
    '''
    if os.path.isfile(folder):
        folder = os.path.dirname(folder)
    if not os.path.isdir(folder):
        raise FileNotFoundError('Nothing found at ' + folder)
    folder = os.path.abspath(folder)
    files = os.listdir(folder)
    # Apply filters or not
    if filters is not None:
        target_files = []
        if not isinstance(filters, list):
            filters = [str(filters)]
        for filter_i in filters:
            filter_i = os.path.basename(filter_i)
            for f in files:
                if filter_i in f:
                    target_files.append(f)
        files = target_files
    if abspath:
        filepaths = []
        for f in files:
            filepaths.append(os.path.join(folder, f))
        files = filepaths
    return files


def copy(
        old:str,
        new:str
    ) -> None:
    '''
    Copies the content of `old` file to `new` file with shutil,
    after making sure that the file exists with `thotpy.file.get()`.
    '''
    original_file_path = get(old)
    file = shutil.copy(original_file_path, new)
    return None


def move(
        old:str,
        new:str
    ) -> None:
    '''
    Moves `old` file to `new` file.
    '''
    original_file_path = get(old)
    file = shutil.move(original_file_path, new)
    return None


def remove(filename:str) -> None:
    '''
    Removes the given file with `filename`.
    '''
    file_path = get(filename)
    if os.path.isdir(file_path):
        shutil.rmtree(file_path)
    elif os.path.isfile(file_path):
        os.remove(file_path)
    else:
        raise FileNotFoundError(f"No such file or directory: '{file_path}'")
    return None


def rename(
        old:str,
        new:str,
        folder=None
    ) -> None:
    '''
    Batch renames files in the given folder, replacing `old` string by `new` string.
    If no `folder` is provided, the current working directory is used.
    '''
    if folder is None:
        files = os.listdir('.')
    elif os.path.isdir(folder):
        files = os.listdir(folder)
    elif os.path.isdir(os.path.join(os.getcwd(), folder)):
        files = os.listdir(os.path.join(os.getcwd(), folder))
    else:
        raise FileNotFoundError('Missing folder at ' + folder + ' or in the CWD ' + os.getcwd())
    for f in files:
        if old in f:
            os.rename(f, f.replace(old, new))
    return None


def rename_on_subfolders(
        old:str,
        new:str,
        folder=None
    ) -> None:
    '''
    Renames the files inside the subfolders in the given `folder`,
    from an `old` string to the `new` string.
    If no `folder` is provided, the current working directory is used.
    '''
    if folder is None:
        things = os.listdir('.')
    elif os.path.isdir(folder):
        things = os.listdir(folder)
    elif os.path.isdir(os.path.join(os.getcwd(), folder)):
        things = os.listdir(os.path.join(os.getcwd(), folder))
    else:
        raise FileNotFoundError('Missing folder at ' + folder + ' or in the CWD ' + os.getcwd())
    for d in things:
        if os.path.isdir(d):
            for f in os.listdir(d):
                if old in f:
                    old_file = os.path.join(d, f)
                    new_file = os.path.join(d, f.replace(old, new))
                    os.rename(old_file, new_file)
    return None


def copy_to_subfolders(
        folder=None,
        extension:str=None,
        strings_to_delete:list=[]
    ) -> None:
    '''
    Copies the files from the `folder` with the given `extension` to individual subfolders.
    The subfolders are named as the original files,
    removing the strings from the `strings_to_delete` list.
    If no `folder` is provided, it runs in the current working directory.
    '''
    if folder is None:
        folder = os.getcwd()
    old_files = get_list(folder, extension)
    if old_files is None:
        raise ValueError('No ' + extension + ' files found in path!')
    for old_file in old_files:
        new_file = old_file
        for string in strings_to_delete:
            new_file = new_file.replace(string, '')
        path = new_file.replace(extension, '')
        os.makedirs(path)
        new_file_path = path + '/' + new_file
        copy(old_file, new_file_path)
    return None


def from_template(
        old:str,
        new:str,
        comment:str=None,
        fixing_dict:dict=None
    ) -> None:
    '''
    Similar to `copy_file()`, but optionally adds a `comment` at the beginning of the new file.
    Also, it optionally corrects the output file with a `fixing_dict` dictionary.
    `old` is the template and `new` is the final file.
    '''
    copy(old, new)
    if comment:
        with open(new, 'r+') as f:
            content = f.read()
            f.seek(0)
            f.write(comment + '\n' + content)
    if fixing_dict:
        with open(new, 'r+') as f:
            content = f.read()
            for key, value in fixing_dict.items():
                content = content.replace(key, value)
            f.seek(0)
            f.write(content)
            f.truncate()
    return None

