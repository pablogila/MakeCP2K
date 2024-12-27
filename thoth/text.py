'''
# Description
Functions to manipulate text files.

# Index
- `insert_at()`
- `insert_under()`
- `replace()`
- `replace_line()`
- `replace_between()`
- `delete_under()`
- `delete_between()`
- `correct_with_dict()`

---
'''


import mmap
from .file import *
from .find import *


def insert_at(text:str,
              file,
              position:int
              ) -> None:
    '''
    Inserts a `text` in the line with `position` index of a given `file`.
    If `position` is negative, starts from the end of the file.
    '''
    file_path = get(file)
    with open(file_path, 'r+') as f:
        lines = f.read().splitlines()
        if position < 0:
            position = len(lines) + position + 1
        if position < 0 or position > len(lines):
            raise IndexError("Position out of range")
        lines.insert(position, text)
        f.seek(0)
        f.write('\n'.join(lines))
        f.truncate()
    return None


def insert_under(text:str,
                 keyword:str,
                 file,
                 number_of_inserts:int=0,
                 skip_lines:int=0,
                 regex:bool=False
                 ) -> None:
    '''
    Inserts the given `text` string under the line(s) containing
    the `keyword` in the given `file`.
    The keyword can be at any position within the line.
    By default all matches are inserted with `number_of_inserts=0`,
    but it can insert only a specific number of matches
    with positive numbers (1, 2...), or starting from the bottom with negative numbers.
    The text can be introduced after a specific number of lines after the match,
    changing the value `skip_lines`. Negative integers introduce the text in the previous lines.
    Regular expressions can be used by setting `regex=True`. 
    '''
    file_path = get(file)
    if regex:
        positions = pos_regex(keyword, file, number_of_inserts)
    else:
        positions = pos(keyword, file, number_of_inserts)
    positions.reverse()  # Must start replacing from the end, otherwise the atual positions may change!
    # Open the file in read-write mode
    with open(file_path, 'r+b') as f:
        with mmap.mmap(f.fileno(), length=0, access=mmap.ACCESS_WRITE) as mm:
            # Get the places to insert the text
            for position in positions:
                start, end = line_pos(position, mm, skip_lines)
                inserted_text = '\n' + text # Ensure we end in a different line
                if end == 0: # If on the first line
                    inserted_text = text + '\n'
                remaining_lines = mm[end:]
                new_line = inserted_text.encode()
                updated_content = new_line + remaining_lines
                mm.resize(len(mm) + len(new_line))
                mm[end:] = updated_content
    return None


def replace(text:str,
            keyword:str,
            file:str,
            number_of_replacements:int=0,
            regex:bool=False
            ) -> None:
    '''
    Replaces the `keyword` string with the `text` string in the given `file`.\n
    The value `number_of_replacements` specifies the number of replacements to perform:
    1 to replace only the first keyword found, 2, 3...
    Use negative values to replace from the end of the file,
    eg. to replace the last found key, use `number_of_replacements=-1`.
    To replace all values, set `number_of_replacements = 0`, which is the value by default.\n
    To search with regular expressions, set `regex=True`.
    ```
    line... keyword ...line -> line... text ...line
    ```
    '''
    file_path = get(file)
    if regex:
        positions = pos_regex(keyword, file, number_of_replacements)
    else:
        positions = pos(keyword, file, number_of_replacements)
    positions.reverse()  # Must start replacing from the end, otherwise the atual positions may change!
    with open(file_path, 'r+') as f:
        content = f.read()
        for start, end in positions:
            content = "".join([content[:start], text, content[end:]])
        f.seek(0)
        f.write(content)
        f.truncate()


def replace_line(text:str,
                 keyword:str,
                 file:str,
                 number_of_replacements:int=0,
                 skip_lines:int=0,
                 additional_lines:int=0,
                 regex:bool=False
                 ) -> None:
    '''
    Replaces the entire line containing the `keyword` string with the `text` string in the given `file`.
    Regular expressions can be used with `regex=True`.\n
    The value `number_of_replacements` specifies the number of lines to replace:
    1 to replace only the first line with the keyword, 2, 3...
    Use negative values to replace from the end of the file,
    e.g., to replace only the last line containing the keyword, use `number_of_replacements = -1`.
    To replace all lines, set `number_of_replacements = 0`, which is the value by default.\n
    The default line to replace is the matching line,
    but it can be any other specific line after or before the matching line;
    this is indicated with `skip_lines` as a positive or negative integer.\n
    More lines can be replaced with `additional_lines` (int).
    Note that the matched line plus the additional lines will be replaced, this is, additional lines +1.
    '''
    file_path = get(file)
    if regex:
        positions = pos_regex(keyword, file, number_of_replacements)
    else:
        positions = pos(keyword, file, number_of_replacements)
    positions.reverse()  # Must start replacing from the end, otherwise the atual positions may change!
    # Open the file in read-write mode
    with open(file_path, 'r+b') as f:
        with mmap.mmap(f.fileno(), length=0, access=mmap.ACCESS_WRITE) as mm:
            for position in positions:
                # Get the positions of the full line containing the match
                line_start, line_end = line_pos(position, mm, skip_lines)
                # Additional lines
                if additional_lines > 0:
                    for _ in range(abs(additional_lines)):
                        line_end = mm.find(b'\n', line_end + 1, len(mm)-1)
                        if line_end == -1:
                            line_end = len(mm) - 1
                            break
                elif additional_lines < 0:
                    for _ in range(abs(additional_lines)):
                        line_start = mm.rfind(b'\n', 0, line_start - 1) + 1
                        if line_start == -1:
                            line_start = 0
                            break
                # Replace the line
                old_line = mm[line_start:line_end]
                new_line = text.encode()
                # Directly modify the memory-mapped region
                if len(new_line) == len(old_line):
                    mm[line_start:line_end] = new_line
                else:  # Adjust content for differing line sizes
                    remaining_content = mm[line_end:]
                    updated_content = new_line + remaining_content
                    mm.resize(len(mm) + len(new_line) - len(old_line))
                    mm[line_start:] = updated_content
    return None
    

def replace_between(text:str,
                    key1:str,
                    key2:str,
                    file:str,
                    reverse_match:bool=False,
                    regex:bool=False
                    ) -> None:
    '''
    Replace lines with a given `text`, between the keywords `key1` and `key2` in a given `file`.
    Keywords can be at any position within the line.
    Regular expressions can be used by setting `regex=True`.\n
    Only the first matches of the keywords are used by default;
    you can use the last ones with `reverse_match=True`.
    ```
    lines...
    key1
    text
    key2
    lines...
    ```
    '''
    index = 1
    if reverse_match:
        index = -1
    delete_between(key1, key2, file, False, reverse_match, regex)
    insert_under(text, key1, file, index, 0, regex)
    return None


def delete_under(keyword:str,
                 file,
                 match_index:int=1,
                 skip_lines:int=0,
                 regex:bool=False
                 ) -> None:
    '''
    Deletes the content under the line containing the `keyword` in the given `file`.
    The keyword can be at any position within the line.
    Regular expressions can be used by setting `regex=True`.\n
    By default the first `match_index` is used; it can be any integer,
    including negative integers to select a match starting from the end of the file.\n
    The content can be deleted after a specific number of lines after the match,
    changing the value `skip_lines`. Negative integers start deleting the content from the previous lines.
    '''
    file_path = get(file)
    if regex:
        positions = pos_regex(keyword, file, match_index)
    else:
        positions = pos(keyword, file, match_index)
    if match_index > 0:  # We only want one match, and should be the last if match_index > 0
        positions.reverse()
    position = positions[0]
    # Open the file in read-write mode
    with open(file_path, 'r+b') as f:
        with mmap.mmap(f.fileno(), length=0, access=mmap.ACCESS_WRITE) as mm:
            # Get the places to insert the text
            start, end = line_pos(position, mm, skip_lines)
            mm.resize(len(mm) - len(mm[end:]))
            mm[end:] = b''
    return None


def delete_between(key1:str,
                   key2:str,
                   file,
                   delete_keys:bool=False,
                   reverse_match:bool=False,
                   regex:bool=False
                   ) -> None:
    '''
    Deletes the content between the line containing the `key1` and `key2` in the given `file`.
    Keywords can be at any position within the line.
    Regular expressions can be used by setting `regex=True`.\n
    Key lines are also deleted if `delete_keys=True`.\n
    Only the first matches of the keywords are used by default;
    you can use the last ones with `reverse_match=True`.
    '''
    file_path = get(file)
    index = 1
    if reverse_match:
        index = -1
    start, end = between_pos(key1, key2, file_path, delete_keys, index, regex)
    with open(file_path, 'r+b') as f:
        with mmap.mmap(f.fileno(), length=0, access=mmap.ACCESS_WRITE) as mm:
            remaining_lines = mm[end+1:]
            mm.resize(len(mm) - len(mm[start:end+1]))
            mm[start:] = remaining_lines
    return None


def correct_with_dict(file:str,
                      fixing_dict:dict
                      ) -> None:
    '''
    Corrects the given text `file` using the `fixing_dict` dictionary.
    '''
    file_path = get(file)
    with open(file_path, 'r+') as f:
        content = f.read()
        for key, value in fixing_dict.items():
            content = content.replace(key, value)
        f.seek(0)
        f.write(content)
    return None

