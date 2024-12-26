'''
# Description
Functions to read and manipulate text.

# Index
- `find_pos()`
- `find_pos_regex()`
- `find_pos_line()`
- `find()`
- `replace()`
- `replace_line()`
- `insert_under()`

The following functions work, but will be updated for faster performance with `find_pos` and `find_pos_regex`:

- `replace_under()`
- `delete_under()`
- `replace_between()`
- `delete_between()`
- `correct_with_dict()`

---
'''


from .file import *
import mmap
import re
from copy import deepcopy


def find_pos(keyword:str,
             file,
             number_of_matches:int=0
             ) -> list:
    '''
    Returns a list of the positions of a `keyword` in a given `file`.\n
    The value `number_of_matches` specifies the max number of matches to return.
    Defaults to 0 to return all possible matches. Set it to 1 to return only one match,
    or to negative integers to start searching from the end of the file upwards.\n
    This method is faster than `find_pos_regex()`, but does not search for regular expressions.
    '''
    file_path = get(file)
    positions = []
    with open(file_path, 'r+b') as f:
        mm = mmap.mmap(f.fileno(), length=0, access=mmap.ACCESS_READ)
    keyword_bytes = keyword.encode()
    if number_of_matches >= 0:
        start = 0
        while number_of_matches == 0 or len(positions) < number_of_matches:
            pos = mm.find(keyword_bytes, start)
            if pos == -1:
                break
            end = pos + len(keyword_bytes)
            positions.append((pos, end))
            start = end
    else:
        start = len(mm)
        while len(positions) < abs(number_of_matches):
            pos = mm.rfind(keyword_bytes, 0, start)
            if pos == -1:
                break
            end = pos + len(keyword_bytes)
            positions.append((pos, end))
            start = pos
        positions.reverse()
    return positions


def find_pos_regex(keyword:str,
                   file,
                   number_of_matches:int=0
                   ) -> list:
    '''
    Returns a list of the positions of a `keyword` in a given `file`.\n
    The value `number_of_matches` specifies the max number of matches to return.
    Defaults to 0 to return all possible matches. Set it to 1 to return only one match,
    or to negative integers to start searching from the end of the file upwards.\n
    This method is slower than `find_pos()`, but it can search for regular expressions.
    '''
    file_path = get(file)
    positions = []
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    if number_of_matches > 0:
        start = 0
        while len(positions) < number_of_matches:
            match = re.search(keyword, content[start:])
            if not match:
                break
            match_start = start + match.start()
            match_end = start + match.end()
            positions.append((match_start, match_end))
            start = match_end
    else:
        all_matches = list(re.finditer(keyword, content))
        if number_of_matches == 0:
            positions = [(match.start(), match.end()) for match in all_matches]
        else:
            positions = [(match.start(), match.end()) for match in all_matches[-abs(number_of_matches):]]
    return positions


def find_pos_line(mm,
                  position:tuple,
                  skip_lines:int=0
                  ) -> tuple:
    '''
    Returns the position of the full line containing the `position` tuple,
    in the given `mm` **memory mapped object**.
    A specific line below can be returned with `skip_lines` being a natural int,
    or previous lines with negative values.
    '''
    start, end = position
    if skip_lines == 0:
        start = mm.rfind(b'\n', 0, start) + 1
        end = mm.find(b'\n', end, len(mm))
    elif skip_lines > 0:
        for i in range(0, abs(skip_lines)):
            start = mm.find(b'\n', end, len(mm)) + 1
            if start == -1:
                start = len(mm)
                end = len(mm)
                break
            end = mm.find(b'\n', start, len(mm))
            if end == -1:
                start = len(mm)
                end = len(mm)
                break
    else:  # previous lines
        for i in range(0, abs(skip_lines)):
            end = mm.rfind(b'\n', 0, start)
            if end == -1:
                start = 0
                end = 0
                break
            start = mm.rfind(b'\n', 0, end) + 1
            if start == -1:
                start = 0
                end = 0
                break
    return start, end


def find(keyword:str,
         file:str,
         number_of_matches:int=0,
         additional_lines:int=0,
         split_additional_lines: bool=False,
         regex:bool=False
         ) -> list:
    '''
    Finds the line(s) containing the `keyword` string in the given `file`,
    returning a list with the matches.\n
    The value `number_of_matches` specifies the max number of matches to be returned.
    Defaults to 0 to return all possible matches. Set it to 1 to return only one match,
    or to negative integers to start the search from the end of the file upwards.\n
    The value `additional_lines` specifies the number of additional lines
    below the target line that are also returned;
    2 to return the found line plus two additional lines below, etc.
    Negative values return the specified number of lines before the target line.
    The original ordering from the file is preserved.
    Defaults to `additional_lines=0`, only returning the target line.
    By default, the additional lines are returned in the same list item as the match separated by a `\\n`,
    unless `split_additional_lines=True`, in which case they are added as additional items in the list.\n
    To use regular expressions in the search, set `regex=True`.
    By default regex search is deactivated, using the faster mmap.find and rfind methods instead.
    '''
    file_path = get(file)
    matches = []
    if regex:
        positions = find_pos_regex(keyword, file, number_of_matches)
    else:
        positions = find_pos(keyword, file, number_of_matches)
    with open(file_path, 'r+b') as f:
        mm = mmap.mmap(f.fileno(), length=0, access=mmap.ACCESS_READ)
    for start, end in positions:
        # Get the positions of the full line containing the match
        line_start = mm.rfind(b'\n', 0, start) + 1
        line_end = mm.find(b'\n', end, len(mm)-1)
        # Default values for the start and end of the line
        if line_start == -1: line_start = 0
        if line_end == -1: line_end = len(mm) - 1
        # Adjust the line_end to add additional lines after the match
        match_start = line_start
        match_end = line_end
        if additional_lines > 0:
            for _ in range(abs(additional_lines)):
                match_end = mm.find(b'\n', match_end + 1, len(mm)-1)
                if match_end == -1:
                    break
        elif additional_lines < 0:
            for _ in range(abs(additional_lines)):
                match_start = mm.rfind(b'\n', 0, match_start - 1) + 1
                if match_start == -1:
                    break
        # Save the matched lines
        matches.append(mm[match_start:match_end].decode())
    if split_additional_lines:
        splitted_matches = []
        for string in matches:
            splitted_match = string.splitlines()
            splitted_matches.extend(splitted_match)
        matches = splitted_matches
    return matches


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
        positions = find_pos_regex(keyword, file, number_of_replacements)
    else:
        positions = find_pos(keyword, file, number_of_replacements)
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
                 regex:bool=False
                 ) -> None:
    '''
    Replaces the entire line containing the `keyword` string with the `text` string in the given `file`.
    The value `number_of_replacements` specifies the number of lines to replace:
    1 to replace only the first line with the keyword, 2, 3...
    Use negative values to replace from the end of the file,
    e.g., to replace only the last line containing the keyword, use `number_of_replacements = -1`.
    The text can be replaced after a specific number of lines after the match,
    changing the value `skip_lines`. Negative integers replace the previous lines.
    To replace all lines, set `number_of_replacements = 0`, which is the value by default.
    ```
    line... keyword ...line -> text
    ```
    '''
    file_path = get(file)
    if regex:
        positions = find_pos_regex(keyword, file, number_of_replacements)
    else:
        positions = find_pos(keyword, file, number_of_replacements)
    positions.reverse()  # Must start replacing from the end, otherwise the atual positions may change!
    # Open the file in read-write mode
    with open(file_path, 'r+b') as f:
        with mmap.mmap(f.fileno(), length=0, access=mmap.ACCESS_WRITE) as mm:
            for position in positions:
                # Get the positions of the full line containing the match
                line_start, line_end = find_pos_line(mm, position, skip_lines)
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
    ```
    '''
    file_path = get(file)
    if regex:
        positions = find_pos_regex(keyword, file, number_of_inserts)
    else:
        positions = find_pos(keyword, file, number_of_inserts)
    positions.reverse()  # Must start replacing from the end, otherwise the atual positions may change!
    # Open the file in read-write mode
    with open(file_path, 'r+b') as f:
        with mmap.mmap(f.fileno(), length=0, access=mmap.ACCESS_WRITE) as mm:
            # Get the places to insert the text
            for position in positions:
                start, end = find_pos_line(mm, position, skip_lines)
                inserted_text = '\n' + text # Ensure we end in a different line
                if end == 0: # If on the first line
                    inserted_text = text + '\n'
                remaining_lines = mm[end:]
                new_line = inserted_text.encode()
                updated_content = new_line + remaining_lines
                mm.resize(len(mm) + len(new_line))
                mm[end:] = updated_content


def replace_under(text:str,
                  keyword:str,
                  file,
                  replace_last:bool=False,
                  skip_lines:int=0,
                  regex:bool=False
                  ) -> None:
    '''
    Replaces the content below the line containing the `keyword`
    by the given `text` string in the given `file`.
    The keyword can be at any position within the line.
    By default the replacement will take place below the first match,
    but it can be replaced from the last match with `replace_last=True`.
    The text can be replaced after a specific number of lines after the match,
    changing the value `skip_lines`. Negative integers replace the text in the previous lines.
    Regular expressions can be used by setting `regex=True`. 
    ```
    '''
    file_path = get(file)
    replace_index = 1
    if replace_last:
        replace_index = -1
    if regex:
        positions = find_pos_regex(keyword, file, replace_index)
    else:
        positions = find_pos(keyword, file, replace_index)
    positions.reverse()  # Must start replacing from the end, otherwise the atual positions may change!
    # Open the file in read-write mode
    with open(file_path, 'r+b') as f:
        with mmap.mmap(f.fileno(), length=0, access=mmap.ACCESS_WRITE) as mm:
            # Get the places to insert the text
            position = positions[0]
            start, end = find_pos_line(mm, position, skip_lines)
            inserted_text = '\n' + text # Ensure we end in a different line
            if end == 0: # If on the first line
                inserted_text = text + '\n'
            remaining_lines = mm[end:]
            new_line = inserted_text.encode()
            mm.resize(len(mm) + len(new_line) - len(remaining_lines))
            mm[end:] = new_line


def delete_under(keyword:str, file:str) -> None:
    '''
    Deletes the lines under the first occurrence of the `keyword` in the given `file`.
    > TODO: IN THE FUTURE SHOULD BE POSITION-AGNOSTIC
    ```
    lines...
    keyword
    (end of file)
    ```
    '''
    file_path = get(file)
    with open(file_path, 'r') as f:
        lines = f.readlines()
    keep = []
    for line in lines:
        if keyword in line:
            break
        else:
            keep.append(line)
    with open(file, 'w') as f:
        f.writelines(keep)
    return None


def replace_between(text:str, key1:str, key2:str, file:str) -> None:
    '''
    Replace lines with a given `text`, between the keywords `key1` and `key2`,
    in a given `file`.
    ```
    lines...
    key1
    text
    key2
    lines...
    ```
    '''
    delete_between(key1, key2, file)
    insert_under(text, key1, file)
    return None


def delete_between(key1:str, key2:str, file:str) -> None:
    '''
    Deletes the lines between two keywords in a given `file`.
    ```
    lines...
    key1
    (lines to be deleted)
    key2
    lines...
    ```
    '''
    file_path = get(file)
    with open(file_path, 'r') as f:
        lines = f.readlines()
    keep = []
    skip = False
    for line in lines:
        if key1 in line:
            skip = True
        if key2 in line:
            skip = False
        if not skip or key1 in line:
            keep.append(line)
    with open(file_path, 'w') as f:
        f.writelines(keep)
    return None


def correct_with_dict(file:str, fixing_dict:dict) -> None:
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

