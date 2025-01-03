'''
# Description
Functions to extract data from raw text strings.

# Index
- `number()`
- `string()`
- `column()`
- `coords()`
- `element()`

---
'''


import re
import maatpy as mt


def number(
        text:str,
        name:str=''
    ) -> float:
    '''
    Extracts the float value of a given `name` variable from a raw `text`.\n
    Example:
    ```python
    >>> text = 'energy =   500.0 Ry'
    >>> thotpy.extract.number(text, 'energy')
    500.0  # float output
    ```
    '''
    if text == None:
        return None
    pattern = re.compile(rf"{name}\s*[:=]?\s*(-?\d+(?:\.\d+)?(?:[eEdD][+\-]?\d+)?)")
    match = pattern.search(text)
    if match:
        return float(match.group(1))
    return None
    

def string(
        text:str,
        name:str='',
        stop:str='',
        strip:bool=True
    ) -> str:
    '''
    Extracts the `text` value of a given `name` variable from a raw string.
    Stops before an optional `stop` string.
    Removes leading and trailing commas by default, change this with `strip=False`.\n
    Example:
    ```python
    >>> text = 'energy =   500.0 Ry were calculated'
    >>> thotpy.extract.string(text, 'energy', 'were')
    '500.0 Ry'  # String output
    ```
    '''
    pattern = re.compile(rf"{name}\s*[:=]?\s*(.*)")
    if stop:
        pattern = re.compile(rf"{name}\s*[:=]?\s*(.*)(?={stop})")
    match = re.search(pattern, text)
    if not match:
        return None
    result = str(match.group(1))
    result = result.strip()
    if strip:
        result = result.strip("'")
        result = result.strip('"')
        result = result.strip()
    return result


def column(
        text:str,
        column:int
    ) -> float:
    '''
    Extracts the desired float `column` of a given `string`.
    '''
    if text is None:
        return None
    columns = text.split()
    pattern = r'(-?\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?)'
    if column < len(columns):
        match = re.match(pattern, columns[column])
        if match:
            return float(match.group(1))
    return None


def coords(text:str) -> list:
    '''
    Returns a list with the float coordinates expressed in a given `text` string.
    '''
    if text is None:
        return None
    columns = re.split(r'[,\s]+', text.strip())
    pattern = r'(-?\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?)'
    matches = []
    for column in columns:
        match = re.match(pattern, column)
        if match:
            matches.append(float(match.group(1)))
    return matches


def element(text:str, index:int=0) -> str:
    '''
    Extract a chemical element from a raw `text` string.
    If there are several elements, you can return a specific `index` match (positive, 0 by default).
    '''
    if text is None:
        return None
    columns = re.split(r'[,\s]+', text.strip())
    pattern = r'\s*([A-Z][a-z]*[a-z]*)\s*'
    matches = []
    for column in columns:
        match = re.match(pattern, column)
        if match:
            matches.append(str(match.group(1)))
    # We have a list with possible matches. Let's determine which are actual elements.
    counter = 0
    last_match = None
    for possible_element in matches:
        possible_element = possible_element.strip()
        if possible_element in mt.atom.keys():
            if counter == abs(index):
                return possible_element
            last_match = possible_element
            counter += 1
    return last_match

