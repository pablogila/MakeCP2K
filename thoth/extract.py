'''
# Description
Functions to extract data from raw text strings.

# Index
- `number()`
- `string()`
- `column()`

---
'''


import re


def number(
        text:str,
        name:str=''
    ) -> float:
    '''
    Extracts the float value of a given `name` variable from a raw `text`.\n
    Example:
    ```python
    >>> text = 'energy =   500.0 Ry'
    >>> thoth.extract.number(text, 'energy')
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
    >>> thoth.extract.string(text, 'energy', 'were')
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
    result.strip()
    if strip:
        result = result.strip("'")
        result = result.strip('"')
        result.strip()
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

