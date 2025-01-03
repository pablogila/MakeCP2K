from thotpy import find


sample = 'tests/samples/sample.txt'
'''For this path to be detected properly, the `pytest` must be executed from `ThotPy`!'''


def test_lines():
    assert find.lines(filepath=sample, key='5', matches=0, additional=0, split=False, regex=False) == ['line5']
    assert find.lines(filepath=sample, key='5', matches=-1, additional=0, split=False, regex=False) == ['line5']
    assert find.lines(filepath=sample, key='5', matches=1, additional=1, split=True, regex=False) == ['line5', 'line6']
    assert find.lines(filepath=sample, key='5', matches=1, additional=-1, split=True, regex=False) == ['line4', 'line5']
    assert find.lines(filepath=sample, key='5', matches=0, additional=0, split=False, regex=True) == ['line5']
    assert find.lines(filepath=sample, key='line', matches=1, additional=0, split=False, regex=True) == ['line1']
    assert find.lines(filepath=sample, key='line', matches=-1, additional=0, split=False, regex=True) == ['line9']
    assert find.lines(filepath=sample, key='line', matches=-2, additional=1, split=True, regex=False) == ['line8', 'line9', 'line9']
    assert find.lines(filepath=sample, key='line', matches=-2, additional=1, split=True, regex=True) == ['line8', 'line9', 'line9']


def test_between():
    assert find.between(filepath=sample, key1='5', key2='7', include_keys=False, match=1, regex=False) == 'line6'
    assert find.between(filepath=sample, key1='5', key2='7', include_keys=True, match=1, regex=False) == 'line5\nline6\nline7'

