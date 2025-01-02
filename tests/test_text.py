import thotpy as th


folder = 'tests/samples/'
sample = folder + 'sample.txt'
sample_copy = folder + 'sample_copy.txt'


def test_insert_at():
    th.file.copy(sample, sample_copy)
    th.text.insert_at(text='MIDDLE', filename=sample_copy, position=1)
    th.text.insert_at(text='START', filename=sample_copy, position=0)
    th.text.insert_at(text='END', filename=sample_copy, position=-1)
    with open(sample_copy, 'r') as f:
        assert f.read() == 'START\nline1\nMIDDLE\nline2\nline3\nline4\nline5\nline6\nline7\nline8\nline9\nEND'
    th.file.remove(sample_copy)


def test_insert_under():
    th.file.copy(sample, sample_copy)
    th.text.insert_under(text='!!!', key='5', filename=sample_copy, skips=0)
    with open(sample_copy, 'r') as f:
        assert f.read() == 'line1\nline2\nline3\nline4\nline5\n!!!\nline6\nline7\nline8\nline9'
    th.file.copy(sample, sample_copy)
    th.text.insert_under(text='!!!', key='5', filename=sample_copy, skips=-1)
    with open(sample_copy, 'r') as f:
        assert f.read() == 'line1\nline2\nline3\nline4\n!!!\nline5\nline6\nline7\nline8\nline9'
    th.file.copy(sample, sample_copy)
    th.text.insert_under(text='!!!', key=r'l[a-z]*5', filename=sample_copy, regex=True)
    with open(sample_copy, 'r') as f:
        assert f.read() == 'line1\nline2\nline3\nline4\nline5\n!!!\nline6\nline7\nline8\nline9'
    th.file.remove(sample_copy)
