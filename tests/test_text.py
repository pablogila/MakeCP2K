import thotpy as th


folder = 'tests/samples/'
sample = folder + 'sample.txt'
sample_copy = folder + 'sample_copy.txt'


def test_insert_at():
    th.file.copy(sample, sample_copy)
    th.text.insert_at(filepath=sample_copy, text='MIDDLE', position=1)
    th.text.insert_at(filepath=sample_copy, text='START', position=0)
    th.text.insert_at(filepath=sample_copy, text='END', position=-1)
    with open(sample_copy, 'r') as f:
        assert f.read() == 'START\nline1\nMIDDLE\nline2\nline3\nline4\nline5\nline6\nline7\nline8\nline9\nEND'
    th.file.remove(sample_copy)


def test_insert_under():
    th.file.copy(sample, sample_copy)
    th.text.insert_under(filepath=sample_copy, key='5', text='!!!', skips=0)
    with open(sample_copy, 'r') as f:
        assert f.read() == 'line1\nline2\nline3\nline4\nline5\n!!!\nline6\nline7\nline8\nline9'
    th.file.copy(sample, sample_copy)
    th.text.insert_under(filepath=sample_copy, key='5', text='!!!', skips=-1)
    with open(sample_copy, 'r') as f:
        assert f.read() == 'line1\nline2\nline3\nline4\n!!!\nline5\nline6\nline7\nline8\nline9'
    th.file.copy(sample, sample_copy)
    th.text.insert_under(filepath=sample_copy, key=r'l[a-z]*5', text='!!!', regex=True)
    with open(sample_copy, 'r') as f:
        assert f.read() == 'line1\nline2\nline3\nline4\nline5\n!!!\nline6\nline7\nline8\nline9'
    th.file.remove(sample_copy)


def test_replace():
    th.file.copy(sample, sample_copy)
    th.text.replace(filepath=sample_copy, key='line5', text='!!!')
    with open(sample_copy, 'r') as f:
        assert f.read() == 'line1\nline2\nline3\nline4\n!!!\nline6\nline7\nline8\nline9'
    th.file.remove(sample_copy)


def test_replace_line():
    th.file.copy(sample, sample_copy)
    th.text.replace_line(filepath=sample_copy, key='line5', text='!!!')
    with open(sample_copy, 'r') as f:
        assert f.read() == 'line1\nline2\nline3\nline4\n!!!\nline6\nline7\nline8\nline9'
    th.file.copy(sample, sample_copy)
    th.text.replace_line(filepath=sample_copy, key='line5', text='')
    with open(sample_copy, 'r') as f:
        assert f.read() == 'line1\nline2\nline3\nline4\nline6\nline7\nline8\nline9'
    th.file.remove(sample_copy)


def test_replace_between():
    th.file.copy(sample, sample_copy)
    th.text.replace_between(filepath=sample_copy, key1='line4', key2='line7', text='!!!')
    with open(sample_copy, 'r') as f:
        assert f.read() == 'line1\nline2\nline3\nline4\n!!!\nline7\nline8\nline9'
    th.file.remove(sample_copy)


def test_remove_between():
    th.file.copy(sample, sample_copy)
    th.text.replace_between(filepath=sample_copy, key1='line4', key2='line7', text='')
    with open(sample_copy, 'r') as f:
        assert f.read() == 'line1\nline2\nline3\nline4\nline7\nline8\nline9'
    th.file.remove(sample_copy)


def test_delete_under():
    th.file.copy(sample, sample_copy)
    th.text.delete_under(filepath=sample_copy, key='5')
    with open(sample_copy, 'r') as f:
        assert f.read() == 'line1\nline2\nline3\nline4\nline5'
    th.file.remove(sample_copy)


def test_correct_with_dict():
    correct = {'line1': 'text', 'line5': ''}
    th.file.copy(sample, sample_copy)
    th.text.correct_with_dict(filepath=sample_copy, replaces=correct)
    with open(sample_copy, 'r') as f:
        assert f.read() == 'text\nline2\nline3\nline4\n\nline6\nline7\nline8\nline9'
    th.file.remove(sample_copy)

