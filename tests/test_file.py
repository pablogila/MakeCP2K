'''For this path to be detected properly, the `pytest` must be executed from `ThotPy`!'''

import thotpy as th


folder = 'tests/'
sample = folder + 'sample.txt'
sample_copy = folder + 'sample_copy.txt'
sample_copy_2 = folder + 'sample_copy_2.txt'
sample_ok = folder + 'sample_ok.txt'
sample_ok_2 = folder + 'sample_ok_2.txt'


def test_get():
    try:
        assert th.file.get(sample) != None
        assert True
    except FileNotFoundError:
        assert False
    try:
        th.file.remove(sample_copy)
    except:
        pass
    try:
        th.file.get(sample_copy)
        assert False
    except FileNotFoundError:
        assert True


def test_template():
    try:
        th.file.get(sample_copy)
        assert False
    except FileNotFoundError:
        assert True
    th.file.from_template(old=sample, new=sample_copy, comment='!!!', fixing_dict={'line':''})
    with open(sample_copy, 'r') as f:
        content = f.read()
        assert content == '!!!\n1\n2\n3\n4\n5\n6\n7\n8\n9'

    th.file.remove(sample_copy)
    try:
        th.file.get(sample_copy)
        assert False
    except FileNotFoundError:
        assert True


def test_copy():
    try:
        th.file.remove(sample_copy)
    except:
        pass
    th.file.copy(sample, sample_copy)
    try:
        assert th.file.get(sample_copy) != None
        assert True
    except FileNotFoundError:
        assert False
    try:
        th.file.move(sample_copy, sample_copy_2)
        assert True
    except:
        assert False

    try:
        th.file.remove(sample_copy_2)
        assert True
    except:
        assert False
    try:
        th.file.remove(sample_copy)
        assert False
    except:
        assert True


def test_rename():
    try:
        th.file.remove(sample_copy)
        th.file.remove(sample_ok)
    except:
        pass
    th.file.copy(sample, sample_copy)
    th.file.rename(old='copy', new='ok', folder=folder)
    try:
        th.file.remove(sample_ok)
        assert True
    except:
        assert False
    try:
        th.file.remove(sample_copy)
        assert False
    except:
        assert True

