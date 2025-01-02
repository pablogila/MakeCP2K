'''For this path to be detected properly, the `pytest` must be executed from `ThotPy`!'''

import thotpy as th


folder = 'tests/samples/'
sample = folder + 'sample.txt'
sample_copy = folder + 'sample_copy.txt'
sample_copy_2 = folder + 'sample_copy_2.txt'
sample_ok = folder + 'sample_ok.txt'
sample_ok_2 = folder + 'sample_ok_2.txt'


def test_get():
    # Clean from previous tests
    try:
        th.file.remove(sample_copy)
    except:
        pass
    # Finds an existing file
    try:
        assert th.file.get(sample) != None
        assert True
    except FileNotFoundError:
        assert False
    # Does not find a non-existing file
    try:
        th.file.get(sample_copy)
        assert False
    except FileNotFoundError:
        assert True
    # get_list, 'tests/sample.txt' in 'fullpath/tests/sample.txt'
    try:
        assert sample in th.file.get_list(folder, filters='sample')[0]
    except:
        assert False


def test_template():
    try:
        th.file.remove(sample_copy)
    except:
        pass
    try:
        th.file.from_template(old=sample, new=sample_copy, comment='!!!', fixing_dict={'line':''})
        with open(sample_copy, 'r') as f:
            content = f.read()
            assert content == '!!!\n1\n2\n3\n4\n5\n6\n7\n8\n9'
    except:
        assert False
    try:
        th.file.remove(sample_copy)
    except:
        pass


def test_copy():
    try:
        th.file.remove(sample_copy)
        th.file.remove(sample_copy_2)
    except:
        pass
    # Copy files
    try:
        th.file.copy(sample, sample_copy)
        assert th.file.get(sample_copy) != None
    except FileNotFoundError:
        assert False
    # Move files
    try:
        th.file.move(sample_copy, sample_copy_2)
        assert th.file.get(sample_copy_2) != None
    except:
        assert False
    try:
        th.file.get(sample_copy)
        assert False
    except FileNotFoundError:
        assert True
    # Remove
    try:
        th.file.remove(sample_copy_2)
    except:
        assert False
    try:
        th.file.get(sample_copy_2)
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


def test_folders():
    try:
        th.file.remove(sample_copy)
        th.file.remove(sample_copy_2)
        th.file.remove(sample_ok)
        th.file.remove(sample_ok_2)
    except:
        pass
    th.file.copy_to_subfolders(folder=folder, extension='.txt', strings_to_delete=['.txt'])
    try:
        assert th.file.get_list(folder=folder+'sample', abspath=False) == ['sample.txt']
    except:
        assert False
    # Check that the folder is deleted
    th.file.remove(folder+'sample')
    try:
        th.file.get_list(folder+'sample')
        assert False
    except FileNotFoundError:
        assert True

