import thotpy as th


def test_extract_number():
    assert th.extract.number(' test = 123 number ', 'test') == 123.0
    assert th.extract.number(' test 123 number ', 'test') == 123.0


def test_extract_string():
    assert th.extract.string(text=' test = "hello" stop ', name='test', stop='stop', strip=True) == 'hello'
    assert th.extract.string(text=' test "hello" stop ', name='test', stop='stop', strip=True) == 'hello'
    assert th.extract.string(text=" test 'hello' stop ", name='test', stop='stop', strip=True) == 'hello'
    assert th.extract.string(text=" test 'hello' stop ", name='test', stop='stop', strip=False) == "'hello'"


def test_extract_column():
    assert th.extract.column(' 123 456.5 789  ', 2) == 789


def test_extract_coords():
    assert th.extract.coords('coordinates: 1.0, 2.0 and 3 these were the coordinates') ==[1.0, 2.0, 3.0]

