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


def test_extract_element():
    string = '  element I Lead Pb Nitrogen H2, Xx2 fake element, O Oxygen, He4 isotope Ag Element '
    assert th.extract.element(text=string, index=0) == 'I'
    assert th.extract.element(text=string, index=1) == 'Pb'
    assert th.extract.element(text=string, index=2) == 'H2'
    assert th.extract.element(text=string, index=3) == 'O'
    assert th.extract.element(text=string, index=4) == 'He4'
    assert th.extract.element(text=string, index=5) == 'Ag'

