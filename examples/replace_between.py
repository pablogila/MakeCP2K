import thotpy as th

th.call.here()

text = '!!!'
key1 = 'key 1'
key2 = 'key 2'
filename = 'sample.txt'
delete_keys = False
from_end = True
regex = False

th.text.replace_between(text, key1, key2, filename, delete_keys, from_end, regex)
