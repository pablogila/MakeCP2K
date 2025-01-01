import thotpy as th

th.call.here()

filename='sample.txt'
key1 = 'key 1'
key2 = 'key 2'
delete_keys = False
from_end = False
regex = False

th.text.replace_between('', key1, key2, filename, delete_keys, from_end, regex)
