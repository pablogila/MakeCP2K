import thoth as th

th.call.here()

file = 'sample.txt'
key1 = 'key 1'
key2 = 'key 2'
delete_keys = False
last_match = False
regex = False

th.text.delete_between(key1, key2, file, delete_keys, last_match, regex)
