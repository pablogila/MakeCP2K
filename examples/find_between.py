import thoth as th

th.call.here()

key1 = 'key'
key2 = 'key'
file = 'sample.txt'
include_keys = True
match_number = 1
regex = False

match = th.text.find_between(key1, key2, file, include_keys, match_number, regex)
print(match)
