import thoth as th

th.call.here()

file = 'sample.txt'
key1 = 'key 1'
key2 = 'key 2'
match_number = 1
regex = False

match = th.text.find_between(key1, key2, file, match_number, regex)
print(match)
