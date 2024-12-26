import thoth as th

th.call.here()

text = '!!!'
file = 'sample.txt'
key1 = 'key 1'
key2 = 'key 2'
last_match = False
regex = False

th.text.replace_between(text, key1, key2, file, last_match, regex)
