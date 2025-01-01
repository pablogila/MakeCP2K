import thotpy as th

th.call.here()

text = '!!!'
key1 = 'key 1'
key2 = 'key 2'
reverse_match = True
file = 'sample.txt'
regex = False

th.text.replace_between(text, key1, key2, file, reverse_match, regex)
