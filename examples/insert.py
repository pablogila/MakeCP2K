import thotpy as th

th.call.here()

file = 'sample.txt'
key = 'key'
key_regex = r'key\s*\d*'
text = '!!!'
number_of_inserts = 1
skip_lines = 2
regex = False

if regex:
    keyword = key_regex
else:
    keyword = key

th.text.insert_under(text, key, file, number_of_inserts, skip_lines, regex)
