import thotpy as th

th.call.here()

file = 'sample.txt'
key = 'key'
key_regex = r'key\s*\d*'
text = '!!!'
inserts = 1  # Number of inserts
skips = 2  # Number of lines to skip
regex = False

if regex:
    keyword = key_regex
else:
    keyword = key

th.text.insert_under(text, key, file, inserts, skips, regex)
