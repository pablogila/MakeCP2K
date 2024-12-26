import thoth as th

th.call.here()

file = 'sample.txt'
key = 'key'
key_regex = r'key\s*\d*'
match = 1
skip_lines = 0
regex = False

if regex:
    keyword = key_regex
else:
    keyword = key

th.text.delete_under(keyword, file, match, skip_lines, regex)
