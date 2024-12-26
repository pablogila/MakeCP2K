import thoth as th

th.call.here()

file = 'sample.txt'
key = 'key'
key_regex = r'key\s*\d*'
text = '!!!testing!!!'
replace_last = False
skip_lines = 0
regex = False

if regex:
    keyword = key_regex
else:
    keyword = key

th.text.replace_under(text, keyword, file, replace_last, skip_lines, regex)
