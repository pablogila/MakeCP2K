import thoth as th

th.call.here()

file = 'sample.txt'
key = 'key'
key_regex = r'key\s*\d*'
text = '!!!testing!!!'
number_of_replacements = -1
skip_lines = 0
regex = False

if regex:
    keyword = key_regex
else:
    keyword = key

th.text.replace_line(text, keyword, file, number_of_replacements, skip_lines, regex)
