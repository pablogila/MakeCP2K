import thotpy as th

th.call.here()

file = 'sample.txt'
key = 'key'
key_regex = r'key\s*\d*'  # 'key' AND whatever number
text = '!!!'
number_of_replacements = 1
regex = False

if regex:
    keyword = key_regex
else:
    keyword = key

th.text.replace(text, keyword, file, number_of_replacements, regex)
