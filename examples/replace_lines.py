import thotpy as th

# Run here
th.call.here()

# Modify these values to test the function
# This will modify sample.txt
file = 'sample.txt'
key = 'key'
key_regex = r'key\s*\d*'
text = '!!!'
replacements = 0 # Number of matches to replace
skips = 0 # Number of lines to skip
additional = 0  # Number of additional lines
regex = False
if regex:
    keyword = key_regex
else:
    keyword = key

# Replace the given line or line plus additional lines
th.text.replace_line(text, keyword, file, replacements, skips, additional, regex)
