import thotpy as th

# We run here the script
th.call.here()

# Modify the parameters as you wish to test the function
# This will modify sample.txt
file = 'sample.txt'
key = 'key'
key_regex = r'key\s*\d*'
match = 2
skip_lines = 0
regex = False
if regex:
    keyword = key_regex
else:
    keyword = key

# Deletes ALL the content under the specified match, up to the end of the file.
th.text.delete_under(file, keyword, match, skip_lines, regex)
