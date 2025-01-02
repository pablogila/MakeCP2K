import thotpy as th

# We run here the script
th.call.here()

# Modify the parameters as you wish to test the function
# This will modify sample.txt
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

# Insert the text under the given keyword
th.text.insert_under(text, key, file, inserts, skips, regex)
