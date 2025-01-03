import thotpy as th

# Run here
th.call.here()

# Modify the values to test the function
# This will modify sample.txt
# You could also delete the text in between by setting text=''
filename = 'sample.txt'
key1 = 'key 1'
key2 = 'key 2'
text = '!!!'
delete_keys = False
from_end = True
regex = False

# Replace between the keywords
th.text.replace_between(filename, key1, key2, text, delete_keys, from_end, regex)
