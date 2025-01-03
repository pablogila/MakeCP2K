import thotpy as th

# We run here the script
th.call.here()

# Modify the parameters as you wish to test the function
# This will modify sample.txt
filename='sample.txt'
key1 = 'key 1'
key2 = 'key 2'
delete_keys = False
from_end = False
regex = False

# If text='', it deletes the text in between!
th.text.replace_between(filename, key1, key2, '', delete_keys, from_end, regex)
