import thotpy as th

# We run here the script
th.call.here()

# Modify the parameters as you wish to test the function
key1 = 'key'
key2 = 'key'
file = 'sample.txt'
include_keys = True
match_number = 1
regex = False

# Find all the text in between the keywords
match = th.find.between(key1, key2, file, include_keys, match_number, regex)
print(match)
