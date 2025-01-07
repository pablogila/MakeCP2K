import thotpy as th

# We run here the script
th.call.here()

# Modify the parameters as you wish to test the function
file = 'sample.txt'
key1 = 'key2'
key2 = 'key3'
include_keys = True
match_number = 1
regex = True

# Find all the text in between the keywords
match = th.find.between(file, key1, key2, include_keys, match_number, regex)
print(match)
