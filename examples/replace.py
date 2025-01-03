import thotpy as th

# Run here
th.call.here()

# Modify this values to test the function
# This will modify sample.txt
filename = 'sample.txt'
key = 'key'
key_regex = r'key\s*\d*'  # 'key' AND whatever number
text = '!!!'
number_of_replacements = 1
regex = False
if regex:
    keyword = key_regex
else:
    keyword = key

# Replace the keywords
th.text.replace(filename, keyword, text, number_of_replacements, regex)
