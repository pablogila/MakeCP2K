import thotpy as th

# We run here the script
th.call.here()

# Modify the parameters as you wish to test the function
# This will modify sample.txt
file = 'sample.txt'
text = '!!!'
position = -1

# Insert the text at the specified position
th.text.insert_at(file, text, position)
