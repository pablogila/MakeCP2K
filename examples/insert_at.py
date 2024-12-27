import thoth as th

th.call.here()

text = '!!!'
file = 'sample.txt'
position = -1

th.text.insert_at(text, file, position)
