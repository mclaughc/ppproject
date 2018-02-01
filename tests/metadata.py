import sys
sys.path.append("..")

from ppproject import Metadata

md = Metadata()
md.set_property("name", "foo")
md.set_property("vals", {'foo': 1, 'bar': 2})

print(md)

print(md.format_string("My name is %name% and my vals are %vals%."))
print(md.format_string("This is a percent symbol %%"))

