from kmtools import structure_tools

s = structure_tools.fetch_structure('4dkl')
structure_tools.process_structure(s)

print(s)

