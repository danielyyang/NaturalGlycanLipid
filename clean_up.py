import os
import glob

keep_lib = ['sugar_sequence.py', 'sugar_utils.py', 'chemical_dictionaries.py', '__init__.py']
keep_scripts = ['main_pipeline.py', '__init__.py']
keep_tests = ['test_linkage.py', '__init__.py']

base_dir = r"d:\Glycan_Database"

for file in glob.glob(os.path.join(base_dir, "lib", "*.py")):
    if os.path.basename(file) not in keep_lib:
        os.remove(file)

for file in glob.glob(os.path.join(base_dir, "scripts", "*.py")):
    if os.path.basename(file) not in keep_scripts:
        os.remove(file)

for file in glob.glob(os.path.join(base_dir, "tests", "*.py")):
    if os.path.basename(file) not in keep_tests:
        os.remove(file)

if os.path.exists(os.path.join(base_dir, "debug_rs.py")):
    os.remove(os.path.join(base_dir, "debug_rs.py"))
