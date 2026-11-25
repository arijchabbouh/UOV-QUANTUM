# Lets the tests import modules from the project root (field.py, linalg.py, ...)
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
