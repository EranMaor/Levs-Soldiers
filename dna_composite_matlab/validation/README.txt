Cross-check of app_calc.m against the Python implementation:
  python3 gen.py      -> writes cases.txt (random channel outputs + Python posteriors)
  MATLAB/Octave: cmp  -> prints max |MATLAB - Python| per case (expected ~1e-16)
