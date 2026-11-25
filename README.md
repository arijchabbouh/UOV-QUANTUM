# Simplified UOV (Unbalanced Oil and Vinegar) signatures

Course project: toy implementation of the UOV multivariate post-quantum signature
over GF(16) and GF(31), compared with the multivariate schemes of liboqs.

> Toy parameters (o = 6..10, v = 12..20): for understanding and measuring only, NOT secure.

## Setup
    python -m venv .venv
    # Windows: .venv\Scripts\activate    macOS/Linux: source .venv/bin/activate
    pip install -r requirements.txt

## Run the tests
    pytest

## Project structure
| File | Content | Status |
|---|---|---|
| field.py | GF(16) and GF(31) with arithmetic tables | done |
| linalg.py | Gaussian elimination, matrix operations | todo |
| keygen.py | secret key O, public key P | todo |
| sign.py | hashing, signing with retry | todo |
| verify.py | verification | todo |
| benchmark.py | timings and sizes, liboqs comparison | todo |
| tests/ | pytest tests | |
| results/ | raw benchmark data and plots | |
