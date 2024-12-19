FPCore optimization of base numeric functions.

## Workflow
TODO: make this more automated, and run it on a server for a whole week or something

- generate python equations using sympy, for example `quadratic_sym.py` makes the
  four branches of the complex quadratic equation
- use `sympy_to_fpcore.py` to transform python equations into FPCore equations
- use `racket -l herbie` to optimize FPCore equations, this will produce a new
  FPCore equation, which can be fed back into herbie again, over and over
- use `FPBench/transform.rkt` to simplify expressions
- use `FPBench/export.rkt` to get a rust expression

Example:
```bash
racket -l herbie -- improve --timeout 999 --num-iters 16 --num-points 512 gen1.fpcore gen2_1.fpcore
racket -y FBench/transform.rkt --cse gen2_1.fpcore gen2_1_cse.fpcore
```

Note that the web version of Herbie is generally better and for some reason provides better
results, but it is less automatable.
