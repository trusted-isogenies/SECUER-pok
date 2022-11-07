# Trusted Curve

<img
 width="33%"
 align="right"
 src="https://user-images.githubusercontent.com/552961/152566766-66e227ff-e239-454e-9705-3e3d7f53b104.png"/>


**Generating a trusted curve as the start of isogeny-based protocols**

Implementations of the proof of isogeny knowledge from "Supersingular Curves You Can Trust".


## Building on Linux

```
cd c-impl;
make;
make test;
```

## Building on Apple M1

```
cd c-impl;
make ARCH=M1;
make test ARCH=M1;
```

## Generating a Proof

The `prove_xxx` executables can be run with the `--initial` argument to start from the curve with j-invariant 1728. If no argument is passed, the program expects a starting curve on `stdin`.

## Verification in Sage

The `verify.sage` script reads a proof on `stdin` and verifies it. The program expects one of the following arguments `--p434`, `--p503`, `--p610`, `--p751` to determine the parameter set. It requires `sage >= 9.7` and the `pycryptodome` package.
