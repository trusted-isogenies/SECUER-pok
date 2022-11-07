# Trusted Curve

<img
 width="33%"
 align="right"
 src="https://user-images.githubusercontent.com/84067835/200272479-3af325ca-7602-4785-b493-fe5152bdc9b8.png"/>


**Generating a trusted curve as the start of isogeny-based protocols**

Implementations of the proof of isogeny knowledge from the paper
*[Supersingular Curves You Can Trust](https://ia.cr/2022/1469)*.


## Building on Linux

```sh
cd c-impl
make
make test
```

## Building on Apple M1

```sh
cd c-impl
make ARCH=M1
make test ARCH=M1
```

## Generating and Verifying a Proof

The `prove_xxx` executables can be run with the `--initial` argument to start from the curve with j‑invariant 1728. If no argument is passed, the program expects a starting curve on `stdin`.

The `verify_xxx` executables expect the output of the corresponding `prove_xxx` executable on `stdin`. In other words, the following sequence of invocations is typical:

```sh
./prove_434 --initial > proof0.txt
./verify_434 < proof0.txt | tail -n1 > curve1.txt
./prove_434 < curve1.txt > proof1.txt
./verify_434 < proof1.txt | tail -n1 > curve2.txt
./prove_434 < curve2.txt > proof2.txt
./verify_434 < proof2.txt | tail -n1 > curve3.txt
# ...
```


## Verification in Sage (slow!)

The `verify.sage` script reads a proof on `stdin` and verifies it. The program expects one of the arguments `--p434`, `--p503`, `--p610`, `--p751` to specify the parameter set. It requires `sage >= 9.7` and the `pycryptodome` package.

