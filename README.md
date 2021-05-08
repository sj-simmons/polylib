# polylib

Provides classes for easily computing with polynomials over a ring or over a field.
After installation, see example usage by typing, in your Python interpreter,
```python
>>> import polylib
>>> help(polylib.Polynomial)
```
or
```python
>>> help(polylib.FPolynomial)
```

Installing this package also installs an executable program, **bernoulli**, that uses
polynomial arithmetic to compute the generating series for bernoulli numbers. Type
```shell
bernoulli
```
or
```shell
~/.local/bin/bernoulli
```
at your command-line for usage details.

## Installation

From your command line:

```shell
pip install polylib --user
```
or, for the latest dev version:
```shell
pip install git+https://github.com/sj-simmons/polylib.git --user
```
