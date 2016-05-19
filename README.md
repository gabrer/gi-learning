# GI-LEARNING LIBRARY

`GI-LEARNING LIBRARY` is a C++ framework for [grammatical inference](https://en.wikipedia.org/wiki/Grammar_induction).

Grammar Induction algorithms:
- RPNI
- EDSM
- Blue*
- L*

Easy instructions:
1. Go inside the project folder "GI-learning"

2. ```
  make distclean
  autoreconf -i
  ./configure
  make
  cd src/
  ./giLearning ../examples/examples_big.txt ../examples/lstar.txt ```

If you obtain some error from the "configure" command, resolve it installing the necessary libraries or tools.
