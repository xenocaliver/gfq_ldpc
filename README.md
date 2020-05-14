
# gfq_ldpc
Non-binary LDPC encode and decode simulation programs

## Preparation
Install following tools.

- gcc, g++
- boost library
- git
- make
- cmake
- [msgpack](https://github.com/msgpack/msgpack-c)

## Downloading the simulator
Clone simulator's repository.

```sh
git clone git@github.com:xenocaliver/gfq_ldpc.git
```

## Build the simulator
Build the simulator as follows.
```sh
cd gfq_ldpc
mkdir build
cd build
cmake ..
make
```
## Make generating matrix
Before running the simulator, you must create generating matrix.
You can convert parity check matrix alist file into generating matrix as follows. 
```sh
./make-gen <parity check alist file name> <generating matrix file name>
```

## Run the simulator
Then, you can run the simulator as follows.

```sh
./gfq_simulator <alist file name> <generating matrix file name> <number of transmission> <sigma for AWGN channel> <Sum-Product iteration limit>
```

## Illustration for implementation
### Finite field calculation
This simulator uses [galois++](https://github.com/wkjarosz/galois) for finite field calculation. [galois++](https://github.com/wkjarosz/galois) wraps finite field element as a class and overrides operators, for example, +,-,* and /. Therefore,
a programmer can code finite field calculations like as mathematical formulae.
Thus parity check and encoding part of the simulator code are very simple.

### Encoding
The simulator's generating matrix is constructed as follows.
Let $\bm{c}$ a codeword(column vector), $H$ a parity check matrix. If $\bm{c}$ is a codeword, 
$$
H\bm{c}=\bm{0}
$$
holds. And we assume we can split $\bm{c}$ into two parts $\bm{p}$ and $\bm{s}$. $\bm{p}$ represents a parity check symbol part and $\bm{s}$ represents input symbol part.
And we devide $H$ into two part i.e.
$$
H=\left(A\vert B\right).
$$
$A$'s number of columns coinides with $\bm{p}$'s dimension. $B$'s number of columns coincides with $\bm{s}$'s dimension. So,
$$
A\bm{p}+B\bm{s}=\bm{0}
$$
holds. Therefore one can calculate $\bm{p}$ by means of following expression:
$$
\bm{p}=-A^{-1}B\bm{s}.
$$
Next, one can calculate $A^{-1}$ by means of $LU$ decompsition. However, one can apply $LU$ decomposition only if $A$ is full rank. However $A$ is not always full rank. So, `make-gen` permuates $H$'s columns and make $A$ be full rank. A column permuation corresponds to a matrix. Let $Q_{i}(i=1,\ldots, m)$ be a matrix corresponding to a column exchange. So we can write
$$
\bm{c}=\left(
    \begin{array}{c}
    -A^{-1}B\\
    I
    \end{array}
    \right)Q_{1}\cdots Q_{m}\bm{s}
$$
where $I$ is an identity matrix and we used $Q_{i}^{-1}=Q_{i}$. Therefore we obtain generating matrix $G$ as follows:
$$
G=\left(
    \begin{array}{c}
    -A^{-1}B\\
    I
    \end{array}
    \right)Q_{1}\cdots Q_{m}.
$$

## Transmitting a symbol over AWGN channel
We assumed that channel is AWGN channel with standard deviation $\sigma$ and assumed
finite field be $\mathbb{F}_{q}, q=2^{m}$. According to [Davey](https://ieeexplore.ieee.org/document/706440), each bit of a symbol is modulated to BPSK and bitwise prior probability for bit $k$ is given by
$$
p_{k}=\frac{1}{1+\exp[2s_{k}y_{k}/\sigma^{2}]}
$$
where $s_{k}$ is BPSK modulated value of the bit $k$ and $y_{k}$ is a observed value of bit $k$. Therefore prior probability according to symbol $g\in\mathbb{F}_{q}$ is given by
$$
f(g) = \prod_{k=1}^{m}\frac{1}{1+\exp[2s_{k}y_{k}/\sigma^{2}]}.
$$

## Decoding

<script async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS_CHTML"></script>
<script type="text/x-mathjax-config">
 MathJax.Hub.Config({
 tex2jax: {
 inlineMath: [["\\(","\\)"] ],
 displayMath: [ ['$$','$$'], ["\\[","\\]"] ]
 }
 });
</script>