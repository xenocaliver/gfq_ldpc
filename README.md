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
- [galois++](https://github.com/wkjarosz/galois)
- [fftw](https://www.fftw.org)

[galois++](https://github.com/wkjarosz/galois) and [fftw](https://www.fftw.org) are downloaded and compiled during this simulator's compilation process.
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
Let $\boldsymbol{c}$ a codeword(column vector), $H$ a parity check matrix. If $\boldsymbol{c}$ is a codeword, 
$$
H\boldsymbol{c}=\boldsymbol{0}
$$
holds. And we assume we can split $\boldsymbol{c}$ into two parts $\boldsymbol{p}$ and $\boldsymbol{s}$. $\boldsymbol{p}$ represents a parity check symbol part and $\boldsymbol{s}$ represents input symbol part.
And we devide $H$ into two part i.e.

$$
H=\left(A\vert B\right).
$$

$A$'s number of columns coinides with $\boldsymbol{p}$'s dimension. $B$'s number of columns coincides with $\boldsymbol{s}$'s dimension. So,

$$
A\boldsymbol{p}+B\boldsymbol{s}=\boldsymbol{0}
$$

holds. Therefore one can calculate $\boldsymbol{p}$ by means of following expression:

$$
\boldsymbol{p}=-A^{-1}B\boldsymbol{s}.
$$

Next, one can calculate $A^{-1}$ by means of $LU$ decompsition. However, one can apply $LU$ decomposition only if $A$ is full rank. However $A$ is not always full rank. So, `make-gen` permuates $H$'s columns and make $A$ be full rank. A column permuation corresponds to a matrix. Let $Q_{i}(i=1,\ldots, l)$ be a matrix corresponding to a column exchange. So we can write

$$
\boldsymbol{c}=\left( \begin{array}{c} -A^{-1}B\\ I \end{array} \right)Q_{1}\cdots Q_{l}\boldsymbol{s}
$$

where $I$ is an identity matrix and we used $Q_{i}^{-1}=Q_{i}$. Therefore we obtain generating matrix $G$ as follows:

$$
G=\left( \begin{array}{c} -A^{-1}B\\ I \end{array} \right)Q_{1}\cdots Q_{l}.
$$

### Transmitting a symbol over AWGN channel
We assumed that channel is AWGN channel with standard deviation $\sigma$ and assumed
finite field be $\mathbb{F}_{q}, q=2^{M}$. According to [Davey](https://ieeexplore.ieee.org/document/706440), each bit of a symbol is modulated to BPSK and bitwise prior probability for bit $i$ is given by

$$
p_{i}=\frac{1}{1+\exp[2s_{i}\vert y_{i}\vert/\sigma^{2}]}
$$

where $s_{i}$ is BPSK modulated value of the bit $i$ and $y_{i}$ is a observed value of bit $i$. Therefore prior probability according to symbol $g\in\mathbb{F}_{q}$ is given by

$$
f(g) = \prod_{i=1}^{M}\frac{1}{1+\exp[2s_{i}\vert y_{i}\vert/\sigma^{2}]}.
$$

### Decoding
We employ probability region Sum-Product algorithm for decoding according to [Davey](https://ieeexplore.ieee.org/document/706440).
Let $q_{mn}(g)$ be a veriable node $n$ to factor node message $m$ where $g$ is a finite field element's value. Let $r_{nm}$ be a factor node $m$ to variable node $n$ message.

### Initialization
For all variable nodes, initialize messages as follows:

$$
q_{mn}(g)=f_{n}(g).
$$

### Factor to variable node message update
Factor to variable node message update processes are given by

$$
r_{nm}(g)=\sum_{x_{i_{1}}}\cdots\sum_{x_{i_{d_{c}-1}}}\boldsymbol{1}[\sum_{k=1}^{d_{c}}h_{mi_{k}}x_{i_{k}}=0|x_{n} = g]\prod_{n^{\prime}\in∂ m\backslash n}q_{mn^{\prime}}(x_{i_{k}})
$$

where $\boldsymbol{1}$ denotes an indicator function i.e.

$$
\boldsymbol{1}[A]=
\begin{cases}
1& (A\quad\mathrm{is\quad true})\\
0& (A\quad\mathrm{is\quad false}).
\end{cases}
$$

However, this update rule's computing complexity is very large and is not practical. Threfore, we use Fourier transformaion according to [Hong](https://ieeexplore.ieee.org/document/6113595/). Another update rule is shown as follows. At first, we apply Fourier transformation $\mathscr{F}$ to variable to factor messages as follows:

$$
\boldsymbol{Q}_{mn} = \mathscr{F}\left[\boldsymbol{q}_{mn}\right]
$$

where \(\boldsymbol{Q}_{mn} = (Q_{mn}(0),Q_{mn}(1),\ldots,Q_{mn}(2^{M}-1))$ and $\boldsymbol{q}_{mn} = (q_{mn}(0), q_{mn}(1),\ldots, q_{mn}(2^{M}-1))\). And we update factor to variable messages in frequency domain $\boldsymbol{R}_{nm}$ as follows:

$$
R_{nm}(g) = \prod_{n^{\prime}\in\partial m\backslash n}Q_{mn^{\prime}}(g)\,\,\,(g = 0, 1,\ldots,2^{M}-1).
$$

Then, we apply inverse Fourier transform to $\boldsymbol{R}_{mn}$ and get variable to factor messages in real domain as follows:

$$
\boldsymbol{r}_{nm} = \mathscr{F}^{-1}\left[\boldsymbol{R}_{nm}\right].
$$

Due to discrete Fourier transformation's property, we must do following renormalization:

$$
r_{nm}(g)\Leftarrow r_{nm}(g)/2^{M+1}.
$$

### Variable to factor message update
Variable to factor node message update processes are given by

$$
q_{mn}(g)=f_{n}(g)\prod_{m^{\prime}\in\partial n\backslash m}r_{nm^{\prime}}(g).
$$

After that normalizing $q_{mn}(g)$ with respect $g$.

### Speculating temporal code word
In order to speculate temporal codeword, calculate probability which each symbols equals to $g$. The probability are given by

$$
p_{n}(g)=f_{n}(g)\prod_{m^{\prime}\in\partial n}r_{nm^{\prime}}(g).
$$

Finally, speculated $n-$ th symbol $\hat{x}_{n}$ is given by

$$
\hat{x}_{n}=\operatorname*{argmax}_{g}p_{n}(g).
$$
