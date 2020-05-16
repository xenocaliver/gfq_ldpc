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
Let $\boldsymbol{c}$ a codeword(column vector), $H$ a parity check matrix. If $\boldsymbol{c}$ is a codeword, 
$$ H\boldsymbol{c}=\boldsymbol{0} $$
holds. And we assume we can split $\boldsymbol{c}$ into two parts $\boldsymbol{p}$ and $\boldsymbol{s}$. $\boldsymbol{p}$ represents a parity check symbol part and $\boldsymbol{s}$ represents input symbol part.
And we devide $H$ into two part i.e.
$$ H=\left(A\vert B\right). $$
$A$'s number of columns coinides with $\boldsymbol{p}$'s dimension. $B$'s number of columns coincides with $\boldsymbol{s}$'s dimension. So,
$$ A\boldsymbol{p}+B\boldsymbol{s}=\boldsymbol{0} $$
holds. Therefore one can calculate $\boldsymbol{p}$ by means of following expression:
$$ \boldsymbol{p}=-A^{-1}B\boldsymbol{s}. $$
Next, one can calculate $A^{-1}$ by means of $LU$ decompsition. However, one can apply $LU$ decomposition only if $A$ is full rank. However $A$ is not always full rank. So, `make-gen` permuates $H$'s columns and make $A$ be full rank. A column permuation corresponds to a matrix. Let $Q_{i}(i=1,\ldots, m)$ be a matrix corresponding to a column exchange. So we can write
$$ \boldsymbol{c}=\left( \begin{array}{c} -A^{-1}B\\ I \end{array} \right)Q_{1}\cdots Q_{m}\boldsymbol{s} $$
where $I$ is an identity matrix and we used $Q_{i}^{-1}=Q_{i}$. Therefore we obtain generating matrix $G$ as follows:
$$ G=\left( \begin{array}{c} -A^{-1}B\\ I \end{array} \right)Q_{1}\cdots Q_{m}. $$

## Transmitting a symbol over AWGN channel
We assumed that channel is AWGN channel with standard deviation $\sigma$ and assumed
finite field be $\mathbb{F}_{q}, q=2^{m}$. According to [Davey](https://ieeexplore.ieee.org/document/706440), each bit of a symbol is modulated to BPSK and bitwise prior probability for bit $k$ is given by
$$ p_{k}=\frac{1}{1+\exp[2s_{k}\vert y_{k}\vert/\sigma^{2}]} $$
where $s_{k}$ is BPSK modulated value of the bit $k$ and $y_{k}$ is a observed value of bit $k$. Therefore prior probability according to symbol $g\in\mathbb{F}_{q}$ is given by
$$f(g) = \prod_{k=1}^{m}\frac{1}{1+\exp[2s_{k}\vert y_{k}\vert/\sigma^{2}]}.$$

## Decoding
We employ probability region Sum-Product algorithm for decoding according to [Davey](https://ieeexplore.ieee.org/document/706440).
Let $q_{mn}(g)$ be a veriable node $n$ to factor node message $m$ where g is a finite field element's value. Let $r_{nm}$ be a factor node $m$ to variable node $n$ message.

### Initialization
For all variable nodes, initialize messages as follows:
$$
q_{mn}(g)=f_{n}(g).
$$

### Factor to variable node message update
Factor to variable node message update processes are given by
$$
r_{nm}(g)=\sum_{x_{i_{1}}}\cdots\sum_{x_{i_{d_{c}-1}}}\boldsymbol{1}[\sum_{k=1}^{d_{c}}h_{mi_{k}}x_{i_{k}}=0]\prod_{i_{k}\in∂ m\backslash n}q_{mn^{\prime}}(x_{i_{k}})
$$
where $\boldsymbol{1}$ denotes an indicator function i.e.
$$
\boldsymbol{1}[A]=
\begin{cases}
1& (A\quad\mathrm{is\quad true})\\
0& (A\quad\mathrm{is\quad false}).
\end{cases}
$$

### Variable to factor message update
Variable to factor node message updte processes are given by

$$
q_{mn}(g)=f_{n}(g)\prod_{m^{\prime}\in∂ n\backslash m}r_{nm^{\prime}}(g).
$$

After that normalizing $q_{mn}(g)$ with respect $g$.

### Speculating tempral code word
In order to speculate temporal codeword, calculate probability which each symbols equals to $g$. The probability are given by
$$
p_{n}(g)=f_{n}(g)\prod_{m^{\prime}\in∂ n}r_{nm}(g)
$$
Finally, speculated $n-$th symbol $\hat{g}_{n}$ is given by
$$
\hat{g}_{n}=\operatorname{argmax}_{g}p_{n}(g).
$$
