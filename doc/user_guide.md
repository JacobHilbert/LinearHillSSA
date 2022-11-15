# User guide


<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=6 orderedList=false} -->

<!-- code_chunk_output -->

- [Model description](#model-description)
- [Model Examples](#model-examples)
  - [Housekeeping gene](#housekeeping-gene)
  - [Self repressor](#self-repressor)
  - [Block matrix notation and generalizable examples](#block-matrix-notation-and-generalizable-examples)
    - [Repressilator](#repressilator)
    - [Loop chains](#loop-chains)
    - [Free chains](#free-chains)
- [Usage examples](#usage-examples)

<!-- /code_chunk_output -->





## Model description

We have

* $S$ chemical species on the state vector $x_S$.

* $N$ reactions, each one with a reaction vector such that $x_S + R_S$ is the new each vector after the reaction, collected as rows in a matrix $R_{NS}$.
    * A full linear dependence term of the propensities on the state vector:
    $$
    A_{NS} x_S +b_N
    $$

* The Hill rate function
$$
\text{Hill}(x,v,\alpha,n) = v \frac{x^n}{x^n + \alpha^n}
$$
    where $n$ is called the hill factor; positive for activation, negative for repression.
    
    * This Hill-type interactions are expressed in matrix notation as
    $$
    \text{Hill}(T_{NS} x_S,v_S, \alpha_S,n_S)
    $$ 
    and the operations are elementwise.
    
So the propensities are calculated as
$$
P(R_{NS}) =  A_{NS} x_S +b_N + \text{Hill}(T_{NS} x_S,v_N, \alpha_N,n_N)
$$

(or if you don't like my matrix notation:)
$$
P(R) = A \vec{x} + \vec{b} + \text{Hill}(T \vec{x},\vec{v},\vec{\alpha},\vec{n})
$$


## Model Examples 

### Housekeeping gene

The equations
$$
\begin{align*}
\frac{\mathrm{d} m\left( t \right)}{\mathrm{d}t} =&\, \gamma_m - k_m m\left( t \right) \\
\frac{\mathrm{d} p\left( t \right)}{\mathrm{d}t} =&\, \gamma_p m\left( t \right) - k_p p\left( t \right)
\end{align*}
$$

correspond to the reactions

$$
\ce{m <=>[\gamma_m][k_m] \emptyset} \qquad \ce{m ->[k_p] m + p} \qquad \ce{p ->[\gamma_p] \emptyset}
$$

Which in matrix notation is:

$$
P\left[\overbrace{\begin{pmatrix}
+1 & 0 \\
-1 & 0 \\
0 & +1 \\
0 & -1
\end{pmatrix}}^{R_{MN}}\right] = 
\overbrace{
\begin{pmatrix}
 0 & 0 \\ 
 \gamma_m & 0 \\ 
 k_m & 0 \\ 
 0 & \gamma_p 
\end{pmatrix}}^{A_{NS}}
\underbrace{\binom{m}{p}}_{X_S} + 
\overbrace{\begin{pmatrix}
 k_m  \\ 0  \\ 0  \\ 0  
\end{pmatrix}}^{b_N} + \text{Hill}\left[ 
\overbrace{\begin{pmatrix}
 0 & 0 \\ 
 0 & 0 \\ 
 0 & 0 \\ 
 0 & 0 
\end{pmatrix}}^{T_{NS}}\underbrace{\binom{m}{p}}_{X_S} ,
\overbrace{\begin{pmatrix}1\\1\\1\\1\\\end{pmatrix}}^{v_N},
\overbrace{\begin{pmatrix}1\\1\\1\\1\\\end{pmatrix}}^{\alpha_N},
\overbrace{\begin{pmatrix}1\\1\\1\\1\\\end{pmatrix}}^{n_N}
\right]
$$

We can recover the differential equations by operating

$$
\frac{\mathrm{d}X}{\mathrm{d} t} = R_{MN}{}^T P(R_{MN})=  \begin{pmatrix}
+1 & -1 & 0 & 0 \\
0 & 0 & +1 & -1
\end{pmatrix} P(R_{MN})
$$

where ${}^T$ is the transpose matrix.

### Self repressor


$$
\begin{align*}
\frac{\mathrm{d} m\left( t \right)}{\mathrm{d}t} =&\, \gamma_m - k_m m\left( t \right)  + v \frac{p^{-n}}{p^{-n}+\alpha^{-n}}\\
\frac{\mathrm{d} p\left( t \right)}{\mathrm{d}t} =&\, \gamma_p m\left( t \right) - k_p p\left( t \right)
\end{align*}
$$

$$
\ce{\emptyset ->[Hillr(p)] m}\qquad \ce{m <=>[\gamma_m][k_m] \emptyset} \qquad \ce{m ->[k_p] m + p} \qquad \ce{p ->[\gamma_p] \emptyset}
$$

$$
P\left[\overbrace{\begin{pmatrix}
+1 & 0 \\
-1 & 0 \\
0 & +1 \\
0 & -1
\end{pmatrix}}^{R_{MN}}\right] = 
\overbrace{
\begin{pmatrix}
 0 & 0 \\ 
 \gamma_m & 0 \\ 
 k_m & 0 \\ 
 0 & \gamma_p 
\end{pmatrix}}^{A_{NS}}
\underbrace{\binom{m}{p}}_{X_S} + 
\overbrace{\begin{pmatrix}
 k_m  \\ 0  \\ 0  \\ 0  
\end{pmatrix}}^{b_N} + \text{Hill}\left[ 
\overbrace{\begin{pmatrix}
 0 & 1 \\ 
 0 & 0 \\ 
 0 & 0 \\ 
 0 & 0 
\end{pmatrix}}^{T_{NS}}\underbrace{\binom{m}{p}}_{X_S} ,
\overbrace{\begin{pmatrix}v\\1\\1\\1\\\end{pmatrix}}^{v_N},
\overbrace{\begin{pmatrix}\alpha\\1\\1\\1\\\end{pmatrix}}^{\alpha_N},
\overbrace{\begin{pmatrix}-n\\1\\1\\1\\\end{pmatrix}}^{n_N}
\right]
$$

### Block matrix notation and generalizable examples

using the language of block matrices and vectors, we adopt the notation
$$
\alpha_{33} =  I_3 \alpha = \begin{pmatrix} 1&0&0\\ 0&1&0\\ 0&0&1 \end{pmatrix} \alpha
\qquad \alpha_3 =  \begin{pmatrix}1\\1\\1\end{pmatrix} \alpha
$$
with $I_3 = 1_{33}$ being the 3x3 identity matrix, and
$$
\begin{pmatrix}
a & b \\ c & d
\end{pmatrix}_{33} = \begin{pmatrix}
a_{33} & b_{33} \\ c_{33} & d_{33}
\end{pmatrix} = \left(
\begin{array}{ccc|ccc}
 a & 0 & 0 & b & 0 & 0 \\
 0 & a & 0 & 0 & b & 0 \\
 0 & 0 & a & 0 & 0 & b \\\hline
 c & 0 & 0 & d & 0 & 0 \\
 0 & c & 0 & 0 & d & 0 \\
 0 & 0 & c & 0 & 0 & d \\
\end{array}
\right)
\qquad 
\begin{pmatrix}a\\b\end{pmatrix}_3 = \begin{pmatrix}a_3\\b_3\end{pmatrix} = \begin{pmatrix}a\\a\\a\\b\\b\\b\end{pmatrix}
$$

It is to note that if, for instance $b$ where to be a matrix itself, $b_{33}$ is defined as $b$, and the same for $b_3\equiv b$ if $b$ is a vector.

#### Repressilator 

With this language, the repressilator can be expressed as

$$
P\left[\overbrace{\begin{pmatrix}
+1 & 0 \\
-1 & 0 \\
0 & +1 \\
0 & -1
\end{pmatrix}_{33}}^{R_{MN}}\right] = 
\overbrace{
\begin{pmatrix}
 0 & 0 \\ 
 \gamma_m & 0 \\ 
 k_m & 0 \\ 
 0 & \gamma_p 
\end{pmatrix}_{33}}^{A_{NS}}
\;X + 
\overbrace{\begin{pmatrix}
 k_m  \\ 0  \\ 0  \\ 0  
\end{pmatrix}_3}^{b_N} + \text{Hill}\left[ 
\overbrace{\begin{pmatrix}
 0 & P_{312} \\ 
 0 & 0 \\ 
 0 & 0 \\ 
 0 & 0 
\end{pmatrix}_{33}}^{T_{NS}}\;X ,
\overbrace{\begin{pmatrix}v\\1\\1\\1\\\end{pmatrix}_3}^{v_N},
\overbrace{\begin{pmatrix}\alpha\\1\\1\\1\\\end{pmatrix}_3}^{\alpha_N},
\overbrace{\begin{pmatrix}-n\\1\\1\\1\\\end{pmatrix}_3}^{n_N}
\right]
$$

where $X=\{m_1,m_2,m_3,p_1,p_2,p_3\}$. The only qualitative difference with the self-repressor is the row-permutation matrix

$$
P_{312} = \begin{pmatrix} 0&0&1\\ 1&0&0\\ 0&1&0 \end{pmatrix}
$$

The expression
$$
\frac{\mathrm{d}X}{\mathrm{d} t} = \begin{pmatrix}
+1 & -1 & 0 & 0 \\
0 & 0 & +1 & -1
\end{pmatrix}_{33} P(R_{MN})
$$

gives the equations
$$
\begin{array}{ccccccc}
m_1{}' & = & k_m & - & m_1\; \gamma _m &+& \text{Hill}\left(p_3,v,\alpha ,-n\right) \\ 
m_2{}' & = & k_m & - & m_2\; \gamma _m &+& \text{Hill}\left(p_1,v,\alpha ,-n\right) \\ 
m_3{}' & = & k_m & - & m_3\; \gamma _m &+& \text{Hill}\left(p_2,v,\alpha ,-n\right) \\ 
p_1{}' & = & m_1 \;k_p &-&p_1 \;\gamma _p & \\ 
p_2{}' & = & m_2 \;k_p &-&p_2 \;\gamma _p & \\ 
p_3{}' & = & m_3 \;k_p &-&p_3 \;\gamma _p & \\ 
\end{array}
$$

Notice that each of the the $\gamma_{m,p},k_{m,p}$ could have been a 3-vector, so this is could have easily been the general repressilator.

#### Loop chains

We can generalize this to a $\ell$-th order cyclic repressilator by replacing all the 3s by $\ell$ and the permutation matrix by 

$$
P_{\text{RotateRight}(\text{Range}(\ell))} = \begin{pmatrix} 
0 & 0 & 0 & 0 & \cdots & 1 \\
1 & 0 & 0 & 0 & \cdots & 0 \\
0 & 1 & 0 & 0 & \cdots & 0 \\
0 & 0 & 1 & 0 & \cdots & 0 \\
\vdots & \ddots & \ddots & \ddots & \ddots & 0 \\
 0 & 0& 0& 0 & 1 & 0
\end{pmatrix}
$$

#### Free chains

If we replace the permutation matrix $P_{312}$ with a -1 diagonal matrix
$$
P_{s-1} = \begin{pmatrix} 
0 & 0 & 0 & 0 & \cdots & 0 \\
1 & 0 & 0 & 0 & \cdots & 0 \\
0 & 1 & 0 & 0 & \cdots & 0 \\
0 & 0 & 1 & 0 & \cdots & 0 \\
\vdots & \ddots & \ddots & \ddots & \ddots & 0 \\
 0 & 0& 0& 0 & 1 & 0
\end{pmatrix}
$$

one obtains a free chain, that does not loop onto itself like the $\ell$-th order chain from before.

***

## Usage examples

(Work in progress here, but head to the [repressilator example](example.py) or the [benchmark](benchmark.py))











