---
title: "XorShift Jump 101, Part 1: Matrix Multiplication"
description: "Jump forward/backward procedures for XorShift RNG explained step by step."
---
# {{page.title}}

{{page.description}}

Published: 2026-06-07.

Originally published: [CodeProject](https://web.archive.org/web/20250819132602/https://www.codeproject.com/Articles/5264513/XorShift-Jump-101-Part-1-Matrix-Multiplication), 2020-04-12.

\[[Source Code](https://github.com/sergeykorop/xorshift-jump/)\]

## Introduction

This article is a sequel to
[`System.Random` and Infinite Monkey Theorem](/monkey-typewriter/),
where we explored the internal structure of the standard random number
generator (RNG) from the .NET Framework. Using the theory of linear
RNGs, we learned how to navigate the sequence of random numbers
jumping $N$ steps forward or backward with logarithmic
complexity. This time, we are going to develop similar functionality
for another family of linear RNGs, known as
[xorshift](https://en.wikipedia.org/wiki/Xorshift).

We will start from building transition matrices for this RNG. Unlike
[lagged Fibonacci generators](https://en.wikipedia.org/wiki/Lagged_Fibonacci_generator),
these matrices have more complex structure so we would prefer to
calculate them automatically. This calculation is performed by a
program written in a powerful but not very popular language,
Haskell. Giving a test drive to Haskell was one of my goals for this
research and it appeared to be worth the efforts. After that, we will
consider another method of jumping based on polynomial arithmetics
\[[2](#ref.2)\] which has been mentioned in the previous article
but not tried yet. Finally, we will compare these two approaches.

Putting all the material in a single article would made it a bit too
long to read at once, so I've decided to split it into two parts. Part&nbsp;1
is dedicated to matrix-based approach, while part&nbsp;2 will consider
the polynomial-based algorithm and performance testing. For the same
reason, I will focus on math background and general ideas rather than
software implementation. The exception is made for code generator in
Haskell: assuming this language is not widely used, I stepped through
each line of the code explaining it in more details.

Why do we need one more RNG algorithm? Unfortunately, `System.Random`
is not perfect. First of all, statistical testing exposed some
problems immanent to the whole family of lagged Fibonacci generators
\[[3](#ref.3), §3.2.2\], \[[4](#ref.4), ch. 7\]. Besides that, a
critical bug has been found in the .NET implementation of this
algorithm: lag values specified in the code don't guarantee the
longest period possible. In fact, I don't know any theoretical
estimate of the period length for these lags. This issue has been
discussed many times by the open source .NET community
([\#5974](https://github.com/dotnet/coreclr/issues/5974),
[\#12746](https://github.com/dotnet/corefx/issues/12746),
[\#23298](https://github.com/dotnet/corefx/issues/23298))
but the code remains as is for the sake of backward
compatibility. (See also an excellent
[research](https://gist.github.com/fuglede/772402ecc3997ada82a03ce65361781b)
made by Søren Fuglede Jørgensen.)

Taking these facts into account, I wouldn't recommend using
`System.Random` if you need not only random numbers but also some kind
of proven estimate of their properties.

Now it would be natural to ask a question: why have we explored
`System.Random` being aware of its bad sides? The answer is simple:
while being bad for the real-world use, this RNG is good to explain
the theory. Indeed, lagged Fibonacci algorithm uses only simple
integer arithmetics and its transition matrices are clear enough to be
written directly without complex calculations. Also, the state vector
of this RNG is sufficiently large to keep a text message, the trick
which brought my previous article to life.

## Background

It is assumed that the reader is familiar with linear RNG theory
basics explained in
[`System.Random` and Infinite Monkey Theorem](/monkey-typewriter/).
Unlike that article, the demo code for this one is written in C++ 11
so the reader should be fluent in this dialect at the intermediate
level.

Polynomial-based approach requires some additional knowledge in the
basics of algebra, such as the notion of polynomial, its properties
and operations involving polynomials. Software implementation is based
on [NTL](http://www.shoup.net/ntl/), a library of number theory
methods developed by Victor Shoup. This library, in turn, needs
[GNU GMP](http://gmplib.org/) and
[GF2X](http://gforge.inria.fr/projects/gf2x/). Versions used for this
article are NTL 11.0.0, GMP 6.1.2, and GF2X 1.2, respectively.

While NTL documentation mentions that it is possible to build this
library with Visual C++, the recommended environment should be
Unix-style, such as G++ (with accompanying toolchain) under Linux or
Cygwin (I used this one).

Some C++ code for this article has been generated automatically by a
program written in Haskell. Source code archive includes ready-to-use
output of this generator so the reader doesn't necessarily have to run
it. Otherwise, some working Haskell system, e.g.,
[Haskell Platform](https://www.haskell.org/platform/windows.html),
would be needed. Readers not familiar with Haskell would also need some
beginner-level guide. I used
[Yet Another Haskell Tutorial](http://users.umiacs.umd.edu/~hal/docs/daume02yaht.pdf)
by Hal Daumé III.

[Google Test](https://github.com/google/googletest) 1.8.1 is used for unit testing infrastructure.

## First Steps

Let's begin from implementing basic RNG functionality for xorshift
algorithm.

First of all, we should choose some particular algorithm from this
family to be used in our studies. Let's implement `xor128` defined by
George Marsaglia in \[[5](#ref.5), p. 5\].

All code related to basic RNG functions would be placed to file
_xs.cc_ with forward declaration of defined entities in header file
_xs.h_. The important design trade-off has been made in this code
which is popular in research programming: we won't follow the
encapsulation principles of the OOP for the sake of maximum
flexibility and expose all internal state of the RNG and other objects
to public scope. Also, I would prefer using &ldquo;bare&rdquo; data structures
processed by non-member functions. Again, this is a research
programming trade-off: when you do a research, you often need to
define some alternative actions on data, and it is convenient to do
right in place where these actions are used instead of adding more and
more class members to some common class definition or bothering with
inheritance.

Marsaglia defines the internal state of `xor128` as four separate
32-bit variables: `x`, `y`, `z`, and `w`, respectively. Since we are
going to apply our matrix-based approach, this state should be treated
as single 128-bit vector. It would be convenient to pack them into
4-element vector of 32-bit unsigned integers to perform matrix
operations uniformly but sometimes we will use `x`, `y`, `z`, and `w`
as aliases for corresponding vector items when they are treated
differently.

```cpp
using state_t = std::array<uint32_t, 4>;
```

This state vector can be initialized by arbitrary value except all
zeroes but originating paper defines some &ldquo;canonical&rdquo; initial seed
which we will use as well.

```cpp
void init(state_t &s)
{
  s[0] = 123456789;
  s[1] = 362436069;
  s[2] = 521288629;
  s[3] = 88675123;
}
```

The next step would be implementing functions to produce random
numbers from this state. This process can be divided into two
sub-problems: we should mutate RNG state to its next value and then
produce the next item of output using this new state. Common practice
is putting these two into one piece of code (Marsaglia's paper follows
it, too) but in our research, we need to be able to separate these
actions. Indeed, we have to be able to either do the next step as
usual or compute the next state using some jump algorithm and then
resume producing random numbers from that state like it has been
obtained during previous access to the generator. So, we will define
two separate functions:

```cpp
void step_forward(state_t &s)
{
  uint32_t &x = s[0], &y = s[1], &z = s[2], &w = s[3];
 
  uint32_t t = x ^ (x << 11);

  x = y; y = z; z = w;
  w = w ^ (w >> 19) ^ (t ^ (t >> 8));
}

uint32_t gen_u32(state_t &s)
{
  step_forward(s);

  return s[3];
}
```

## Matrix-Based Approach

After defining `step_forward()`, it's good time now to discuss why
xorshift RNGs belong to linear family. Section
[Checking for linearity](#checking-for-linearity) involves some amount
of math so readers more interested in programming can skip it till the
final conclusion. Next section,
[Calculating transition matrix](#calculating-transition-matrix),
provides step by step explanation on how to
generate needed transition matrix programmatically.

### Checking for Linearity

Our goal will be to show that single step of xorshift method can be
expressed as $S^{(i+1)} = TS^{(i)}$, where $T$ is a constant matrix
known as _transition matrix_ and $S^{(i)}$ is a state vector at step
$i$. For xorshift both $S$ and $T$ are defined over $\mathbb{F}_2$:
$S \in \mathbb{F}_2^{N}$, and $T \in \mathbb{F}_2^{N\times N}$, where $N$
is a size of state vector ($N=128$ for the rest of this article). In
programmers' parlance, it means that items of these vectors or
matrices are single bits and allowed arithmetical operations are
multiplication (bitwise &ldquo;and&rdquo;) and addition modulo $2$
(bitwise &ldquo;xor&rdquo;).

To simplify our formulas, let's omit the step index and denote
$S'=S^{(i+1)}$ and $S=S^{(i)}$. Our transition formula then would be
$S' = TS$. Also, let's use lowercase letters with indices for vector
and matrix items, such as $s_i$ or $t_{i,j}$, $0\le i < N$,
$0\le j < N$.

Given that $S'=TS$, every component of $S'$ is

$$
\label{eq:component} s'_i = \bigoplus_{0\le j < N} t_{i,j}s_j, \qquad 0 \le i < N. \tag{1}
$$

Let's begin from parts of updated state vector named `x`, `y`, `z`.
They cover state vector bits from $s_0$ till $s_{95}$ and are
obtained by copying other 32-bit chunks as is:

$$
s'_i = s_{i+32}, \qquad 0 \le i < 96.
$$

Therefore, $75\%$ of our transition matrix rows would be like this:

$$
t_{i,j} = \begin{cases}
 1, & j = i + 32,\\
 0, &\mbox{otherwise}\\
\end{cases}, \qquad 0 \le i < 96, 0 \le j < 128,
$$

satisfying the requirements.

On the other hand, chunk `w` is calculated in a more complicated
manner: we can see a composition of bitwise &ldquo;xor&rdquo;
operations with bitwise shifts. Writing exact formulas for each
component would be time-consuming so we would apply induction:

1. The innermost elements of the formulas are state vector components,
so they match the structure of ($\ref{eq:component}$).

2. Given subexpressions matching ($\ref{eq:component}$), we can apply
addition modulo 2 and get the sum satisfying the same condition (this
follows from properties of addition).

3. Given a vector of subexpressions matching ($\ref{eq:component}$),
we can see that shift operation either replaces subexpression for one
component by subexpression for another, keeping property
($\ref{eq:component}$) intact, or replaces it with zero which is also
conformant to that condition.

Putting it all together, we can conclude that transformations defined
in `step_forward()` are indeed suitable for matrix-based jump
algorithms described in my previous article.

### Calculating Transition Matrix

As mentioned before, it is possible to calculate transition matrix $T$
from code implementing xorshift RNG. To do that, we need to write
formulas for each operation in the code componentwise, i.e., for every
bit of state vector $S$. This is tedious and error-prone. Even worse,
since xorshift is a family of RNGs, we will have to calculate matrix
$T$ for every kind of xorshift from the very beginning. The solution
is simple: let's write a program to do that for us automatically.

This program should take some formulas describing given kind of
xorshift and transform those formulas using rules of arithmetics into
the form ($\ref{eq:component}$). Doing that for every of $N$ bits of
the state vector $S$, we will get $N\times N$ items of the matrix $T$.

Symbolic manipulations with formulas require appropriate language to
implement them. The main data structure to process would be a
tree-like structure whose nodes represent each constant, variable or
operation and their interaction in formulas. Functional language
[Haskell](https://haskell.org/) seems be an appropriate tool due to
its first-class support for recursive data structures, algebraic data
types and pattern matching.

Also, Haskell provides a powerful feature to help with making
experiments and rapid prototyping, known as
[REPL](https://en.wikipedia.org/wiki/Read%E2%80%93eval%E2%80%93print_loop) (Read-eval-print-loop).
Using REPL, we can implement our matrix generator interactively,
specifying type and function definitions one by one and use them
immediately from the interactive Haskell shell to see how they work.

Source code archive provided with the article contains file
_xsTgen.hs_. I'm going to explain how it works step by step. To do
that, we will load this file into the Haskell interpreter, named
_ghci_ if you are using the recommended Haskell Platform. After that,
we will enter some Haskell expressions to the REPL and examine their
output. Note that development process was very similar: I wrote
function definitions one by one and tested them immediately. (This
process is much easier if you use Haskell-aware editor, e.g. Emacs.)
Such feedback is very useful when you do exploratory programming.

To start the interpreter, open command-line prompt window, change
current directory to be one containing our source file, _xsTgen.hs_
and run _ghci_. If everything went fine, you should see this kind of
start-up message:

```text
GHCi, version 8.6.3: http://www.haskell.org/ghc/  :? for help
Prelude>
```

Load source code using command `:l`:

```text
Prelude> :l xsT
[1 of 1] Compiling Main             ( xsTgen.hs, interpreted )
Ok, one module loaded.
*Main>
```

Windows users may also start _WinGHCi_ as GUI application and load
source file via menu.

Let's begin our tour from data type representing formulas of xorshift
transformations:

```haskell
data Expr =
  Zero
  | Var String Integer
  | Xor [Expr]
  deriving (Eq, Show)
```

Here, we have a [tagged union](https://en.wikipedia.org/wiki/Tagged_union)
which may store either constant zero or variable referenced by name
with index or a set of arbitrary expressions xor'ed together. The
clause `deriving` tells Haskell to create some helper code for our
data type automatically.

We can see now the information about defined data type:

```text
*Main> :i Expr
data Expr = Zero | Var String Integer | Xor [Expr]
        -- Defined at xsTgen.hs:1:1
instance [safe] Show Expr -- Defined at xsTgen.hs:5:17
instance [safe] Eq Expr -- Defined at xsTgen.hs:5:13
```

and type some expressions to see how it works:

```text
*Main> Zero
Zero
*Main> Var "x" 1
Var "x" 1
*Main> Xor [Var "x" 1, Var "y" 2]
Xor [Var "x" 1,Var "y" 2]
*Main> Xor [Var "x" 1, Xor[Zero, Var "y" 2]]
Xor [Var "x" 1,Xor [Zero,Var "y" 2]]
```

These expressions define internal representation for such formulas as:
$0$, $x_1$, $x_1\oplus y_2$, $x_1\oplus (0\oplus y_2)$,
respectively. Note that our `Expr` is not enough to represent the
complete $\mathbb{F}_2$ arithmetics: we don't have constant $1$ or
multiplication to write things like $y_1\oplus 1$ or
$x_1\cdot y_2$. This is a bit restrictive but we don't need such ability to
solve our current problem.

Now it's time to teach our code to process vectors of expressions so
we won't have to write formulas for each bit of RNG state. If we look
at the _xsTgen.hs_, data type `Expr` is followed by this set of
definitions:

```haskell
bus name cntr = map (Var name) [0 .. cntr - 1]

zeroes = Zero : zeroes

xor x y = zipWith (\x y -> Xor [x, y]) x y

crop x y = zipWith (\ x y -> x) x y

shl x n = crop (take n zeroes ++ x) x

shr x n = crop (drop n x ++ zeroes) x
```

These equations define some simple functions. Let's see how they
work. Function `bus` is used to create a vector of indexed variables
given their name and size:

```text
*Main> bus "x" 3
[Var "x" 0,Var "x" 1,Var "x" 2]
*Main> bus "y" 4
[Var "y" 0,Var "y" 1,Var "y" 2,Var "y" 3]
```

Helper function `zeroes` produces an infinite sequence of zero
constants. Therefore, we can't evaluate it as is but we can consume
some specified number of zeroes when needed:

```text
*Main> take 5 zeroes
[Zero,Zero,Zero,Zero,Zero]
*Main> take 7 zeroes
[Zero,Zero,Zero,Zero,Zero,Zero,Zero]
```

Function `Xor` takes two expression vectors and produces a new one
xoring them component-wise:

```text
*Main> xor (bus "x" 3) (bus "y" 3)
[Xor [Var "x" 0,Var "y" 0],Xor [Var "x" 1,Var "y" 1],Xor [Var "x" 2,Var "y" 2]]
```

One more helper function `crop` takes two expression vectors and
returns the components of the first vector taken up to the length of
the second one:

```text
*Main> crop (bus "x" 4) (bus "y" 3)
[Var "x" 0,Var "x" 1,Var "x" 2]
```

Using this helper, we are now able to define two shift functions
operating on expression vectors moving their components:

```text
*Main> shl (bus "x" 4) 2
[Zero,Zero,Var "x" 0,Var "x" 1]
*Main> shr (bus "x" 4) 2
[Var "x" 2,Var "x" 3,Zero,Zero]
```

Note that function names seem counter-intuitive if we look at their
results as expression vectors per se but they would be more clear if
we interpret vector components as bits in the word so shifting them
left indeed means shifting towards the most-significant bit.

Using these primitives, we can define xorshift-like transformations
with formulas quite close to the original C or C++ code:

```text
*Main> let x = bus "x" 3
*Main> let y = bus "y" 3
*Main> xor x (shl y 2) ++ y
[Xor [Var "x" 0,Zero],Xor [Var "x" 1,Zero],Xor [Var "x" 2,Var "y" 0],
Var "y" 0,Var "y" 1,Var "y" 2]
```

We can also use one more feature of Haskell, the ability to treat
functions like binary operations enclosing their names in backquotes:

```text
*Main> x `xor` (y `shl` 2) ++ y
[Xor [Var "x" 0,Zero],Xor [Var "x" 1,Zero],Xor [Var "x" 2,Var "y" 0],
Var "y" 0,Var "y" 1,Var "y" 2]
```

In fact, we created a kind of
mini-[DSL](https://en.wikipedia.org/wiki/Domain-specific_language)
for our problem.

Now it's time to learn how to transform our expressions. The ultimate
goal would be a linear form ($\ref{eq:component}$) which is directly
convertible to transition matrix.

Looking at the definition of type `Expr`, we can see that its values
are trees where leaves can be either zeroes or variable components and
all inner nodes are `Xor` operations applied to subtrees. In math, we
express this with formulas using nested subexpressions:
$x_1 \oplus (x_2 \oplus (y_3 \oplus \ldots))$
or $x_1 \oplus ((x_2 \oplus y_3) \oplus \ldots))$.

Taking into account properties of addition modulo 2, we can convert
these expressions to the linear form which would be the same for both
of them: $x_1 \oplus x_2 \oplus y_3 \oplus \ldots$ Function `flatten`
works exactly this way:

```haskell
flatten::Expr -> [Expr]
flatten (Xor (x:xs)) = flatten x ++ (concatMap flatten xs)
flatten (Xor [])   = []
flatten x          = [x]
```

It takes an expression tree and returns a sequence of primitives
(constants or variables) which are implicitly xor'ed together. To do
that, we traverse the expression tree recursively and for each `Xor`
operation node we flatten each subtree and concatenate them
together. Other kinds of nodes are collected into this list as is.

Let's apply this function to expression trees described above:

```text
*Main> flatten (Xor [Var "x" 1, Xor [Var "x" 2, Xor [Var "y" 3, Zero]]])
[Var "x" 1,Var "x" 2,Var "y" 3,Zero]
*Main> flatten (Xor [Var "x" 1, Xor [Xor[Var "x" 2, Var "y" 3], Zero]])
[Var "x" 1,Var "x" 2,Var "y" 3,Zero]
```

(Since our expression trees are not suited to describe incomplete
formulas, I replaced the ellipsis with zero constant to make them
well-formed.) Also note that we can get duplicates after flattening:

```text
*Main> flatten (Xor [Var "x" 1, Xor [Xor[Var "x" 2, Var "x" 1], Zero]])
[Var "x" 1,Var "x" 2,Var "x" 1,Zero]
```

We've got a significant progress on the way towards our final goal! We
are able now to define xorshift transformations in a concise form and
convert formulas to linear sequence of involved variables. Our next
step would be to convert those lists of variables to vectors of
coefficients which, in turn, will become rows of transition matrix.

Let's define one more function:

```haskell
matrix e vars = zipWith (\e v -> (v, makeRow e)) e vars where
  makeRow e = map (\v -> length (filter (== v) e) `mod` 2) vars
```

Parameter `e` is a list of expression lists `[[Expr]]`. Each list item
corresponds to single bit of state vector and this item is a list of
`Var` or `Zero` items produced by `flatten`. Second parameter, `vars`,
is a list of `Var` where $i$-th item denotes variable corresponding to
the $i$-th item of the state vector. We have already introduced `x`
and `y` to be our building blocks and we can evaluate them in REPL to
recall their values. Now it's time to define a formula we are going to
process:

```text
*Main> x
[Var "x" 0,Var "x" 1,Var "x" 2]
*Main> y
[Var "y" 0,Var "y" 1,Var "y" 2]
*Main> let e = x `xor` (y `shl` 2) ++ y
*Main> let v = x ++ y
*Main> e
[Xor [Var "x" 0,Zero],Xor [Var "x" 1,Zero],Xor [Var "x" 2,Var "y" 0],Var "y" 0,
Var "y" 1,Var "y" 2]
*Main> v
[Var "x" 0,Var "x" 1,Var "x" 2,Var "y" 0,Var "y" 1,Var "y" 2]
```

Here, we have two vector variables, `x` and `y`. These three-bit
vectors are combined to state vector
$S=\langle x_0, x_1, x_2, y_0, y_1, y_2\rangle$.
Expression `e` is a xorshift-like transformation over these
vectors. To get the first argument of `matrix`, we should apply
`flatten` to expression `e` componentwise. Variable `v` is just a
collection of all components from `x` and `y`, listed in one sequence
as needed by the second argument of `matrix`.

```text
*Main> matrix (map flatten e) v
[(Var "x" 0,[1,0,0,0,0,0]),(Var "x" 1,[0,1,0,0,0,0]),(Var "x" 2,[0,0,1,1,0,0]),
(Var "y" 0,[0,0,0,1,0,0]),(Var "y" 1,[0,0,0,0,1,0]),(Var "y" 2,[0,0,0,0,0,1])]
```

Great! The output looks exactly as specified by
($\ref{eq:component}$): for each bit of the state vector we have a
vector of coefficients.

To see how `matrix` works, let's choose some particular state vector
item to examine. Let it be the $x_2$ since its coefficient vector has
2 non-zero items while others have only one. We should prepare it
first:

```text
*Main> let e2 =  flatten (e!!2)
*Main> e2
[Var "x" 2,Var "y" 0]
```

Since we are going to consider the single state bit, calculating
`matrix` would effectively be the same as calculating `makeRow` with
appropriate arguments, which in turn may be replaced by its
definition:

```text
*Main> map (\v' -> length (filter (== v') e2) `mod` 2) v
[0,0,1,1,0,0]
```

(Note that I had to replace function argument name `v` to `v'` to
avoid conflict with the top-level name used to keep list of involved
variables.)

Now let's see what will happen if we apply the innermost expression:

```text
*Main> v
[Var "x" 0,Var "x" 1,Var "x" 2,Var "y" 0,Var "y" 1,Var "y" 2]
*Main> map (\v' -> filter (== v') e2) v
[[],[],[Var "x" 2],[Var "y" 0],[],[]]
```

As you can see, each variable mentioned in `v` maps to a list
containing itself and the rest of variables is mapped to empty list.

Adding one more level, we wrap the innermost `filter` with `length`:

```text
*Main> map (\v' -> length (filter (== v') e2)) v
[0,0,1,1,0,0]
```

Looks like what we needed! But we have also taken a modulo $2$ as the
final stroke in `matrix` definition. We have to do that to cover the
case when some variable is referred more than once in the same
expression. Variables with even number of occurrences should disappear
due to mutual cancellation of these terms and taking modulo $2$ does
exactly that: all even counts become zeroes and all odd counts become
ones.

**A word of caution:** Looking at our code, we can see that each
expression list out of $N$ lists (where $N$ is a state vector length)
is matched against the list of variables which is also of length
$N$. This gives us the resulting complexity of $O(N^2)$. We could
probably improve it using sorted expression lists but since this code
is not intended to be run often, let's trade its efficiency for
simplicity.

The problem of getting transition matrix for given xorshift
transformation has been solved. Some technical issues still need their
resolution, however. First of all, we need to pack our bit matrix to
integers for future use in C++ code. This goal is achieved with the
following function:

```haskell
packMatrix m vars = map (\(v, b) -> (v, tail (packRow b vars "" 0))) m where
  packRow xb@(b:bs) xv@((Var n i):vs) vname acc | n == vname = packRow bs vs vname (acc + b*2^i)
                                                | otherwise  = acc : packRow xb xv n 0
  packRow [] _ _ acc = [acc]
  packRow _ [] _ acc = [acc]
```

Let's test it:

```text
*Main> packMatrix (matrix (map flatten e) v) v
[(Var "x" 0,[1,0]),(Var "x" 1,[2,0]),(Var "x" 2,[4,1]),(Var "y" 0,[0,1]),
(Var "y" 1,[0,2]),(Var "y" 2,[0,4])]
```

Function `packMatrix` takes each row of the transition matrix
constructed with function `matrix` and groups bits corresponding to
the same family of variables. That is, if we have two families of
variables in our expression, $x_i$ and $y_i$, we will have two integer
words per matrix row. Within the same word, bit weight corresponds to
variable index, e.g., bit row item for $x_i$ takes the $i$-th bit
within word representing all $x$ components. For example, in matrix
row for `Var "x" 2` we see `[4, 1]`. First word is for $x$, and its
the only non-zero bit has position $2$ (counting from zero). Indeed,
the expression for $x_2$ includes $x_2$ itself. Doing the same with
the second word, we can see that another term would be $y_0$.

Next function will take a packed transition matrix and generate C++
code for array initializer:

```haskell
matrixCode pm = header ++ (concat $ intersperse ",\n" 
(map (\(v, r) -> "  {" ++ hexRow r ++ "}") pm)) ++ "\n};\n" where
  header = printf "#include <cstdint>\n\nuint32_t T[%d][%d] = \
           \{\n" (length pm) (length $ snd $ head pm)
  hexRow r = concat $ intersperse ", " (map (printf "0x%08XU") r)
```

Its result looks like this:

```text
*Main> matrixCode $ packMatrix (matrix (map flatten e) v) v
"#include <cstdint>\n\nuint32_t T[6][2] = {\n  {0x00000001U, 0x00000000U},\n  {0
x00000002U, 0x00000000U},\n  {0x00000004U, 0x00000001U},\n  {0x00000000U, 0x0000
0001U},\n  {0x00000000U, 0x00000002U},\n  {0x00000000U, 0x00000004U}\n};\n"
```

We are almost done! Let's define `main` function:

```haskell
main =
  do
    putStr (matrixCode $ packMatrix (matrix (map flatten xs128) vars) vars)
  where x = bus "x" 32 
        y = bus "y" 32
        z = bus "z" 32
        w = bus "w" 32
        vars = (x++y++z++w)
        t = x `xor` (x `shl` 11)
        xs128 = y ++ z ++ w ++ (w `xor` (w `shr` 19) `xor` t `xor` (t `shr` 8))
```

Here, we have a `xor128` definition written very similar to its
original C code. We apply our transformations as described above and
put the resulting string to the standard output. That's all!

Special rule should also be added to _Makefile_ to automate using this
generator:

```make  
xsT.cc: xsTgen.hs
	$(HS) xsTgen.hs > xsT.cc
```

where makefile variable `HS` points to _runhaskell_ program from
Haskell Platform or its equivalent from other Haskell implementation
of your choice. It means that _xsT.cc_ will only be generated when
missing or when generator itself has been modified. Source code
archive for this article contains pre-generated _xsT.cc_ so you will
be able to build C++ code without installing Haskell.

### Forward Step

Jumping $k$ steps forward from RNG state $S$ if you know transition
matrix $T$ is as easy as calculating $T^kS$. (See my [previous
article](/monkey-typewriter/) for details.) So, we need implementing
these arithmetical primitives first. Let's define transition matrix
representation in _xs.h_ as:

```cpp
using tr_t = uint32_t[STATE_SIZE_EXP][STATE_SIZE];
```

and our next chunk of code in _xs.cc_ would be:

```cpp
inline uint32_t get_bit(const tr_t A, size_t row, size_t col)
...
inline void set_bit(tr_t A, size_t row, size_t col, uint32_t val)
...
void mat_mul(tr_t C, const tr_t A, const tr_t B)
...
void mat_pow(tr_t B, const tr_t A, uint64_t n)
```

This code is more or less straightforward as long as you are familiar
with matrix arithmetics,
[exponentiation by squaring](/https://en.wikipedia.org/wiki/Exponentiation_by_squaring),
and bit manipulations in C++. Therefore, I would just briefly describe
some implementation trade-offs made to keep it more clear than
efficient:

* Our matrix $T$ uses packed representation, where each integer
word corresponds to many matrix items, one bit per item. Therefore,
using bitwise logical operations in C++, we could process many items
in parallel, but in my code for this article, I traded efficiency for
clarity. That is, I simply defined access functions for one-bit matrix
items, `get_bit` and `set_bit`, and implemented trivial $O(N^3)$
matrix multiplication on top of these accessors.

* One more tradeoff was using internal temporary storage in matrix
multiplication. Without these temporary variables, caller would be
responsible for avoiding aliases, i.e., function calls like this:
`mat_mul(B, A, B)`, there multiplication result is written over one of
the arguments. Putting the caller on duty for allocating non-aliased
buffers could probably be more efficient but less clear.

* Function `mat_pow` calculating matrix power takes its exponent as
64-bit integer. This is not enough to cover all possible step lengths
for `xor128` but handling multiple-word numbers would make code less
clear.

### Backward Step

Stepping backward with transition matrices uses one more special
matrix which is inverse to transition matrix $T$. That is, we need to
calculate $T^{-1}$ such that $TT^{-1}=I$, where $I$ is [identity
matrix](https://en.wikipedia.org/wiki/Identity_matrix).

While it is possible to obtain $T^{-1}$ from $T$ solving matrix
equations or with some other matrix inversion algorithm (which could
be quite difficult to do in $\mathbb{F}_2$), there is simpler way
based on xorshift properties. We know that period of this RNG is
$2^N-1$, where $N$ is the number of state bits. It means that starting
from some state $S$ and doing $2^N-1$ steps, we should get that state
$S$ back. In our matrix notation, this can be written as $T^PS = S$,
where $P = 2^N-1$. Since $IS = S$ by definition of identity matrix, we
can conclude that $T^P = I$. We can also write the latter equation as
$TT^{P-1} = I$ and finally conclude that $T^{-1} = T^{P-1}$.

The same result can be obtained without this math. Indeed, if doing
$P$ steps brings us back to initial step $S$, what should happen if we
do one step less? We should end up with a state $S^\ast$ such that the
next step will make it $S$, i.e. $TS^\ast=S$. So, this state $S^\ast$ is a
direct precedent to $S$, as needed, and we can jump to it immediately
using $T^{P-1}$ as transition matrix.

After figuring out how $T^{-1}$ could be calculated, the
implementation is obvious. We could use our `mat_pow` function to get
an appropriate power of $T$, but as you can remember, we limited the
exponent value to be 64-bit for simplicity. So, I've implemented
separate function, `calculate_bwd_transition_matrix`, which
effectively calculates matrix exponentiation but with hardcoded
exponent value $2^N-2$. This value has simple binary representation:
$N-1$ ones followed by single zero. Therefore, code in `mat_pow`
checking for next bit being one can be replaced by test of bit
position.

This initialization should be somehow performed before we do any
backward step. In production code, we should use precalculated
constant value like we do for matrix $T$, but in our proof-of-concept
code, I'm calculating it before the first usage. This makes first
backward step significantly slower than forward one.

### Generic Matrix-based Transition Functions

For convenience, there are two functions wrapping all matrix
manipulations into simple interface:

```cpp
enum class Direction {FWD, BWD};

void prepare_transition(tr_t T_k, uint64_t k, Direction dir);

void do_transition(state_t &s, const tr_t T);
```

First function calculates transition matrix for $k$ steps using
appropriate single-step matrix depending from direction specified as
parameter. It looks attractive to make parameter $k$ signed integer
and pass negative values to jump backward, but this idea has its own
drawback: even with 128-bit integers, we will need full unsigned
precision to represent all possible jump sizes, so it would be better
to pass direction as separate parameter.

Second function, `do_transition`, simply takes matrix prepared by
`prepare_transition` and applies it to given RNG state.

Making these two transition phases separate is needed if you are going
to perform multiple jumps with the same number of steps. Since matrix
exponentiation is more expensive than applying transition matrix to
RNG state, we can precalculate transition matrix once and re-use it
many times later.

## Pit Stop

First part of the article comes to its end. We have refreshed our
background in matrix-based jump algorithms and applied it to one more
family of RNGs. Using functional language Haskell to generate
transition matrix has been very unusual but interesting
adventure. With Haskell, convenient DSL-like notation describing
xorshift formula can be developed in a few lines of code. This code
may look quite esoteric but using read-eval-print loop, the popular,
if not standard, feature of Haskell implementations, you can develop
and try it part by part and step by step.

In the second part, we will dive deeper in the mysteries of algebra
and learn how polynomials can be used to solve this problem in a
completely different way.

## References

<a id="ref.1">[1]</a> Hiroshi Haramoto, Makoto Matsumoto, and Pierre L’Ecuyer.
  A fast jump ahead algorithm for linear recurrences in a polynomial space.
  In _Proceedings of the 5th International Conference on Sequences and Their Applications_,
  SETA ’08, pages 290–298, Berlin, Heidelberg, 2008. Springer-Verlag.
  (Available [online](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/ARTICLES/jump-seta-lfsr.pdf).)

<a id="ref.2">[2]</a> Hiroshi Haramoto, Makoto Matsumoto, Takuji Nishimura, François Panneton, and Pierre L’Ecuyer.
  Efficient jump ahead for $\mathbb{F}\_2$-linear random number generators. _INFORMS Journal on Computing,_ 20(3)
:385–390, 2008.
  (Available [online](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/ARTICLES/jumpf2-printed.pdf).)


<a id="ref.3">[3]</a> Donald E. Knuth. _The Art of Computer Programming, Volume 2 (3rd Ed.): Seminumerical Algorithms._
  Addison-Wesley Longman Publishing Co., Inc., Boston, MA, USA, 1997.

<a id="ref.4">[4]</a> P. L’Ecuyer and R. Simard. TestU01: A C library for empirical testing of random number generators.
  _ACM Transactions on Mathematical Software,_ 33(4): Article 22, August 2007.

<a id="ref.5">[5]</a> George Marsaglia. Xorshift RNGs.
  _Journal of Statistical Software, Articles,_ 8(14):1–6, 2003.
  (Available [online](http://www.jstatsoft.org/v08/i14/paper).)

## History
  * 12<sup>th</sup> April, 2020: Initial version
