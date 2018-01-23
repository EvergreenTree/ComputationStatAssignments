$$\newcommand{\dd}{\mathrm d}\newcommand{\N}{\mathbb N} \newcommand{\R}{\mathbb R} \newcommand{\Q}{\mathbb Q} \newcommand{\P}{\mathbb P} \newcommand{E}{\mathbb E} \newcommand{\B}{\mathcal B} \newcommand{\F}{\mathcal F} \newcommand{A}{\mathcal A}$$
# Probability, Stochastic Processes, Lecture 1

giambattista.giacomin@math.univ-paris-diderot.fr

Time: Thusday (Changing room), Thursday (SG)

### content:

- convergence in law 
- conditional expectation
- martinales
- markov chains

Recall:

$(\Omega,\A,\P)$, where $\A \subset \mathcal{P}(\Omega)$ is a $\sigma$-algebra and $card(\mathcal P(\Omega))=2^{card(\Omega)}$. $A\in \A$ is an event.$\P$ is a probability mesure

$\A$ is a $\sigma$-algebra if 
- $\emptyset \in \A$
- $A \in \A \implies A^C \in \A$
- $\{A_i\}_{i \in \N}, A_{i} \in \A \implies \bigcup_{j \in \N} A_i \in \A$ 

###definition of "$\lim A_n$"
(reasonable because $\lim_{n\rightarrow\infty}𝟙_A(\omega)\in \{0,1\}$)

- $A_n \nearrow \implies \lim A_n = \bigcup A_n$
- $A_n \searrow \implies \lim A_n = \bigcap A_n$

####observation: 

$\A_j \prec \A$ ($\A_j$ is sub $\sigma$-algebra of $\A$) then $\bigcap{\A}$ is a $\sigma$-algebra.

$$\sigma(\mathcal{C}) = \bigcap_{\tilde{\A} \prec \A, \mathcal{C}\in \tilde \A} \tilde \A$$

###Monotone class Thm: 

$\B(\Omega) = \sigma(\Pi)$

where $\B(\Omega) = \sigma(\mathcal{O}(\Omega))$ and $\Pi$ is stable by intersection

###Random Variables:

$X: \Omega \longrightarrow E$ is a r.v. if $X^{-1}(B) \in \A, \forall B \in \mathcal{E}$ or $X^{-1}: \mathcal{E}\longrightarrow \A$

### Distribution Function:

$F_X: \R \longrightarrow [0,1]: F_X(x)= \P(X\le x) = \P(-\infty,x))$

- if $X$ is a r.v, $F$ satisfies 3 conditions (non-decreasing, 0-1, "CADLAG",右连续左极限)
- if $F$ satisfies the 3 conditions, then $\exists X$ s.t. $F = F_X$


*proof of the 2nd part:*

suppose $U\sim \mathcal{U}(0,1)$, $\P_U = \lambda$, $ F^\leftarrow(U)$defines a r.v. on $[0,1)]$

$\P(F^\leftarrow(U)\le x) = \lambda(\{w\in \Omega:F^\leftarrow(U)\le x)\})$



*exercise:* $U_1,U_2 \sim \mathcal{U}(0,1)$ then $\P(U_1=U_2) = 0$

- *CHALLENGE*: construct $U_1,U_2 \sim \mathcal{U}(0,1)$ independent defined on the same space (actually countably many can be defined)

- *CHALLENGE*: prove $|\int f(x)\mu(dx)| < \int |f(x)|\mu(dx)$

  proof:

   $\begin{align}|\int f(x)\mu(dx)| &=\alpha\int f(x)\mu(dx)\\&=\int \alpha f(x)\mu(dx)\end{align}$

###Inverse

- $F^\leftarrow (w) = \inf \{z:F(z)\ge w\}$ left-continuous
- $F^\rightarrow (w) = \inf \{z:F(z)> w\}$ right-continuous

are equivalent because discontinuous point

### Characteristic Function

- satisfies 4 conditions ($\le 1$, unif $C_0$, $\Phi_{-X} = \bar\Phi_{X}$)
- if  $\Phi$ satisfies the 4 conditions (weaker condition for 2), then it's a Characteristic Function of some r.v.

Inverse formula:  

- $\forall a<b$,

  $$\begin{align}\lim_{T\rightarrow\infty}\frac{1}{2\pi}\int_{-T}^T\frac{e^{-ita}-e^{itb}}{it}\Phi_X(t)\dd t&=\P(X\in(a,b))+\frac{1}{2}\P(X=a)+\frac{1}{2}\P(X=b)\\&=F_X(b)-F_X(a)+\frac{1}{2}(\P(X=a)-\P(X=b))\end{align}$$

  (i.e. The characteristic function characters the density)

- $\forall \Phi_X\in L^1, \exists p(\cdot)$ s.t.$\frac{1}{2\pi}\int e^{-itx}\Phi_X(t)\dd t = p(x)$

*<u>proof</u>:* 

$$\begin{align}Q(T)&=\frac{1}{2\pi}\int_{-T}^T\frac{e^{-ita}-e^{itb}}{it}\Phi_X(t)\dd t\\&=\frac{1}{2\pi}\int_{-T}^T\frac{e^{-ita}-e^{itb}}{it}\int e^{itx}\P_X(\dd x)\dd t\\&\stackrel{\text{Fubini}}{=}\frac{1}{2\pi}\int_\R\underbrace{\int_{-T}^{T}\frac{e^{i(x-a)t}-e^{i(x-b)t}}{it}\dd t}_{q_T(a,b;x)} \P_X(\dd x)\end{align}$$

$$\begin{align}q_T(a,b;x)&\stackrel{\text{devide by Re and Im}}{=}2\int_{0}^{T}\frac{\sin(x-a)t-\sin(x-b)t}{t}\dd t\\&=2\begin{cases}-S(|X-a|T)+S(|X-b|T)&\text{if }x<a\\S(|X-a|T)+S(|X-b|T)&\text{if }a\le x\le b\\S(|X-a|T)-S(|X-b|T)&\text{if }x>b\end{cases}\\&=\begin{cases}\pi&\text{if }x=a\\2\pi&\text{if }a<x<b\\\pi&\text{if }x=b\\0&\text{else}\end{cases}\end{align}$$

fgor $\alpha > 0$, define $$\text{Si}(\alpha x)=\int_{0}^{T}\frac{\sin(\alpha x)}{x}\dd x=\int_{0}^{T}\text{sinc}\ \!{\alpha x}\dd x$$, $\text{Si}$ is an odd function 

if we take the limit $T\to 0$, $$Q(T)=\lim_{T\to 0}\frac{1}{2\pi}\int q_T(a,b;x)\P(\dd x)\stackrel{\|S\|_\infty<\infty}{=}\frac{1}{2\pi}\int\lim_{T\to 0}q_T(a,b;T)\P(\dd x)=\P((a,b))+\frac{1}{2}\P(\{a\})+\frac{1}{2}\P(\{b\})$$, thus conclude.$\square$

*<u>more:</u>* we proved that for all continuous $x$ of $F_X$, $$F_X(b)-F_X(a)=\P((a,b]) = $$

furthermore, $F_X$ is $\mathcal{C}^1$ (for the density exists)

$$\lim_{b\searrow a}\frac{F(b)-F(a)}{b-a}=\lim_{b\searrow a}\frac{1}{2\pi}\int_\R \frac{e^{-ita}-e^{iab}}{it(b-a)}\Phi_X(t)\dd t \stackrel{\text{Dominate cv}}{=} \frac{1}{2\pi}\int_\R e^{-ita}\Phi_X(t)\dd t$$, which is continuous for $\Phi$ is continuous.



c.f. the proof in the note: are they equivalent? Yes.



### Convergence in Law

- Alternative: cv of distribution function $F_X$: not working!

  *<u>counter example:</u>* $X_n=\frac{1}{n}$, then $F_{X_n}=\mathbb{1}_{[\frac{1}{n},\infty)}\to \mathbb{1}_{(0,\infty)}$, which is not a distribution function.

  but as we have seen last semester, cv in law is equivalent to $\forall\alpha\in ${continuous pts of $F_X$} (alternatively, ),$F_{X_n}(\alpha)\to F_X(\alpha)$

- Why must we use $h\in\mathcal{C}_b^0(\R,\R)$ as test function , and must not exclude continuous condition

  e*<u>counter example:</u>* same as above, take $h = \mathbb{1}_{(-\infty,0]}$, we have $\E h(X_n)\nrightarrow\E h(X)$




$\implies$:$F_X(t-\varepsilon)<F_{X_n}(t)<F_X(t+\varepsilon)$

