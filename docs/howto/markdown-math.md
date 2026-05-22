# Writing Equations in Markdown

This page is a practical reference for writing mathematical equations on this
documentation site. If you are familiar with Word's equation editor but new to
LaTeX, start here.

Equations are written in **LaTeX** syntax, which is the standard for
typesetting mathematics in scientific publishing. The site renders them
automatically wherever you write them.

## Inline vs Display Equations

There are two ways to include an equation:

| Type | Syntax | When to use |
|------|--------|-------------|
| **Inline** | `$...$` | Part of a sentence, e.g. "where $x$ is…" |
| **Display** | `$$...$$` | Standalone, centred equation on its own line |

**Inline example** — the source:

```markdown
The leaf area index $L$ is the one-sided area of leaf tissue per unit
ground area (m² m⁻²).
```

renders as:

The leaf area index $L$ is the one-sided area of leaf tissue per unit
ground area (m² m⁻²).

**Display example** — the source:

```markdown
$$
E = mc^2
$$
```

renders as:

$$
E = mc^2
$$

!!! tip
    Use display equations for anything that has multiple terms, fractions, or
    needs to be easy to refer to. Use inline for short symbols and simple
    expressions that flow naturally in a sentence.

## Subscripts and Superscripts

| What | Syntax | Result |
|------|--------|--------|
| Subscript | `x_i` | $x_i$ |
| Superscript | `x^2` | $x^2$ |
| Both | `x_i^2` | $x_i^2$ |
| Multi-character subscript | `x_{ij}` | $x_{ij}$ |
| Multi-character superscript | `x^{n+1}` | $x^{n+1}$ |

Use curly braces `{}` whenever the subscript or superscript is more than one
character.

```markdown
$$
A_{ij} = \sum_{k=1}^{n} B_{ik} \, C_{kj}
$$
```

$$
A_{ij} = \sum_{k=1}^{n} B_{ik} \, C_{kj}
$$

## Greek Letters

Write `\` followed by the letter name. Capitalise the first letter for the
uppercase version.

| Symbol | Syntax | Symbol | Syntax |
|--------|--------|--------|--------|
| $\alpha$ | `\alpha` | $\Alpha$ / A | (use `A`) |
| $\beta$ | `\beta` | $\Beta$ / B | (use `B`) |
| $\gamma$ | `\gamma` | $\Gamma$ | `\Gamma` |
| $\delta$ | `\delta` | $\Delta$ | `\Delta` |
| $\epsilon$ | `\epsilon` | — | — |
| $\varepsilon$ | `\varepsilon` | — | — |
| $\zeta$ | `\zeta` | — | — |
| $\eta$ | `\eta` | — | — |
| $\theta$ | `\theta` | $\Theta$ | `\Theta` |
| $\kappa$ | `\kappa` | — | — |
| $\lambda$ | `\lambda` | $\Lambda$ | `\Lambda` |
| $\mu$ | `\mu` | — | — |
| $\nu$ | `\nu` | — | — |
| $\xi$ | `\xi` | $\Xi$ | `\Xi` |
| $\pi$ | `\pi` | $\Pi$ | `\Pi` |
| $\rho$ | `\rho` | — | — |
| $\sigma$ | `\sigma` | $\Sigma$ | `\Sigma` |
| $\tau$ | `\tau` | — | — |
| $\phi$ | `\phi` | $\Phi$ | `\Phi` |
| $\varphi$ | `\varphi` | — | — |
| $\chi$ | `\chi` | — | — |
| $\psi$ | `\psi` | $\Psi$ | `\Psi` |
| $\omega$ | `\omega` | $\Omega$ | `\Omega` |

## Fractions

```markdown
$\frac{numerator}{denominator}$
```

Examples:

| Syntax | Result |
|--------|--------|
| `\frac{a}{b}` | $\frac{a}{b}$ |
| `\frac{dA}{dt}` | $\frac{dA}{dt}$ |
| `\frac{1}{2}` | $\frac{1}{2}$ |

In display math, fractions are larger and easier to read:

```markdown
$$
R_n = \frac{G_s \, I_{dir}}{1 + k_x \cdot L}
$$
```

$$
R_n = \frac{G_s \, I_{dir}}{1 + k_x \cdot L}
$$

## Square Roots and nth Roots

```markdown
$\sqrt{x}$          % square root
$\sqrt[n]{x}$       % nth root
```

| Syntax | Result |
|--------|--------|
| `\sqrt{2}` | $\sqrt{2}$ |
| `\sqrt{x^2 + y^2}` | $\sqrt{x^2 + y^2}$ |
| `\sqrt[3]{x}` | $\sqrt[3]{x}$ |

## Sums, Products, and Integrals

These are written with limits using `_` (lower) and `^` (upper):

```markdown
$$
\sum_{i=1}^{n} x_i
$$

$$
\prod_{i=1}^{n} x_i
$$

$$
\int_{0}^{\infty} f(x) \, dx
$$
```

$$
\sum_{i=1}^{n} x_i
\qquad
\prod_{i=1}^{n} x_i
\qquad
\int_{0}^{\infty} f(x) \, dx
$$

For double integrals: `\iint`, triple: `\iiint`, and surface/contour:
`\oint`.

## Common Functions

Use `\` before function names so they are typeset upright (not italic):

| Syntax | Result | Syntax | Result |
|--------|--------|--------|--------|
| `\sin x` | $\sin x$ | `\cos x` | $\cos x$ |
| `\tan x` | $\tan x$ | `\exp(x)` | $\exp(x)$ |
| `\ln x` | $\ln x$ | `\log x` | $\log x$ |
| `\log_{10} x` | $\log_{10} x$ | `\max` | $\max$ |
| `\min` | $\min$ | `\lim_{x \to 0}` | $\lim_{x \to 0}$ |

## Operators and Relations

| Symbol | Syntax | Symbol | Syntax |
|--------|--------|--------|--------|
| $\times$ | `\times` | $\cdot$ | `\cdot` |
| $\div$ | `\div` | $\pm$ | `\pm` |
| $\leq$ | `\leq` | $\geq$ | `\geq` |
| $\neq$ | `\neq` | $\approx$ | `\approx` |
| $\propto$ | `\propto` | $\equiv$ | `\equiv` |
| $\in$ | `\in` | $\notin$ | `\notin` |
| $\rightarrow$ | `\rightarrow` | $\leftarrow$ | `\leftarrow` |
| $\infty$ | `\infty` | $\partial$ | `\partial` |
| $\nabla$ | `\nabla` | $\forall$ | `\forall` |
| $\exists$ | `\exists` | $\therefore$ | `\therefore` |

## Spacing

LaTeX removes whitespace inside math mode. Use these commands to add space:

| Command | Space |
|---------|-------|
| `\,` | thin space (most common) |
| `\;` | medium space |
| `\quad` | large space |
| `\qquad` | double large space |

```markdown
$f(x) \, dx$       % thin space before dx — conventional
$A \quad B$        % gap between terms
```

## Brackets and Delimiters

Plain brackets `( )` and `[ ]` work directly. For tall expressions, use
`\left` and `\right` to scale them automatically:

```markdown
$$
\left( \frac{a}{b} \right)^2
$$
```

$$
\left( \frac{a}{b} \right)^2
$$

Other delimiters: `\left[ \right]`, `\left\{ \right\}`,
`\left| \right|` (absolute value), `\left\| \right\|` (norm).

## Multi-line Equations with `align`

To align a derivation or a system of equations at the `=` sign, use the
`align` environment. Use `&` to mark the alignment point and `\\` to end
each line:

```markdown
$$
\begin{align}
A &= \pi r^2 \\
C &= 2 \pi r
\end{align}
$$
```

$$
\begin{align}
A &= \pi r^2 \\
C &= 2 \pi r
\end{align}
$$

A multi-step derivation:

```markdown
$$
\begin{align}
E_a &= E_{a,0} \cdot f(T) \\
    &= E_{a,0} \cdot \exp\!\left( \frac{T - T_0}{10} \right)
\end{align}
$$
```

$$
\begin{align}
E_a &= E_{a,0} \cdot f(T) \\
    &= E_{a,0} \cdot \exp\!\left( \frac{T - T_0}{10} \right)
\end{align}
$$

## Equation Numbering

Display equations in the `equation` environment are numbered automatically:

```markdown
$$
\begin{equation}
GPP = A_{canopy} \cdot \Delta t
\end{equation}
$$
```

To add a manual tag to any display equation, use `\tag{}`:

```markdown
$$
NPP = GPP - R_a \tag{1}
$$
```

$$
NPP = GPP - R_a \tag{1}
$$

## Matrices

Use the `pmatrix` (parentheses), `bmatrix` (square brackets), or `matrix`
(no brackets) environments:

```markdown
$$
\mathbf{A} =
\begin{pmatrix}
a_{11} & a_{12} \\
a_{21} & a_{22}
\end{pmatrix}
$$
```

$$
\mathbf{A} =
\begin{pmatrix}
a_{11} & a_{12} \\
a_{21} & a_{22}
\end{pmatrix}
$$

Bold upright symbols for vectors and matrices: `\mathbf{A}` → $\mathbf{A}$,
`\mathbf{v}` → $\mathbf{v}$.

## Text Inside Equations

Use `\text{...}` to include upright text inside an equation (units,
conditions, labels):

```markdown
$$
R_s \quad [\text{W m}^{-2}]
$$

$$
f(x) = \begin{cases}
1 & \text{if } x > 0 \\
0 & \text{otherwise}
\end{cases}
$$
```

$$
f(x) = \begin{cases}
1 & \text{if } x > 0 \\
0 & \text{otherwise}
\end{cases}
$$

## Quick Reference

| What | Syntax | Example result |
|------|--------|----------------|
| Inline math | `$x$` | $x$ |
| Display math | `$$...$$` | centred block |
| Subscript | `x_i` | $x_i$ |
| Superscript | `x^2` | $x^2$ |
| Fraction | `\frac{a}{b}` | $\frac{a}{b}$ |
| Square root | `\sqrt{x}` | $\sqrt{x}$ |
| Sum | `\sum_{i=1}^{n}` | $\sum_{i=1}^{n}$ |
| Integral | `\int_0^\infty` | $\int_0^\infty$ |
| Greek | `\alpha`, `\beta`, `\Gamma` | $\alpha$, $\beta$, $\Gamma$ |
| Multiply dot | `\cdot` | $\cdot$ |
| Multiply cross | `\times` | $\times$ |
| Approx | `\approx` | $\approx$ |
| Infinity | `\infty` | $\infty$ |
| Partial deriv. | `\partial` | $\partial$ |
| Auto-scaling brackets | `\left( \right)` | scales to content |
| Multi-line aligned | `\begin{align} ... \end{align}` | — |
| Manual tag | `\tag{3}` | equation label (3) |
| Bold (matrix/vector) | `\mathbf{A}` | $\mathbf{A}$ |
| Upright text | `\text{units}` | upright |

## Further Reading

A complete reference for LaTeX math commands is available at
[LaTeX Math Symbols](https://oeis.org/wiki/List_of_LaTeX_mathematical_symbols)
or the
[Wikibooks LaTeX Mathematics](https://en.wikibooks.org/wiki/LaTeX/Mathematics)
page. MathJax supports essentially the full AMS-LaTeX math command set.
