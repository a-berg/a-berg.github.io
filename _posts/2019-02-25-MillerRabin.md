---
layout: post
title: "Implementación del test de Miller-Rabin en C++"
author: Adrián Berges Enfedaque
date: 25-02-2019
math: true
---

$$
\newcommand{\complexity}[1]{\mathcal{O}(#1)\ }
\newcommand{\linear}{n}
\newcommand{\logn}{\log n}
\newcommand{\loglinear}{n\,\log n}
\newcommand{\polylog}[1]{\log^{#1} n}
\newcommand{\Int}{\mathbb{Z}}
\newcommand{\groupZnZ}[1][n]{(\Int/#1\Int)^*}
$$

<!-- Contexto relativo a por qué me dio por implementar esto. -->
Recientemente empecé a leer un libro de C++ orientado a proyectos y ejercicios ([The Modern C++ Challenge][cpp-chall]).
El ejercicio 4 es un ejercicio sobre números primos, que pide:

> Dado un número natural $N$, encontrar el mayor primo $p$ tal que $p<N$.

Normalmente, este tipo de ejercicios se resuelve mediante la [criba de
eratóstenes][wiki-criba] (la solución que da el libro también es esta); dado que
suele ser un problema introductorio y no se plantea para la búsqueda de la mejor
solución posible.

No obstante, en este caso me propuse implementar algún algoritmo eficiente de
test de primalidad, ya que es un problema frecuente tanto en colecciones de ejercicios
(como el [Proyecto Euler][euler-p7]), como en aplicaciones de seguridad.

> Notación. **Complejidad asintótica**: se utiliza la notación
> $\complexity{f(n)}$ para indicar que el número de operaciones de un algoritmo
> está dominado por el término $f(n)$. También puede usarse para indicar la
> memoria que ocupa un algoritmo, en cuyo caso se indicará que es _complejidad en memoria_.

El algoritmo de Eratostenes, si se implementa de forma eficiente, tiene una
complejidad de $\complexity{n\log\log n}$. Esto lo hace inviable para rangos de
números de un tamaño normal en ordenadores actuales (enteros de 32 y 64 bits), y
ya no digamos para aplicaciones comerciales donde las claves son de 256 bits.

Tras una búsqueda rápida, encontré varios algoritmos:

* Test de Fermat
* Test de Miller-Rabin
* Test AKS
* Test por curva elíptica

El test de Fermat es muy sencillo y bastante eficiente, con un a complejidad de
$\complexity{\polylog{3}}$, pero falla para ciertos números (los [números de
Carmichael][wiki-carmichael]). Los tests de curva elíptica y AKS eran bastante
complejos para una implementación rápida y hecha por un lego; por tanto me
decanté por el test de Miller-Rabin, el cual es más robusto e igual de rápido
que el de Fermat, sin ser demasiado difícil de implementar.

Adicionalmente, me auto-impuse la meta de conseguir un algoritmo que funcionase
para enteros sin signo de 64 bits sin overflow. Esto, que se explicará más
adelante en el artículo, constituyó en sí mismo un reto bastante grande, aunque
la solución, como suele pasar a menudo, es bastante sencilla. También, ya que el
objetivo es aprender y volverme un mejor programador, quise implementar tests
con la librería [Catch2][catch2] y construir el proyecto usando Cmake. El
resultado es un pequeño [repositorio][link-github] que usaré para seguir
subiendo cosas relacionadas con algoritmos de carácter numérico en C++.

<!-- Hablar del teorema de numeros primos -->

[cpp-chall]: https://www.packtpub.com/application-development/modern-c-challenge
[wiki-criba]: https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
[euler-p7]: https://projecteuler.net/problem=7
[wiki-carmichael]: https://en.wikipedia.org/wiki/Carmichael_number
[wiki-fermat]: https://en.wikipedia.org/wiki/Fermat_primality_test
[catch2]: https://github.com/catchorg/Catch2
[link-github]: https://github.com/a-berg/cpp-numeric

# Qué es el test de Miller-Rabin
<!-- Comentar que es un test de primalidad, probabilístico. -->

El test de Miller-Rabin es un test de primalidad probabilístico, es decir: dado
un _candidato_ $n$, nos dirá si el número es compuesto con total seguridad, pero
si el algoritmo resuelve que $n$ es primo sólo podemos estar _bastante seguros_
de que este resultado es cierto.

La base del test se encuentra en el Teorema 2.1 de [este artículo][05rene]. El
teorema establece algunas de las claves del algoritmo; siendo $n$ el número que
queremos comprobar (asumiremos que es un número impar, ya que no hay primos
pares mayores que 2):

* Se trabaja en el grupo multiplicativo $\groupZnZ$, es decir el grupo de los enteros módulo $n$.
* Descomponer $n-1$ como $m2^k$, para algún $k\geq1$ y $m$ impar.
* Para un $a\in\groupZnZ$ dado, comprobar si pertenece al conjunto $B = \{a\in\groupZnZ~:~a^m = 1 \vee a^{m\,2^i} = -1~\Finv\,i\in[0,k)\}$.

<!-- , donde el símbolo $\Finv$ se ha usado como abreviación de _para algún_. -->

> Nota: siempre que se trate de números dentro del grupo $\groupZnZ​$, se
> entiende que hablar del producto $ab = c​$ es equivalente a decir $ab \equiv c\mod n​$.
> Asimismo, cualquier número $a​$ que pertenezca a este grupo cumple que $0<a<n​$.

La explicación del test se detalla a continuación.

Para cualquier $a\in\groupZnZ$, si $n$ es primo y $a^2 = 1$, entonces $a=1$ o
$a=-1$. Esto es así porque si $n$ divide a $a^2-1 = (a-1)(a+1)$, entonces se
pueden dar dos casos: que $n$ divida a $(a-1)$, en cuyo caso se tiene que
$(a-1)=0$ y por tanto $a=1$; o que $n$ divida a $(a+1)$, que de forma similar
resulta en que $a=-1$.

Por el contrapositivo de este principio, si encontramos que $a^2 = 1$ con 
$a \neq\pm 1$, entonces $n$ no puede ser primo.

El conjunto $B$ sería entonces el conjunto de los números $a$ que "_atestiguan_"
que $n$ es compuesto (por eso se llama "testigo" o "_witness_" a $a$ en la
implementación). El resultado del Teorema 2.1 establece un límite superior para
la probabilidad de un falso positivo; sin embargo, el [artículo original de Rabin][rabin-orig]
es algo más claro en su resultado (Teorema 1).

El Teorema 1 de Rabin establece que la cardinalidad de $B$ es _al menos_ tres
cuartas partes de la cardinalidad de $\groupZnZ$, y esto significa que tomando
un $a$ cualquiera, la posibilidad de que $n$ pase el test a pesar de ser
compuesto es de como mucho $0.25$. Normalmente, la probabilidad es mucho menor.

En general, para números enteros de 64 bits, la probabilidad de que un número
compuesto pase 128 rondas es tremendamente pequeña, en concreto menor que
$n^{-1}$. Si queremos estar totalmente seguros, hacer pasar $2\log^2{n}$ rondas
al número reduce esta probabilidad a $n^{-\log n}$: para números medianamente
grandes, existen más posibilidades de que se de un resultado erróneo por
problemas de fábrica del ordenador, que por el propio algoritmo.

[05rene]: http://www.math.leidenuniv.nl/~psh/ANTproc/05rene.pdf
[rabin-orig]: https://ac.els-cdn.com/0022314X80900840/1-s2.0-0022314X80900840-main.pdf?_tid=e1893542-ba98-44c6-948b-b66e2e5de342&acdnat=1550343584_cc62a6fdf0434b15d7095eb2310953b4
[mr-deterministic-proof]: http://www.ams.org/journals/mcom/2014-83-290/S0025-5718-2014-02830-5/home.html

# Boceto del algoritmo
<!-- Tras la introducción, examinar cómo se llega al algoritmo -->

Vamos a suponer que ya se ha descompuesto $n-1$ en $d2^s$, y que se ha
seleccionado un $a$ en un paso anterior. Sea:

* $x = a^d$ y
* La iteración $x \mapsto x^2$, que genera la serie $a^{d2^i}$.

Para un único test de Miller-Rabin, será necesario comprobar los siguientes
casos:

* Si $x = \pm 1$, el test es positivo ($n$ es primo).
* Durante la iteración $x \mapsto x^2$, si $x^2 = 1$ el test falla (porque ya
  eliminamos el caso $x = 1$ anteriormente).
* Por el contrario, si $x = -1$, el test ha pasado.
* Finalmente, si se llega al final de la iteración y $x = a^{n-1} \neq -1$ se incumple la condición de Fermat y el test falla.

Un esbozo del algoritmo en pseudocódigo puede ser:

```
miller-rabin-test (UInt :: n, d, s, a) -> Bool

x = a^d % n
if (x == 1 || x == -1)
    return true;

for i = 1 .. s-1:
    x = x^2 % n;

    if (x == 1) // test failed, n is not prime
        return false;

    if (x == (number - 1)) // test passed, run next one
        return true;

if (x != - 1) // test failed, n is not prime
    return false;
```

Este es el núcleo del algoritmo, no obstante, en la implementación se deberán
tener en cuenta las siguientes particularidades:

* Conviene hacer una comprobación inicial para eliminar números pares o primos
  conocidos.
* Se debe descomponer el número $n-1$.
* Hay que definir e inicializar un generador de números aleatorios para obtener $a$.
* Se debe repetir el test k veces.
* Hay que prestar atención a la multiplicación y exponenciación modulares, para
  realizarlas de forma correcta (sin overflow) y eficiente.
* Usaremos enteros sin signo, por lo que no es viable el concepto de $-1$. No
  obstante, en aritmética modular $a \equiv -1 \mod m \iff a \equiv m - 1 \mod m$.

# Implementación

Empecemos con las partes básicas: descartar primos obvios y múltiplos de éstos y
factorizar $n-1$:

```cpp
bool miller_rabin_test(const uint64_t &number, const uint64_t &acc)
{
    // initial check
    if (number == 2 || number == 3)
        return true;
    if (number % 2 == 0 || number % 3 == 0)
        return false;
    // factorise number to 2^s * d
    uint64_t d = number - 1;
    uint64_t s = 0;
    while ((d & 1) == 0 && d > 0)
    {
        d >>= 1;
        s++;
    }

    /* TODO: Miller-Rabin */
};
```

Las únicas dificultades para entender el código de arriba son el uso de
operadores bit a bit: los de desplazamiento `>>` y `<<` y los lógicos (`&`).

Los desplazamientos de bit equivalen a dividir entre o multiplicar por
**potencias de** 2; por ejemplo: `5 << 2` transformaría la secuencia `101` en
`10100`, lo que equivale a multiplicar por 4: `10100` en binario es 20 en
decimal. Igualmente, `22 >> 1` daría como resultado 11 (`10110 -> 1011`).
Por último, el operador `c >>= 1` es equivalente a `c = c >> 1`.

Por otro lado, el **y-lógico** equivale a hacer una operación modular, por
ejemplo: `1010 & 0001 = 0` (10 módulo 2 es 0). Esto es porque se hace un `AND`
bit a bit. Como 1 tiene un único bit no cero (tiene 1 dígito), se compara con el
último bit del número.

La razón para usar estas operaciones es que en números muy grandes resulta más
seguro (si no más eficiente). Durante el debug del programa me di cuenta de que
ciertos inputs muy grandes podían llegar a hacer overflow con las operaciones normales.

El siguiente paso es iniciar la distribución uniforme y hacer el bucle externo
del test.

```cpp
bool miller_rabin_test(const uint64_t &number, const uint64_t &acc)
{
    /* snip */

    // initialize distribution
    std::default_random_engine generator;
    std::uniform_int_distribution<uint64_t> distribution(2, number - 2);
    // actual miller-rabin here
    uint64_t witness;
    uint64_t x;
    for (int k = 0; k < acc; ++k) // run k miller-rabin tests
    {
        /* TODO: inner loop */
    }
};
```

Para el bucle interno, se hacen los pasos descritos en el apartado anterior

```cpp
/* Inner loop contents */
witness = distribution(generator);
x = powmod(witness, d, number);
if (x == 1 || x == number - 1) // test passed, run next one
    continue;

for (int i = 0; i < s - 1; ++i)
{
    x = powmod(x, 2, number); // <- evil modular exponentiation

    if (x == 1) // test failed, n is not prime
        return false; // break from function

    if (x == (number - 1)) // test passed, run next one
        break; // break from inner for loop into outer for loop
}

if (x != number - 1) // test failed, n is not prime
    return false;
```

Parece sencillo. Excepto que los números en la función `powmod` pueden hacerse
muy, muy grandes. Echemos un vistazo a la implementación "ingenua" de un
algoritmo clásico de exponenciación modular:

```cpp
uint64_t powmod(uint64_t base, uint64_t exponent, uint64_t modulus)
{
    if (modulus == 1)
        return 0;

    uint64_t result = 1;
    base = base % modulus;

    while (exponent > 0)
    {
        if (exponent & 1 == 0)
            result = (result * base) % modulus;
        exponent >>= 1;
        base = (base * base) % modulus; // overflows!
    }

    return result;
}
```

El algoritmo se conoce como [algoritmo de exponenciación por cuadratura][squaring-exp].

Este algoritmo se basa en el hecho de que $a^{b+c} = a^b\cdot a^c$. Además, se
observa que si un exponente $n$ es par, se podrá descomponer en
$\frac{n}{2}+\frac{n}{2}$ y si es impar, en $\frac{n}{2}+\frac{n}{2}+1$. Por
tanto el algoritmo sólo requiere ir elevando al cuadrado la base y multiplicarla
por el resultado cuando el exponente se vuelve impar.

Para adaptarlo a aritmética modular, se usa que
$a \cdot b \mod m \equiv \left((a \mod m)\cdot(b\mod m)\right)\mod m$.


[squaring-exp]: https://en.wikipedia.org/wiki/Exponentiation_by_squaring

## El problema de la multiplicación

Cuando la base es un número mayor de 32 bits (o también cuando
`log2(result*base)>64`), se produce un overflow: se intenta meter un número `c`
de más de 64 bits en el campo de 64 bits y lo que se obtiene es `c % 2^64`.

No conozco ningún teorema o resultado que permita transformar `(c % 2^64) % m`
en el valor de `c % m`, por lo que aparece un problema si se rehúsa a restringir
los valores de entrada a enteros de 32 bits o hacer una conversión de tipos en
la operación `base*base % modulus` a enteros de 128 bits.

La solución es una adaptación del [algoritmo ruso de multiplicación][russian-peasant-algo],
el cual explota el hecho de que $a\cdot b = \frac{a}{2}\cdot2\,b$. En
particular, se construyen las siguientes series:

$$
\begin{aligned}
a_i &= \lfloor a_{i-1}/2 \rfloor\\
b_i &= 2\,b_{i-1},
\end{aligned}
$$

y el conjunto 
$C=\{b_i ~:~ a_i \text{ es impar } \forall i = 1,\ldots, \lfloor\log a \rfloor\}$.
Si se realiza el sumatorio de todos los elementos de $C$, se obtiene el producto
$ab$. Esto es porque el algoritmo realmente se basa en descomponer $a$ en una
suma de potencias de 2 (su forma binaria). El algoritmo tiene una complejidad de
$\complexity{\log n}$ operaciones.

Hagamos primero la transformación del algoritmo a aritmética modular. Para ello,
se establece la siguiente igualdad:

$$
\left(\sum_{i=0}{c_i\mod m}\right)\mod m =
\left(\sum_{i=0}{c_i}\right)\mod m,
$$

donde los $c_i$ son los elementos de $C$ (pero vale para cualquier sumatorio de
números). Con esto, se tendría el siguiente algoritmo de multiplicación, que
todavía desborda:

```cpp
uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m) {
    uint64_t res = 0;
    uint64_t temp_b;

    while (a > 0) {
        if (a & 1)
            (res += b) %= m; // ≡ res = (res + b) % m

        a >>= 1; // halve a
        (b <<= 1) %= m; // double b, modulo m
    }
    return res;
}
```

Para hacerlo robusto frente a desbordamientos, se necesitan unas pequeñas
modificaciones.

Inicialmente (y sin conocer todavía el algoritmo ruso de multiplicación),
intenté resolver el problema intentando dividir el problema en múltiples casos
dependiendo de si los números que entraban al algoritmo tenían un tamaño grande
($2^{64}-\varepsilon$). Funcionaba muy bien cuando uno de los múltiplos o el
módulo eran grandes, pero no cuando ambos operandos rondaban los 32 bits.

Finalmente, la solución (como pasa a menudo) resultó ser sencilla. He de decir
que [no es mía][so-multiplication-answer], aunque la he mejorado ligeramente.

Para adaptar el algoritmo, se usa que la suma de dos números -- llamémoslos en
este caso $i<m$ y $j<m$ -- módulo $m$ será siempre $i+j$ si $i+j<m$, y $i+(j-m)$
si $i+j>m$ (respectivamente, $(i-m)+j$). Dado que siempre estamos doblando $b$ y
dividiendo $a$, el único número que debe preocuparnos en su desbordamiento es
$b$. Además, dado que siempre se toma módulo $m$, tanto el resultado $res$ como
$2b$ serán menores que $m$ en todo momento durante la ejecución del bucle.

Supongamos que queremos comprobar si la suma de dos números desborda, obviamente
no podemos realizar la operación y luego hacer la comprobación; pero sí se puede
comprobar que se cumpla la desigualdad: $i > m-j$. Sin embargo, falta hacer una
comprobación extra: siempre dependemos de que inicialmente, $b<m$. Por tanto,
hágase la sustitución de las operaciones que involucran a $b$ por estas
alternativas:

```cpp
uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m) {
    uint64_t res = 0;
    uint64_t temp_b;

    b %= m; // now b will be less than m

    while (a != 0) {
        if (a & 1) {
            // equivalent to: res = (res + b) % m
            if (b >= m - res)
                res -= m;
            res += b;
        }
        a >>= 1;

        // equivalent to b = 2*b % m
        temp_b = b;
        if (b >= m - b)
            temp_b -= m;
        b += temp_b;
    }
    return res;
}
```

La operación `b %= m` no afecta a los cálculos, porque
$a \cdot b \mod m \equiv \left((a \mod m)\cdot(b\mod m)\right)\mod m$.

[russian-peasant-algo]: https://en.wikipedia.org/wiki/Ancient_Egyptian_multiplication#Russian_peasant_multiplication
[so-multiplication-answer]: https://stackoverflow.com/a/18680280/6560267

# Complejidad

Para concluir, se va a hacer un pequeño análisis del coste computacional de la
solución.

Como se explica en el [primer artículo][05rene], la complejidad es
$\complexity{k\polylog{1+\mu}}$, con $\mu = 2$ para un algoritmo de
multiplicación normal, y $\mu = 1+\varepsilon$ cuando se usan algoritmos
rápidos. En este caso no estoy seguro de que $\mu$ corresponde al algoritmo ruso
de multiplicación, pero voy a asumir que $\mu = 2$, y que los algoritmos rápidos
son los basados en la [transformación rápida de Fourier][FFT]. Es resumen, se
tendría $\complexity{k\polylog{3}}$.

Para tener una probabilidad de error baja, se fija $k=2\log n$ y por tanto la
complejidad del algoritmo será del orden de $\complexity{\polylog{4}}$.

Adicionalmente, para encontrar el mayor primo $p$ tal que $p\leq N$, si se
empieza en $N$ y se desciende desde ahí, está claro que como mucho se deberán
hacer $\sqrt{n}$ comprobaciones. Esto hace que el tiempo total del algoritmo sea
de $\complexity{\sqrt{n}\polylog{4}}$.

Un resultado más interesante se obtiene si se aplica el Teorema de los Números
Primos ([TNP][TNP-wiki]), el cual dice que la distribución de densidad de los
números primos tiende asintóticamente a $\frac{n}{\log n}$. Esto significa que, con
números grandes, podemos esperar un promedio de $\log n$ comprobaciones antes de
encontrar el número primo buscado, y por tanto se espera un tiempo de ejecución
promedio proporcional a $\complexity{\polylog{5}}$: no está mal si se considera
que para números de 64 bits se harían unas $10^9$ operaciones, comparadas con
las $10^{21}$ operaciones que resultarían de aplicar la criba de Eratóstenes;
como curiosidad, en números de 256 bits la diferencia es de 68 órdenes de
magnitud ($10^{12}$ vs $10^{80}$).

Se dice que la complejidad del test de Miller-Rabin es polinómica en $w$ (la
_longitud_ del número), mientras que la criba es **exponencial** respecto a este
parámetro

[FFT]: http://www.ams.org/journals/mcom/1965-19-090/S0025-5718-1965-0178586-1/S0025-5718-1965-0178586-1.pdf
[TNP-wiki]: https://en.wikipedia.org/wiki/Prime_number_theorem

<!-- # Verificación

La verificación (o _testing_) del algoritmo se realizó con un test dividido en
varios casos, usando Catch2. -->
