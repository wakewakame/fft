# Cooley-Tukey FFT

DFT (離散フーリエ変換) の定義は以下である。

```math
X\left[k\right] = \sum_{n=0}^{N-1} x\left[n\right] \exp\left(-2\pi i \frac{n k}{N}\right)
```

ここで $$X$$ の要素数 $$N$$ が 2 つの整数 $$N = N\_1 \times N\_2$$ で表せる時、 DFT は 2 つの DFT に分割することができる。

上式に

```math
\begin{align}
N &= N_1 \times N_2 & \\
k &= N_2 k_1 + k_2 & \left( 0 \le k_1 < N_1, 0 \le k_2 < N_2 \right)
\end{align}
```

を代入して整理してみると、

```math
\begin{align}
                & X\left[k\right]             &= \sum_{n=0}^{N-1} x\left[n\right] \exp\left(-2\pi i \frac{n k}{N} \right) \\
\Leftrightarrow & X\left[N_2 k_1 + k_2\right] &= \sum_{n_2=0}^{N_2-1} \sum_{n_1=0}^{N_1-1} x\left[N_1 n_2 + n_1\right] \exp\left(-2\pi i \frac{(N_1 n_2 + n_1) (N_2 k_1 + k_2)}{N_1 N_2} \right) \\
                &                             &= \sum_{n_2=0}^{N_2-1} \left( \sum_{n_1=0}^{N_1-1} \left( x\left[N_1 n_2 + n_1\right] \exp\left(-2\pi i \frac{n_1 k_2}{N_1 N_2} \right) \right) \exp\left(-2\pi i \frac{n_1 k_1}{N_1} \right) \right) \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) \exp\left(-2\pi i n_2 k_1 \right) \\
                &                             &= \sum_{n_2=0}^{N_2-1} \left( \sum_{n_1=0}^{N_1-1} \left( x\left[N_1 n_2 + n_1\right] \exp\left(-2\pi i \frac{n_1 k_2}{N_1 N_2} \right) \right) \exp\left(-2\pi i \frac{n_1 k_1}{N_1} \right) \right) \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) \\
\end{align}
```

と式変形でき、内側の DFT と外側の DFT に分割することができる。

...これだと分かりづらいので、 $$N\_1 = 4, N\_2 = 2$$ のように具体的な数を入れて考えてみる。

```math
\begin{align}
                & X\left[k\right]         &= \sum_{n=0}^{8-1} x\left[n\right] \exp\left(-2\pi i \frac{n k}{8} \right) \\
\Leftrightarrow & X\left[2 k_1 + 0\right] &= \sum_{n_2=0}^{2-1} \left( \sum_{n_1=0}^{4-1} x\left[4 n_2 + n_1\right] \exp\left(-2\pi i \frac{n_1 0}{8} \right) \exp\left(-2\pi i \frac{n_1 k_1}{4} \right) \right) \exp\left(-2\pi i \frac{n_2 0}{2} \right) \\
                &                         &= \sum_{n_2=0}^{2-1} \left( \sum_{n_1=0}^{4-1} x\left[4 n_2 + n_1\right] \exp\left(-2\pi i \frac{n_1 k_1}{4} \right) \right) \\
                &                         &= \left( \sum_{n_1=0}^{4-1} x\left[n_1\right] \exp\left(-2\pi i \frac{n_1 k_1}{4} \right) \right) + \left( \sum_{n_1=0}^{4-1} x\left[4 + n_1\right] \exp\left(-2\pi i \frac{n_1 k_1}{4} \right) \right) \\
                &                         &= \sum_{n_1=0}^{4-1} \left( x\left[n_1\right] + x\left[4 + n_1\right] \right) \exp\left(-2\pi i \frac{n_1 k_1}{4} \right) \\
                & X\left[2 k_1 + 1\right] &= \sum_{n_2=0}^{2-1} \left( \sum_{n_1=0}^{4-1} x\left[4 n_2 + n_1\right] \exp\left(-2\pi i \frac{n_1 1}{8} \right) \exp\left(-2\pi i \frac{n_1 k_1}{4} \right) \right) \exp\left(-2\pi i \frac{n_2 1}{2} \right) \\
                &                         &= \sum_{n_2=0}^{2-1} \left( \sum_{n_1=0}^{4-1} x\left[4 n_2 + n_1\right] \exp\left(-2\pi i \frac{n_1}{8} \right) \exp\left(-2\pi i \frac{n_1 k_1}{4} \right) \right) \exp\left(-2\pi i \frac{n_2}{2} \right) \\
                &                         &= \left( \sum_{n_1=0}^{4-1} x\left[n_1\right] \exp\left(-2\pi i \frac{n_1}{8} \right) \exp\left(-2\pi i \frac{n_1 k_1}{4} \right) \right) - \left( \sum_{n_1=0}^{4-1} x\left[4 + n_1\right] \exp\left(-2\pi i \frac{n_1}{8} \right) \exp\left(-2\pi i \frac{n_1 k_1}{4} \right) \right) \\
                &                         &= \sum_{n_1=0}^{4-1} \left( x\left[n_1\right] - x\left[4 + n_1\right] \right) \exp\left(-2\pi i \frac{n_1}{8} \right) \exp\left(-2\pi i \frac{n_1 k_1}{4} \right) \\
\end{align}
```

以上の分割操作を繰り返すことにより DFT の計算量を $$O(N^2)$$ から $$O(N \log(N))$$ に減らすことができる。

この操作をグラフにすると、かの有名なバタフライ演算の図となる。

```
x[0] ---+-------------------+-------------+---------> X[0]
       ^|                  ^|            ^v
x[1] --||--+---------------||--+---------+--(*-W0)--> X[4]
       || ^|               |v ^|
x[2] --||-||--+------------+--||-(*-W0)---+---------> X[2]
       || || ^|               |v         ^v
x[3] --||-||-||--+------------+--(*-W2)--+--(*-W0)--> X[6]
       |v || || ^|
x[4] --+--||-||-||-(*-W0)---+-------------+---------> X[1]
          |v || ||         ^|            ^v
x[5] -----+--||-||-(*-W1)--||--+---------+--(*-W0)--> X[5]
             |v ||         |v ^|
x[6] --------+--||-(*-W2)--+--||-(*-W0)---+---------> X[3]
                |v            |v         ^v
x[7] -----------+--(*-W3)-----+--(*-W2)--+--(*-W0)--> X[7]
```

以上をコードで表すと以下のようになる。

```
TODO
```

参考: [Cooley-Tukey\_FFT\_algorithm - Wikipedia](https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm#Variations)

# Stockham FFT

TODO

# FFT による畳み込み

TODO

# Bluestein's FFT

TODO
