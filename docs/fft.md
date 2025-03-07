# Cooley-Tukey FFT

DFT (離散フーリエ変換) の定義は以下である。

```math
X\left[k\right] =
  \sum_{n=0}^{N-1} x\left[n\right]
  \exp\left(-2\pi i \frac{n k}{N}\right)
```

ここで $$X$$ の要素数 $$N$$ が 2 つの整数の積 $$N = N\_1 \times N\_2$$ で表せる時、 DFT は 2 つの DFT に分割して計算することができる。

上式に

```math
\begin{align}
N &= N_1 \times N_2 & \\
k &= N_1 k_2 + k_1 & \left( 0 \le k_1 < N_1, 0 \le k_2 < N_2 \right)
\end{align}
```

を代入して整理してみると、

```math
\begin{align}
& X\left[k\right] &=
  \sum_{n=0}^{N-1}
    x\left[n\right]
    \exp\left(-2\pi i \frac{n k}{N} \right) \\
& &=
  \sum_{n_1=0}^{N_1-1}
    \sum_{n_2=0}^{N_2-1}
      x\left[N_2 n_1 + n_2\right]
      \exp\left(-2\pi i \frac{(N_2 n_1 + n_2) k}{N_1 N_2} \right) \\
\Leftrightarrow & X\left[N_1 k_2 + k_1\right] &=
  \sum_{n_1=0}^{N_1-1} 
    \sum_{n_2=0}^{N_2-1}
      x\left[N_2 n_1 + n_2\right]
      \exp\left(-2\pi i \frac{(N_2 n_1 + n_2) (N_1 k_2 + k_1)}{N_1 N_2} \right) \\
& &=
  \sum_{n_1=0}^{N_1-1}
    \sum_{n_2=0}^{N_2-1}
      \left(
        x\left[N_2 n_1 + n_2\right]
        \exp\left(-2\pi i \frac{n_2 k_1}{N_1 N_2} \right)
      \right)
      \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right)
    \exp\left(-2\pi i \frac{n_1 k_1}{N_1} \right)
    \exp\left(-2\pi i n_1 k_2 \right) \\
& &=
  \sum_{n_1=0}^{N_1-1}
    \sum_{n_2=0}^{N_2-1}
      \left(
        x\left[N_2 n_1 + n_2\right]
        \exp\left(-2\pi i \frac{n_2 k_1}{N_1 N_2}\right)
      \right)
      \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right)
    \exp\left(-2\pi i \frac{n_1 k_1}{N_1} \right) \\
\end{align}
```

と式変形でき、内側の DFT と外側の DFT に分割することができる。

...これだと分かりづらいので、 $$X$$ を偶数と奇数に分けて離散フーリエ変換する例を考えてみる。  
上式に $$N\_1 = 2, k\_1 = 0, 1$$ を代入し整理してみると、

```math
\begin{align}
X\left[2 k_2 + 0\right] &=&
  \sum_{n_1=0}^{2-1}
    \sum_{n_2=0}^{N_2-1}
      \left(
        x\left[N_2 n_1 + n_2\right]
        \times 1
      \right)
      \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right)
    \times 1 \\
&=&
  \sum_{n_2=0}^{N_2-1}
    x\left[n_2\right]
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) + \\
&&
  \sum_{n_2=0}^{N_2-1}
    x\left[N_2 + n_2\right]
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) \\
&=&
  \sum_{n_2=0}^{N_2-1}
    \left(x\left[n_2\right] + x\left[N_2 + n_2\right]\right)
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) \\
X\left[2 k_2 + 1\right] &=&
  \sum_{n_1=0}^{2-1}
    \sum_{n_2=0}^{N_2-1}
      \left(
        x\left[N_2 n_1 + n_2\right]
        \exp\left(-2\pi i \frac{n_2}{2 N_2}\right)
      \right)
      \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right)
    \exp\left(-2\pi i \frac{n_1}{2} \right) \\
&=&
  \sum_{n_2=0}^{N_2-1}
    \left(
      x\left[n_2\right]
      \exp\left(-2\pi i \frac{n_2}{2 N_2}\right)
    \right)
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) \times 1 + \\
&&
  \sum_{n_2=0}^{N_2-1}
    \left(
      x\left[N_2 + n_2\right]
      \exp\left(-2\pi i \frac{n_2}{2 N_2}\right)
    \right)
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) \times -1 \\
&=&
  \sum_{n_2=0}^{N_2-1}
    \left(
      \left(x\left[n_2\right] - x\left[N_2 + n_2\right]\right)
      \exp\left(-2\pi i \frac{n_2}{2 N_2}\right)
    \right)
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right)\\
\end{align}
```

となり、要素数 $$N$$ の DFT が要素数 $$N/2$$ の DFT 2 回分に分割することができた。  
この操作をグラフにすると、かの有名なバタフライ演算の図となる。

```
                                                                                  ___
x[0] --+-(+)---------------------------------------------------------------------|   |-> X[0]
       |  ^                                                                      | D |
x[1] --|--|------------+-(+)-----------------------------------------------------|   |-> X[2]
       +-----------+   |  ^                                                      | F |
x[2] -----|--------|---|--|------------+-(+)-------------------------------------|   |-> X[4]
          |        |   +-----------+   |  ^                                      | T |
x[3] -----|--------|------|--------|---|--|------------+-(+)---------------------|___|-> X[6]
          |        v      |        |   +-----------+   |  ^                       ___
x[4] -----+-(*-1)-(+)-----|--------|------|--------|---|--|------------(*W_8^0)--|   |-> X[1]
                          |        v      |        |   +-----------+             | D |
x[5] ---------------------+-(*-1)-(+)-----|--------|------|--------|---(*W_8^1)--|   |-> X[3]
                                          |        v      |        |             | F |
x[6] -------------------------------------+-(*-1)-(+)-----|--------|---(*W_8^2)--|   |-> X[5]
                                                          |        v             | T |
x[7] -----------------------------------------------------+-(*-1)-(+)--(*W_8^3)--|___|-> X[7]
```

以上の分割操作を繰り返すことにより DFT の計算量を $$O(N^2)$$ から $$O(N \log N)$$ に減らすことができる。

最後にサンプルコードを示す。

```
TODO
```

参考

- [Cooley-Tukey\_FFT\_algorithm - Wikipedia](https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm#Variations)
- [Cooley-Tukey のアルゴリズム - OTFFT: High Speed FFT Library](http://wwwa.pikara.ne.jp/okojisan/stockham/cooley-tukey.html)

# Stockham FFT

TODO

# FFT による畳み込み

TODO

# Bluestein's FFT

Cooley-Tukey FFT では要素数が $$N = N\_1 \times N\_2$$ などの合成数になる配列しか計算量を削減することができなかった。  
Bluestein's FFT は Cooley-Tukey FFT より遅いものの、要素数が素数など任意の配列長を計算量 $$O(N \log N)$$ でフーリエ変換することができる。

まず、DFT (離散フーリエ変換) の定義は以下である。

```math
X\left[k\right] =
  \sum_{n=0}^{N-1} x\left[n\right]
  \exp\left(-2\pi i \frac{n k}{N}\right)
```

ここで $$nk$$ を

```math
\begin{align}
nk &=
  \frac{-(k-n)^2}{2} +
  \frac{n^2}{2} + \frac{k^2}{2} \\
\end{align}
```

のように展開し整理すると、

```math
\begin{align}
& X\left[k\right] &=
  \sum_{n=0}^{N-1}
    x\left[n\right]
    \exp\left(-2\pi i \frac{n k}{N} \right) \\
& &=
  \exp\left(-\pi i \frac{k^2}{N} \right)
  \sum_{n=0}^{N-1}
    \left(
      x\left[n\right]
      \exp\left(-\pi i \frac{n^2}{N} \right)
    \right)
  \exp\left(\pi i \frac{(k-n)^2}{N} \right) \\
\end{align}
```

となる。さらに

```math
\begin{align}
a\left[n\right] &=
  x\left[n\right]
  \exp\left(-\pi i \frac{n^2}{N} \right) \\
b\left[n\right] &=
  \exp\left(\pi i \frac{n^2}{N} \right)
\end{align}
```

と置き換えると

```math
\begin{align}
& X\left[k\right] &=
  \left(b\left[k\right]^*\right)
  \left(
    \sum_{n=0}^{N-1}
      a\left[n\right]
      b\left[k-n\right]
  \right)
\end{align}
```

となり、これは畳み込みの形と同じとなる。
よって、これは Cooley-Tukey FFT の畳み込みを用いて計算量 $$O(N \log N)$$ で計算することができる。

参考

- [Chirp Z-transform - Wikipedia](https://en.wikipedia.org/wiki/Chirp_Z-transform#Bluestein's_algorithm)
