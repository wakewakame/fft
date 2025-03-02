# Cooley-Tukey FFT

DFT (離散フーリエ変換) の定義は以下である。

```math
X\left[k\right] = \sum_{n=0}^{N-1} x\left[n\right] \exp\left(-2\pi i \frac{n k}{N}\right)
```

ここで $$X$$ の要素数 $$N$$ が 2 つの整数の積 $$N = N\_1 \times N\_2$$ で表せる時、 DFT は 2 つの DFT に分割して計算することができる。

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

...これだと分かりづらいので、 $$X$$ を偶数と奇数に分けて離散フーリエ変換する例を考えてみる。  
上式に $$N\_2 = 2, k\_2 = 0, 1$$ を代入し整理してみると、

```math
\begin{align}
X\left[2 k_1 + 0\right]     &= \sum_{n_2=0}^{2-1} \left( \sum_{n_1=0}^{N_1 - 1} \left( x\left[N_1 n_2 + n_1\right] \exp\left(-2\pi i \frac{n_1 0}{2 N_1} \right) \right) \exp\left(-2\pi i \frac{n_1 k_1}{N_1} \right) \right) \exp\left(-2\pi i \frac{n_2 0}{2} \right) \\
                            &= \sum_{n_2=0}^{2-1} \left( \sum_{n_1=0}^{N_1 - 1} x\left[N_1 n_2 + n_1\right] \exp\left(-2\pi i \frac{n_1 k_1}{N_1} \right) \right) \\
                            &= \left( \sum_{n_1=0}^{N_1 - 1} x\left[n_1\right] \exp\left(-2\pi i \frac{n_1 k_1}{N_1} \right) \right) + \left( \sum_{n_1=0}^{N_1 - 1} x\left[N_1 + n_1\right] \exp\left(-2\pi i \frac{n_1 k_1}{N_1} \right) \right) \\
                            &= \sum_{n_1=0}^{N_1 - 1} \left( x\left[n_1\right] + x\left[N_1 + n_1\right] \right) \exp\left(-2\pi i \frac{n_1 k_1}{N_1} \right) \\
X\left[2 k_1 + 1\right]     &= \sum_{n_2=0}^{2-1} \left( \sum_{n_1=0}^{N_1 - 1} \left( x\left[N_1 n_2 + n_1\right] \exp\left(-2\pi i \frac{n_1 1}{2 N_1} \right) \right) \exp\left(-2\pi i \frac{n_1 k_1}{N_1} \right) \right) \exp\left(-2\pi i \frac{n_2 1}{2} \right) \\
                            &= \sum_{n_2=0}^{2-1} \left( \sum_{n_1=0}^{N_1 - 1} \left( x\left[N_1 n_2 + n_1\right] \exp\left(-2\pi i \frac{n_1}{2 N_1} \right) \right) \exp\left(-2\pi i \frac{n_1 k_1}{N_1} \right) \right) \exp\left(-2\pi i \frac{n_2}{2} \right) \\
                            &= \left( \sum_{n_1=0}^{N_1 - 1} \left( x\left[n_1\right] \exp\left(-2\pi i \frac{n_1}{2 N_1} \right) \right) \exp\left(-2\pi i \frac{n_1 k_1}{N_1} \right) \right) - \left( \sum_{n_1=0}^{N_1 - 1} \left( x\left[N_1 + n_1\right] \exp\left(-2\pi i \frac{n_1}{2 N_1} \right) \right) \exp\left(-2\pi i \frac{n_1 k_1}{N_1} \right) \right) \\
                            &= \sum_{n_1=0}^{N_1 - 1} \left( \left( x\left[n_1\right] - x\left[N_1 + n_1\right] \right) \exp\left(-2\pi i \frac{n_1}{2 N_1} \right) \right) \exp\left(-2\pi i \frac{n_1 k_1}{N_1} \right) \\
\end{align}
```

となり、要素数 $$N$$ の DFT が要素数 $$N/2$$ の DFT 2 回分に分割することができた。  
この操作をグラフにすると、かの有名なバタフライ演算の図となる。

```
                                                  ___
x[0] -----(+)------------------------------------|   |-> X[0]
         | ^                                     | D |
x[1] ----|-|------(+)----------------------------|   |-> X[2]
       +---+     | ^                             | F |
x[2] --|-|-------|-|------(+)--------------------|   |-> X[4]
       | |     +---+     | ^                     | T |
x[3] --|-|-----|-|-------|-|------(+)------------|___|-> X[6]
       | v     | |     +---+     | ^              ___
x[4] ---(-)----|-|-----|-|-------|-|---(*W_8^0)--|   |-> X[1]
               | v     | |     +---+             | D |
x[5] -----------(-)----|-|-----|-|-----(*W_8^1)--|   |-> X[3]
                       | v     | |               | F |
x[6] -------------------(-)----|-|-----(*W_8^2)--|   |-> X[5]
                               | v               | T |
x[7] ---------------------------(-)----(*W_8^3)--|___|-> X[7]
```

以上の分割操作を繰り返すことにより DFT の計算量を $$O(N^2)$$ から $$O(N \log N)$$ に減らすことができる。

最後にサンプルコードを示す。

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
