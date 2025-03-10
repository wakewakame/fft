# Cooley-Tukey FFT

配列長 $N$ の配列 $x$ があったとき、その DFT (離散フーリエ変換) の結果 $X$ は

$$
\begin{aligned}
& X[k] =
  \sum_{n=0}^{N-1} x[n]
  \exp\left(-2\pi i \frac{n k}{N}\right) &
& (0 \le k < N) &
\end{aligned}
$$

で定義される。  
ここで $x$ の配列長 $N$ が 2 つの整数の積 $N = N_1 \times N_2$ で表せる時、 DFT は 2 つの DFT に分割して計算することができる。

上式に

$$
\begin{aligned}
& N &=& N_1 \times N_2 & & & \\
& k &=& N_1 k_2 + k_1 &
& ( 0 \le k_1 < N_1, 0 \le k_2 < N_2 ) &
\end{aligned}
$$

を代入して整理すると

$$
\begin{aligned}
& X[k] &=&
  \sum_{n=0}^{N-1}
    x[n]
    \exp\left(-2\pi i \frac{n k}{N} \right) & \\
& &=&
  \sum_{n=0}^{N_1 N_2 - 1}
    x[n]
    \exp\left(-2\pi i \frac{n k}{N_1 N_2} \right) & \\
& &=&
  \sum_{n_2=0}^{N_2-1}
    \sum_{n_1=0}^{N_1-1}
      x[N_2 n_1 + n_2]
      \exp\left(-2\pi i \frac{(N_2 n_1 + n_2) k}{N_1 N_2} \right) & \\
\Leftrightarrow & X[N_1 k_2 + k_1] &=&
  \sum_{n_2=0}^{N_2-1}
    \sum_{n_1=0}^{N_1-1}
      x[N_2 n_1 + n_2]
      \exp\left(-2\pi i \frac{(N_2 n_1 + n_2) (N_1 k_2 + k_1)}{N_1 N_2} \right) & \\
& &=&
  \sum_{n_2=0}^{N_2-1}
    \sum_{n_1=0}^{N_1-1}
      \left(
        x[N_2 n_1 + n_2]
        \exp\left(-2\pi i \frac{n_2 k_1}{N_1 N_2} \right)
      \right)
      \exp\left(-2\pi i \frac{n_1 k_1}{N_1} \right)
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right)
    \exp(-2\pi i n_1 k_2 ) & \\
& &=&
  \sum_{n_2=0}^{N_2-1}
    \sum_{n_1=0}^{N_1-1}
      \left(
        x[N_2 n_1 + n_2]
        \exp\left(-2\pi i \frac{n_2 k_1}{N_1 N_2}\right)
      \right)
      \exp\left(-2\pi i \frac{n_1 k_1}{N_1} \right)
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) & \\
\end{aligned}
$$

と式変形でき、よく見ると外側の DFT と内側の DFT に分割された形になっていることがわかる。

...これだと分かりづらいので、 $X$ を偶数個目と奇数個目に分けて離散フーリエ変換する例を考えてみる。  
上式に $N_1 = 2, k_1 = 0, 1$ を代入し整理してみると、

$$
\begin{aligned}
& X[2 k_2 + 0] &=&
  \sum_{n_2=0}^{N_2-1}
    \sum_{n_1=0}^{2-1}
      \left(
        x[N_2 n_1 + n_2]
        \exp\left(-2\pi i \frac{0}{2 N_2}\right)
      \right)
      \exp\left(-2\pi i \frac{0}{N_1} \right)
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) & \\
& &=&
  \sum_{n_2=0}^{N_2-1}
    \left(
      \sum_{n_1=0}^{2-1}
        x[N_2 n_1 + n_2]
    \right)
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) & \\
& &=&
  \sum_{n_2=0}^{N_2-1}
    (x[n_2] + x[N_2 + n_2])
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) & \\
& X[2 k_2 + 1] &=&
  \sum_{n_2=0}^{N_2-1}
    \sum_{n_1=0}^{2-1}
      \left(
        x[N_2 n_1 + n_2]
        \exp\left(-2\pi i \frac{n_2}{2 N_2}\right)
      \right)
      \exp\left(-2\pi i \frac{n_1}{2} \right)
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) & \\
& &=&
  \sum_{n_2=0}^{N_2-1}
    \left(
      \sum_{n_1=0}^{2-1}
        \left(
          x[N_2 n_1 + n_2]
        \right)
        \exp\left(-2\pi i \frac{n_1}{2} \right)
    \right)
    \exp\left(-2\pi i \frac{n_2}{2 N_2}\right)
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) & \\
& &=&
  \sum_{n_2=0}^{N_2-1}
    \left(
      x[n_2] \exp(0) +
      x[N_2 + n_2] \exp(-\pi i)
    \right)
    \exp\left(-2\pi i \frac{n_2}{2 N_2}\right)
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) & \\
& &=&
  \sum_{n_2=0}^{N_2-1}
    \left(
      x[n_2] -
      x[N_2 + n_2]
    \right)
    \exp\left(-2\pi i \frac{n_2}{2 N_2}\right)
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) & \\
\end{aligned}
$$

となり、配列長 $N$ の DFT が配列長 $N/2$ の DFT 2 回分に分割することができた。  
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

以上の分割操作を繰り返すことにより DFT の計算量を $O(N^2)$ から $O(N \log N)$ に減らすことができる。

最後にサンプルコードを示す。

```
TODO
```

参考

- [Cooley-Tukey\_FFT\_algorithm - Wikipedia](https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm#Variations)
- [Cooley-Tukey のアルゴリズム - OTFFT: High Speed FFT Library](http://wwwa.pikara.ne.jp/okojisan/stockham/cooley-tukey.html)

# Stockham FFT

TODO

# FFT による巡回畳み込み

配列長 $N$ の配列 $a$ と $b$ があったとき、その巡回畳み込み $c$ は

$$
\begin{aligned}
& c[n] =
  \sum_{m=0}^{N - 1}
    a[m] \cdot b[(n-m) \mod N] &
  & ( 0 \le n < N ) &
\end{aligned}
$$

で定義される。  
そして、このときの $a, b, c$ は離散フーリエ変換 $\mathcal{F}$ を用いて以下が成り立つことが知られている。

$$
\mathcal{F}(c) = \mathcal{F}(a) \cdot \mathcal{F}(b)
$$

証明:

$$
\begin{aligned}
& \mathcal{F}(c) &=&
  \sum_{n=0}^{N-1}
    c[n]
    \exp\left(-2\pi i \frac{n k}{N}\right) & \\
& &=&
  \sum_{n=0}^{N-1}
    \left(
    \sum_{m=0}^{N - 1}
      a[m] \cdot b[(n-m) \mod N]
    \right)
    \exp\left(-2\pi i \frac{n k}{N}\right) & \\
& &=&
  \sum_{m=0}^{N-1}
    a[m]
    \sum_{n=0}^{N - 1}
      b[(n-m) \mod N]
    \exp\left(-2\pi i \frac{n k}{N}\right) & \\
& & & \downarrow l = n - m\ と置くと & \\
& &=&
  \sum_{m=0}^{N-1}
    a[m]
  \sum_{l=0}^{N - 1}
    b[l]
    \exp\left(-2\pi i \frac{(m+l) k}{N}\right) & \\
& &=&
  \sum_{m=0}^{N-1}
    a[m]
    \exp\left(-2\pi i \frac{m k}{N}\right)
  \sum_{l=0}^{N - 1}
    b[l]
    \exp\left(-2\pi i \frac{l k}{N}\right) & \\
& &=& \mathcal{F}(a) \cdot \mathcal{F}(b) &
\end{aligned}
$$

これを利用して FFT を使えば巡回畳み込みを $O(N \log N)$ で計算できる。

$$
c[n] = \mathcal{F}^{-1}(\mathcal{F}(a) \cdot \mathcal{F}(b))[n]
$$

ただし、Cooley-Tukey FFT では配列長が合成数になっている必要があった。  
計算結果は変えずに配列長を $N'$ に伸ばしたい場合、 $a'$ は $a$ の末尾に $0$ を追加するだけで良いが、 $b'$ はやや複雑である。  
$b'$ は $b'[(n-m) \mod N']\ (-(N-1) \le n-m \le N-1)$ の範囲が元々の計算に利用されていた。  
言い換えると $b'$ の末尾 $N-1$ 個と先頭 $N$ 個の範囲が利用されるため、この要素は変更せず中央に $0$ を追加することとなる。  

$$
\begin{aligned}
& a &= [&a_{0}&, &a_{1}&, &a_{2}&] \\
& b &= [&b_{0}&, &b_{1}&, &b_{2}&] \\
& a' &= [&a_{0}&, &a_{1}&, &a_{2}&, &0&, &...& &&&&, &0&] \\
& b' &= [&b_{0}&, &b_{1}&, &b_{2}&, &0&, &...&, &0&, &b_{1}&, &b_{2}&]
\end{aligned}
$$

つまり $N'$ は最低でも $2N-1$ 以上にする必要があるということになる。

参考

- [Convolution - Wikipedia](https://en.wikipedia.org/wiki/Convolution)

# Bluestein's FFT

Cooley-Tukey FFT では配列長が素数などの場合は計算量を削減することができなかった。  
Bluestein's FFT は Cooley-Tukey FFT より数倍遅いものの、どのような配列長であっても計算量 $O(N \log N)$ でフーリエ変換することができる。

配列長 $N$ の配列 $x$ があったとき、その DFT (離散フーリエ変換) の結果 $X$ は

$$
\begin{aligned}
& X[k] =
  \sum_{n=0}^{N-1} x[n]
  \exp\left(-2\pi i \frac{n k}{N}\right) &
& (0 \le k < N) &
\end{aligned}
$$

で定義される。  
ここで $nk$ を

$$
nk =
  \frac{-(k-n)^2}{2} +
  \frac{n^2}{2} + \frac{k^2}{2}
$$

のように展開し整理すると

$$
\begin{aligned}
& X[k] &=&
  \sum_{n=0}^{N-1}
    x[n]
    \exp\left(-2\pi i \frac{n k}{N} \right) & \\
& &=&
  \exp\left(-\pi i \frac{k^2}{N} \right)
  \sum_{n=0}^{N-1}
    x[n]
    \exp\left(-\pi i \frac{n^2}{N} \right)
  \exp\left(\pi i \frac{(k-n)^2}{N} \right) & \\
\end{aligned}
$$

となる。さらに

$$
\begin{aligned}
& a[n] &=&
  x[n]
  \exp\left(-\pi i \frac{n^2}{N} \right) & \\
& b[n] &=&
  \exp\left(\pi i \frac{n^2}{N} \right) &
\end{aligned}
$$

と置き換えると

$$
\begin{aligned}
& X[k] =
  b[k]^*
  \times
  \left(
    \sum_{n=0}^{N-1}
      a[n]
      b[(k-n) \mod 2N]
  \right) &
& (b[k]^*\ は\ b[k]\ の複素共役) &
\end{aligned}
$$

となり、巡回畳み込みと同じ形となる。  
よって、これは Cooley-Tukey FFT を使った畳み込み手法を用いて計算量 $O(N \log N)$ で計算することができる。  

ただし、 Cooley-Tukey で計算するにあたり 2 つ注意点がある。

- $a, b$ を FFT に入力するため、配列サイズが 2 の累乗となるよう末尾に 0 を追加する必要がある
- $b$ はそのままだとマイナス範囲へのアクセスがあるので、マイナス範囲を末尾に移動させる必要がある

**TODO: 以下の文章は書きかけ。もう少しわかりやすく直したい。**

...だた、これだと $a, b$ の配列サイズが合成数ではないケースもあるので、そのままだと Cooley-Tukey FFT を使った畳み込みができない。
そのため、配列の数が $2^N$ などになるよう、配列末尾に $0$ を追加して配列サイズを稼ぐ。

...また、これだと $b[k-n]$ は $-(N-1) \le k-n < N$ の範囲で利用されるので、マイナス範囲へのアクセスがコードに落とし込みづらい。  
$b$ のマイナスの範囲を末尾に移動させた $b'$ を考える。

$$
\begin{aligned}
& b[n] &=& [b_{-2}, b_{-1}, &b_{0}&, b_{1}, b_{2}, 0, 0, 0] & \\
& b'[n] &=& [               &b_{0}&, b_{1}, b_{2}, 0, 0, 0, b_{-2}, b_{-1}] & \\
& &=& [                     &b_{0}&, b_{1}, b_{2}, 0, 0, 0, b_{2}, b_{1}] &
\end{aligned}
$$

こうすることで、 $b[n] = b'[n \mod N_{b'}]$ と表現でき、配列のマイナスアクセスが防げるようになる。

$$
\begin{aligned}
& X[k] =
  b[k]^*
  \times
  \left(
    \sum_{n=0}^{N-1}
      a[n]
      b'[(k-n) \mod N_b]
  \right) &
& (b[k]^*\ は\ b[k]\ の複素共役) &
\end{aligned}
$$

最後にサンプルコードを示す。

```
TODO
```

参考

- [Chirp Z-transform - Wikipedia](https://en.wikipedia.org/wiki/Chirp_Z-transform#Bluestein's_algorithm)
