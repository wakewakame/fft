# Cooley-Tukey FFT

配列長 $N$ の配列 $x$ があったとき、その DFT (離散フーリエ変換) の結果 $X$ は

```math
\begin{aligned}
& X[k] =
  \sum_{n=0}^{N-1} x[n]
  \exp\left(-2\pi i \frac{n k}{N}\right) &
& (0 \le k < N) &
\end{aligned}
```

で定義される。  
ここで $x$ の配列長 $N$ が 2 つの整数の積 $N = N_1 \times N_2$ で表せる時、 DFT は 2 つの DFT に分割して計算することができる。

上式に

```math
\begin{aligned}
& N &=& N_1 \times N_2 & & & \\
& k &=& N_1 k_2 + k_1 &
& ( 0 \le k_1 < N_1, 0 \le k_2 < N_2 ) &
\end{aligned}
```

を代入して整理すると

```math
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
```

と式変形でき、よく見ると外側の DFT と内側の DFT に分割された形になっていることがわかる。

...これだと分かりづらいので、 $X$ を偶数個目と奇数個目に分けて離散フーリエ変換する例を考えてみる。  
上式に $N_1 = 2, k_1 = 0, 1$ を代入し整理してみると、

```math
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
```

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

```rust
fn fft(x: &Vec<(f64, f64)>, inverse: bool) -> Vec<(f64, f64)> {
    // 事前に w を計算しておく
    let sign = if inverse { 1.0 } else { -1.0 };
    let w: Vec<(f64, f64)> = (0..x.len())
        .map(|i| {
            let theta = sign * std::f64::consts::TAU * i as f64 / x.len() as f64;
            (theta.cos(), theta.sin())
        })
        .collect();

    // 出力用 (x1 と x2 を交互に使い回す)
    let mut x1 = x.clone();
    let mut x2 = vec![(0f64, 0f64); x1.len()];

    // バタフライ演算
    let mut stride = 1;
    let mut len = x.len();
    while len > 1 {
        // x2 を 0 で初期化
        for k in 0..x2.len() {
            x2[k] = (0.0, 0.0);
        }

        // len を 2 つの整数の積に分解
        let len1 = {
            const FACTORS: [usize; 6] = [7, 6, 5, 4, 3, 2];
            FACTORS.into_iter().find(|f| len % f == 0).unwrap_or(len)
        };
        let len2 = len / len1;

        // バタフライ演算
        for k1 in 0..len1 {
            for n1 in 0..len1 {
                for n2 in 0..len2 {
                    let k = len1 * n2 + k1;
                    let n = len2 * n1 + n2;
                    let w = w[(stride * n * k1) % w.len()];
                    for offset in 0..stride {
                        let k = stride * k + offset;
                        let n = stride * n + offset;
                        x2[k].0 += x1[n].0 * w.0 - x1[n].1 * w.1;
                        x2[k].1 += x1[n].0 * w.1 + x1[n].1 * w.0;
                    }
                }
            }
        }

        // x1 と x2 を入れ替えて次のステップへ
        std::mem::swap(&mut x1, &mut x2);
        stride *= len1;
        len = len2;
    }

    // 逆変換の場合は正規化
    if inverse {
        for n in 0..x.len() {
            x1[n].0 /= x.len() as f64;
            x1[n].1 /= x.len() as f64;
        }
    }

    x1
}

fn main() {
    let len = 2 * 2 * 2 * 3 * 3 * 7 * 13;

    // 1 Hz の正弦波
    let x: Vec<(f64, f64)> = (0..len)
        .map(|n| std::f64::consts::TAU * (n as f64) / (len as f64))
        .map(|n| (n.cos(), n.sin()))
        .collect();

    let y = fft(&x, false);

    for v in y.iter().take(5) {
        println!("{:.3}, {:.3}", v.0 / len as f64, v.1 / len as f64);
        // 0.000, 0.000
        // 1.000, 0.000
        // 0.000, 0.000
        // 0.000, 0.000
        // 0.000, 0.000
        // ...
    }
}
```

参考

- [Cooley-Tukey\_FFT\_algorithm - Wikipedia](https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm#Variations)
- [Cooley-Tukey のアルゴリズム - OTFFT: High Speed FFT Library](http://wwwa.pikara.ne.jp/okojisan/stockham/cooley-tukey.html)

# FFT による巡回畳み込み

配列長 $N$ の配列 $a$ と $b$ があったとき、その巡回畳み込み $c$ は

```math
\begin{aligned}
& c[n] =
  \sum_{m=0}^{N - 1}
    a[m] \cdot b[(n-m) \mod N] &
  & ( 0 \le n < N ) &
\end{aligned}
```

で定義される。  
そして、このときの $a, b, c$ は離散フーリエ変換 $\mathcal{F}$ を用いて以下が成り立つことが知られている。

```math
\mathcal{F}(c) = \mathcal{F}(a) \cdot \mathcal{F}(b)
```

証明:

```math
\begin{aligned}
& \mathcal{F}(c)[k] &=&
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
& &=& \mathcal{F}(a)[k] \cdot \mathcal{F}(b)[k] &
\end{aligned}
```

これに FFT を利用すると巡回畳み込みを $O(N \log N)$ で計算できる。

```math
c[n] = \mathcal{F}^{-1}(\mathcal{F}(a) \cdot \mathcal{F}(b))[n]
```

ただし、Cooley-Tukey FFT では配列長が合成数になっている必要があった。  
そのため、配列長が素数などの場合は適当に長さを調整したい。  
計算結果は変えずに配列長を $N \rightarrow N'$ に伸ばした配列 $a', b'$ を作れるか考える。  
巡回畳み込みの定義をもう一度見返してみる。

```math
\begin{aligned}
& c[n] =
  \sum_{m=0}^{N - 1}
    a[m] \cdot b[(n-m) \mod N] &
  & ( 0 \le n < N ) & \\
\Rightarrow & c[n] =
  \sum_{m=0}^{N' - 1}
    a'[m] \cdot b'[(n-m) \mod N'] &
  & ( 0 \le n < N ) &
\end{aligned}
```


$a'$ は $a$ の末尾を $0$ で埋めると良さそうだが、 $b'$ はやや複雑である。  
$b'$ は $b'[n \mod N']\ (-N < n < N)$ の範囲が元々の計算に利用される。  
言い換えると $b'$ の末尾 $N-1$ 個と先頭 $N$ 個の範囲が利用されるため、この要素は変更せず中央に $0$ を追加することとなる。  
つまり $N'$ は最低でも $2N-1$ 以上にする必要があるということになる。

```math
\begin{aligned}
& a &= [&a_{0}&, &a_{1}&, &a_{2}&] \\
& b &= [&b_{0}&, &b_{1}&, &b_{2}&] \\
&&& \downarrow & \\
& a' &= [&a_{0}&, &a_{1}&, &a_{2}&, &0&, &...& &&&&, &0&] \\
& b' &= [&b_{0}&, &b_{1}&, &b_{2}&, &0&, &...&, &0&, &b_{1}&, &b_{2}&]
\end{aligned}
```

このような配置にすることで任意の配列長の巡回畳み込みを Cooley-Tukey FFT で行うことができる。

```math
c[n] = \mathcal{F}^{-1}(\mathcal{F}(a') \cdot \mathcal{F}(b'))[n]
```

参考

- [Convolution - Wikipedia](https://en.wikipedia.org/wiki/Convolution)

# Bluestein's FFT

Cooley-Tukey FFT では配列長が素数などの場合は計算量を削減することができなかった。  
Bluestein's FFT は Cooley-Tukey FFT より数倍遅いものの、どのような配列長であっても計算量 $O(N \log N)$ でフーリエ変換することができる。

配列長 $N$ の配列 $x$ があったとき、その DFT (離散フーリエ変換) の結果 $X$ は

```math
\begin{aligned}
& X[k] =
  \sum_{n=0}^{N-1} x[n]
  \exp\left(-2\pi i \frac{n k}{N}\right) &
& (0 \le k < N) &
\end{aligned}
```

で定義される。  
ここで $nk$ を

```math
nk =
  \frac{-(k-n)^2}{2} +
  \frac{n^2}{2} + \frac{k^2}{2}
```

のように展開し整理すると

```math
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
```

となる。さらに

```math
\begin{aligned}
& a[n] &=&
  x[n]
  \exp\left(-\pi i \frac{n^2}{N} \right) & \\
& b[n] &=&
  \exp\left(\pi i \frac{n^2}{N} \right) &
\end{aligned}
```

と置き換えると

```math
\begin{aligned}
& X[k] =
  b[k]^* \times
  \left(
    \sum_{n=0}^{N-1}
      a[n]
      b[k-n]
  \right) &
& (b[k]^*\ は\ b[k]\ の複素共役) &
\end{aligned}
```

となり、これは畳み込みと同じ形となる。

つまり、これは「FFT による巡回畳み込み」で紹介した手法と似た要領で計算量 $O(N \log N)$ で解くことができる。  

まずは $a, b$ の長さが $2N-1$ 以上になるよう $0$ 埋めで調節する。  
また $b$ はマイナスインデックスの範囲を末尾に移動する。

```math
\begin{aligned}
& a &= [&&&&&a_{0}&, &a_{1}&, &a_{2}&] \\
& b &= [&b_{-2}&, &b_{-1}&, &b_{0}&, &b_{1}&, &b_{2}&] \\
&&&&&&& \downarrow \\
& a' &= [&&&&&a_{0}&, &a_{1}&, &a_{2}&, &0&, &...&, &0&, &0&, &0&] \\
& b' &= [&b_{-2}&, &b_{-1}&, &b_{0}&, &b_{1}&, &b_{2}&, &0&, &...&, &0&] \\
&&= [&&&&&b_{0}&, &b_{1}&, &b_{2}&, &0&, &...&, &0&, &b_{-2}&, &b_{-1}&] \\
&&= [&&&&&b_{0}&, &b_{1}&, &b_{2}&, &0&, &...&, &0&, &b_{2}&, &b_{1}&]
\end{aligned}
```

このような配置にすると巡回畳み込みと同じ形にできる。

```math
\begin{aligned}
& X[k] =
  b'[k]^* \times
  \left(
    \sum_{n=0}^{N'-1}
      a'[n]
      b'[(k-n) \mod {N'}]
  \right) &
\end{aligned}
```

あとは

```math
X[k] = b'[k]^* \times \mathcal{F}^{-1}(\mathcal{F}(a') \cdot \mathcal{F}(b'))[k]
```

を計算すれば良いこととなる。

最後にサンプルコードを示す。

```
TODO
```

参考

- [Chirp Z-transform - Wikipedia](https://en.wikipedia.org/wiki/Chirp_Z-transform#Bluestein's_algorithm)
