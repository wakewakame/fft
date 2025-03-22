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

## 具体例

$X$ を偶数個目と奇数個目に分けて離散フーリエ変換する例を考えてみる。  
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

## サンプルコード

最後にサンプルコードを示す。  
このコードは 0BSD ライセンス (出典明記不要で自由に利用可能) の元に公開する。

[playground](https://play.rust-lang.org/?version=stable&mode=release&edition=2024&gist=bde7a1d497e134169c0188b2f5d7d5f9)

```rust
fn main() {
    // DFT と Cooley-Tukey FFT の計算時間を比較
    println!("len\tdft (sec)\tfft_cooley_tukey (sec)");
    let mut rand = random();
    let len_patterns: Vec<usize> = (0..=10).map(|x| 1 << x).collect();
    for len in len_patterns {
        let x: Vec<(f64, f64)> = (0..len).map(|_| (rand(), rand())).collect();

        let start = std::time::Instant::now();
        let y1 = dft(&x, false);
        let dft_dur = start.elapsed().as_secs_f64();

        let start = std::time::Instant::now();
        let y2 = fft_cooley_tukey(&x, false);
        let cooley_tukey_dur = start.elapsed().as_secs_f64();

        std::hint::black_box((y1, y2));
        println!("{}\t{:.9}\t{:.9}", len, dft_dur, cooley_tukey_dur);
    }
}

fn dft(x: &Vec<(f64, f64)>, inverse: bool) -> Vec<(f64, f64)> {
    // 事前に w を計算しておく
    let sign = if inverse { 1.0 } else { -1.0 };
    let w: Vec<(f64, f64)> = (0..x.len())
        .map(|n| {
            let theta = sign * std::f64::consts::TAU * n as f64 / x.len() as f64;
            (theta.cos(), theta.sin())
        })
        .collect();

    // DFT
    let mut y = vec![(0f64, 0f64); x.len()];
    for k in 0..y.len() {
        for n in 0..x.len() {
            let w_nk = w[(n * k) % w.len()];
            y[k].0 += x[n].0 * w_nk.0 - x[n].1 * w_nk.1;
            y[k].1 += x[n].0 * w_nk.1 + x[n].1 * w_nk.0;
        }
        if inverse {
            y[k].0 /= y.len() as f64;
            y[k].1 /= y.len() as f64;
        }
    }

    return y;
}

fn fft_cooley_tukey(x: &Vec<(f64, f64)>, inverse: bool) -> Vec<(f64, f64)> {
    // 事前に w を計算しておく
    let sign = if inverse { 1.0 } else { -1.0 };
    let w: Vec<(f64, f64)> = (0..x.len())
        .map(|n| {
            let theta = sign * std::f64::consts::TAU * n as f64 / x.len() as f64;
            (theta.cos(), theta.sin())
        })
        .collect();

    // 出力用 (y1 と y2 を交互に使い回す)
    let mut y1 = x.clone();
    let mut y2 = vec![(0f64, 0f64); x.len()];

    // FFT
    let mut len = y1.len();
    let mut stride = 1;
    while len >= 2 {
        // len を 2 つの整数の積に分解
        let len1 = (2..len).find(|&n| len % n == 0).unwrap_or(len);
        let len2 = len / len1;

        // バタフライ演算
        y2.iter_mut().for_each(|y| *y = (0.0, 0.0));
        for k1 in 0..len1 {
            for n1 in 0..len1 {
                for n2 in 0..len2 {
                    let k = len1 * n2 + k1;
                    let n = len2 * n1 + n2;
                    let w = w[(stride * n * k1) % w.len()];
                    for offset in 0..stride {
                        let k = stride * k + offset;
                        let n = stride * n + offset;
                        y2[k].0 += y1[n].0 * w.0 - y1[n].1 * w.1;
                        y2[k].1 += y1[n].0 * w.1 + y1[n].1 * w.0;
                    }
                }
            }
        }

        // y1 と y2 を入れ替えて次のステップへ
        std::mem::swap(&mut y1, &mut y2);
        len = len2;
        stride *= len1;
    }

    // 逆変換の場合は正規化
    if inverse {
        for k in 0..y1.len() {
            y1[k].0 /= y1.len() as f64;
            y1[k].1 /= y1.len() as f64;
        }
    }

    return y1;
}

// 疑似乱数生成 (xorshift32)
fn random() -> impl FnMut() -> f64 {
    let mut state = 1u32;
    move || {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        return (state - 1) as f64 / std::u32::MAX as f64;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dft() {
        // 1 Hz の正弦波を用意
        let x: Vec<(f64, f64)> = (0..8)
            .map(|n| std::f64::consts::TAU * n as f64 / 8.0)
            .map(|n| (n.cos(), n.sin()))
            .collect();

        // 期待する DFT の結果
        let expect_fft: Vec<(f64, f64)> = vec![0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            .into_iter()
            .map(|n| (n, 0.0))
            .collect();

        // DFT
        let actual_fft = dft(&x, false);
        assert_eq!(actual_fft.len(), expect_fft.len());
        for (a, b) in expect_fft.iter().zip(actual_fft.iter()) {
            assert!((a.0 - b.0).abs() < 1e-10, "{} != {}", a.0, b.0);
            assert!((a.1 - b.1).abs() < 1e-10, "{} != {}", a.1, b.1);
        }

        // iDFT
        let actual_ifft = dft(&actual_fft, true);
        assert_eq!(actual_ifft.len(), x.len());
        for (a, b) in x.iter().zip(actual_ifft.iter()) {
            assert!((a.0 - b.0).abs() < 1e-10, "{} != {}", a.0, b.0);
            assert!((a.1 - b.1).abs() < 1e-10, "{} != {}", a.1, b.1);
        }
    }

    #[test]
    fn test_fft() {
        // 乱数列を用意
        let mut rand = random();
        let len = 2 * 3 * 5 * 7;
        let x: Vec<(f64, f64)> = (0..len).map(|_| (rand(), rand())).collect();

        // FFT の結果が DFT と一致することを確認する
        let expect_fft = dft(&x, false);

        // FFT
        let actual_fft = fft_cooley_tukey(&x, false);
        assert_eq!(actual_fft.len(), expect_fft.len());
        for (a, b) in expect_fft.iter().zip(actual_fft.iter()) {
            assert!((a.0 - b.0).abs() < 1e-10, "{} != {}", a.0, b.0);
            assert!((a.1 - b.1).abs() < 1e-10, "{} != {}", a.1, b.1);
        }

        // iFFT
        let actual_ifft = fft_cooley_tukey(&actual_fft, true);
        assert_eq!(actual_ifft.len(), x.len());
        for (a, b) in x.iter().zip(actual_ifft.iter()) {
            assert!((a.0 - b.0).abs() < 1e-10, "{} != {}", a.0, b.0);
            assert!((a.1 - b.1).abs() < 1e-10, "{} != {}", a.1, b.1);
        }
    }

    #[test]
    fn test_random() {
        let mut rand = random();
        for _ in 0..9999 {
            rand();
        }
        assert_eq!(rand(), 1799336687.0 / std::u32::MAX as f64);
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

## サンプルコード

最後にサンプルコードを示す。  
このコードは 0BSD ライセンス (出典明記不要で自由に利用可能) の元に公開する。

[playground](https://play.rust-lang.org/?version=stable&mode=release&edition=2024&gist=88b5e5c692e0ebff8514f910863dd3e4)

```rust
fn main() {
    // Cooley-Tukey FFT と Bluestein's FFT の計算時間を比較
    println!("len\tfft_cooley_tukey (sec)\tfft_bluestein (sec)");
    let mut rand = random();
    let len_patterns: Vec<usize> = (0..=10).map(|x| 1 << x).collect();
    for len in len_patterns {
        let x: Vec<(f64, f64)> = (0..len).map(|_| (rand(), rand())).collect();

        let start = std::time::Instant::now();
        let y1 = fft_cooley_tukey(&x, false);
        let cooley_tukey_dur = start.elapsed().as_secs_f64();

        let start = std::time::Instant::now();
        let y2 = fft_bluestein(&x, false);
        let bluestein_dur = start.elapsed().as_secs_f64();

        std::hint::black_box((y1, y2));
        println!("{}\t{:.9}\t{:.9}", len, cooley_tukey_dur, bluestein_dur);
    }
}

fn fft_bluestein(x: &Vec<(f64, f64)>, inverse: bool) -> Vec<(f64, f64)> {
    // 事前に w を計算しておく
    let sign = if inverse { 1.0 } else { -1.0 };
    let w: Vec<(f64, f64)> = (0..x.len())
        .map(|n| {
            let theta = sign * std::f64::consts::PI * (n * n) as f64 / x.len() as f64;
            (theta.cos(), theta.sin())
        })
        .collect();

    // Bluestein's FFT
    let len = (x.len() * 2 - 1).next_power_of_two();
    let a: Vec<(f64, f64)> = (0..len)
        .map(|n| match n {
            _ if n < x.len() => (
                x[n].0 * w[n].0 - x[n].1 * w[n].1,
                x[n].0 * w[n].1 + x[n].1 * w[n].0,
            ),
            _ => (0.0, 0.0),
        })
        .collect();
    let b: Vec<(f64, f64)> = (0..len)
        .map(|n| match n {
            _ if n < x.len() => (w[n].0, -w[n].1),
            _ if len - n < x.len() => (w[len - n].0, -w[len - n].1),
            _ => (0.0, 0.0),
        })
        .collect();
    let a_fft = fft_cooley_tukey(&a, false);
    let b_fft = fft_cooley_tukey(&b, false);
    let y_fft: Vec<(f64, f64)> = (a_fft.iter().zip(b_fft.iter()))
        .map(|(a, b)| (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0))
        .collect();
    let mut y = fft_cooley_tukey(&y_fft, true);
    y.truncate(x.len());

    for k in 0..y.len() {
        y[k] = (
            y[k].0 * b[k].0 - y[k].1 * -b[k].1,
            y[k].0 * -b[k].1 + y[k].1 * b[k].0,
        );
        if inverse {
            y[k].0 /= y.len() as f64;
            y[k].1 /= y.len() as f64;
        }
    }

    return y;
}

fn fft_cooley_tukey(x: &Vec<(f64, f64)>, inverse: bool) -> Vec<(f64, f64)> {
    // 事前に w を計算しておく
    let sign = if inverse { 1.0 } else { -1.0 };
    let w: Vec<(f64, f64)> = (0..x.len())
        .map(|n| {
            let theta = sign * std::f64::consts::TAU * n as f64 / x.len() as f64;
            (theta.cos(), theta.sin())
        })
        .collect();

    // 出力用 (y1 と y2 を交互に使い回す)
    let mut y1 = x.clone();
    let mut y2 = vec![(0f64, 0f64); x.len()];

    // FFT
    let mut len = y1.len();
    let mut stride = 1;
    while len >= 2 {
        // len を 2 つの整数の積に分解
        let len1 = (2..len).find(|&n| len % n == 0).unwrap_or(len);
        let len2 = len / len1;

        // バタフライ演算
        y2.iter_mut().for_each(|y| *y = (0.0, 0.0));
        for k1 in 0..len1 {
            for n1 in 0..len1 {
                for n2 in 0..len2 {
                    let k = len1 * n2 + k1;
                    let n = len2 * n1 + n2;
                    let w = w[(stride * n * k1) % w.len()];
                    for offset in 0..stride {
                        let k = stride * k + offset;
                        let n = stride * n + offset;
                        y2[k].0 += y1[n].0 * w.0 - y1[n].1 * w.1;
                        y2[k].1 += y1[n].0 * w.1 + y1[n].1 * w.0;
                    }
                }
            }
        }

        // y1 と y2 を入れ替えて次のステップへ
        std::mem::swap(&mut y1, &mut y2);
        len = len2;
        stride *= len1;
    }

    // 逆変換の場合は正規化
    if inverse {
        for k in 0..y1.len() {
            y1[k].0 /= y1.len() as f64;
            y1[k].1 /= y1.len() as f64;
        }
    }

    return y1;
}

// 疑似乱数生成 (xorshift32)
fn random() -> impl FnMut() -> f64 {
    let mut state = 1u32;
    move || {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        return (state - 1) as f64 / std::u32::MAX as f64;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fft_bluestein() {
        // 乱数列を用意
        let mut rand = random();
        let len = 2 * 3 * 5 * 7;
        let x: Vec<(f64, f64)> = (0..len).map(|_| (rand(), rand())).collect();

        // Bluestein's FFT の結果が Cooley-Tukey FFT と一致することを確認する
        let expect_fft = fft_cooley_tukey(&x, false);

        // FFT
        let actual_fft = fft_bluestein(&x, false);
        assert_eq!(actual_fft.len(), expect_fft.len());
        for (a, b) in expect_fft.iter().zip(actual_fft.iter()) {
            assert!((a.0 - b.0).abs() < 1e-10, "{} != {}", a.0, b.0);
            assert!((a.1 - b.1).abs() < 1e-10, "{} != {}", a.1, b.1);
        }

        // iFFT
        let actual_ifft = fft_bluestein(&actual_fft, true);
        assert_eq!(actual_ifft.len(), x.len());
        for (a, b) in x.iter().zip(actual_ifft.iter()) {
            assert!((a.0 - b.0).abs() < 1e-10, "{} != {}", a.0, b.0);
            assert!((a.1 - b.1).abs() < 1e-10, "{} != {}", a.1, b.1);
        }
    }
}
```

参考

- [Chirp Z-transform - Wikipedia](https://en.wikipedia.org/wiki/Chirp_Z-transform#Bluestein's_algorithm)
