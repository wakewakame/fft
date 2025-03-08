# Cooley-Tukey FFT

DFT (離散フーリエ変換) の定義は以下である。

$$
\begin{aligned}
& X\left[k\right] =
  \sum_{n=0}^{N-1} x\left[n\right]
  \exp\left(-2\pi i \frac{n k}{N}\right) &
& \left(0 \le k < N\right) &
\end{aligned}
$$

ここで $x$ の要素数 $N$ が 2 つの整数の積 $N = N_1 \times N_2$ で表せる時、 DFT は 2 つの DFT に分割して計算することができる。

上式に

$$
\begin{aligned}
& N &=& N_1 \times N_2 & & & \\
& k &=& N_1 k_2 + k_1 &
& \left( 0 \le k_1 < N_1, 0 \le k_2 < N_2 \right) &
\end{aligned}
$$

を代入して整理すると

$$
\begin{aligned}
& X\left[k\right] &=&
  \sum_{n=0}^{N-1}
    x\left[n\right]
    \exp\left(-2\pi i \frac{n k}{N} \right) & \\
& &=&
  \sum_{n=0}^{N_1 N_2 - 1}
    x\left[n\right]
    \exp\left(-2\pi i \frac{n k}{N_1 N_2} \right) & \\
& &=&
  \sum_{n_2=0}^{N_2-1}
    \sum_{n_1=0}^{N_1-1}
      x\left[N_2 n_1 + n_2\right]
      \exp\left(-2\pi i \frac{(N_2 n_1 + n_2) k}{N_1 N_2} \right) & \\
\Leftrightarrow & X\left[N_1 k_2 + k_1\right] &=&
  \sum_{n_2=0}^{N_2-1}
    \sum_{n_1=0}^{N_1-1}
      x\left[N_2 n_1 + n_2\right]
      \exp\left(-2\pi i \frac{(N_2 n_1 + n_2) (N_1 k_2 + k_1)}{N_1 N_2} \right) & \\
& &=&
  \sum_{n_2=0}^{N_2-1}
    \sum_{n_1=0}^{N_1-1}
      \left(
        x\left[N_2 n_1 + n_2\right]
        \exp\left(-2\pi i \frac{n_2 k_1}{N_1 N_2} \right)
      \right)
      \exp\left(-2\pi i \frac{n_1 k_1}{N_1} \right)
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right)
    \exp\left(-2\pi i n_1 k_2 \right) & \\
&&& \downarrow \textrm{末尾の}\ \exp\left(-2\pi i n_1 k_2 \right)\ \textrm{は常に}\ 1\ \textrm{となるので削除できる} & \\
& &=&
  \sum_{n_2=0}^{N_2-1}
    \sum_{n_1=0}^{N_1-1}
      \left(
        x\left[N_2 n_1 + n_2\right]
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
& X\left[2 k_2 + 0\right] &=&
  \sum_{n_2=0}^{N_2-1}
    \sum_{n_1=0}^{2-1}
      \left(
        x\left[N_2 n_1 + n_2\right]
        \exp\left(-2\pi i \frac{0}{2 N_2}\right)
      \right)
      \exp\left(-2\pi i \frac{0}{N_1} \right)
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) & \\
& &=&
  \sum_{n_2=0}^{N_2-1}
    \left(
      \sum_{n_1=0}^{2-1}
        x\left[N_2 n_1 + n_2\right]
    \right)
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) & \\
& &=&
  \sum_{n_2=0}^{N_2-1}
    \left(x\left[n_2\right] + x\left[N_2 + n_2\right]\right)
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) & \\
& X\left[2 k_2 + 1\right] &=&
  \sum_{n_2=0}^{N_2-1}
    \sum_{n_1=0}^{2-1}
      \left(
        x\left[N_2 n_1 + n_2\right]
        \exp\left(-2\pi i \frac{n_2}{2 N_2}\right)
      \right)
      \exp\left(-2\pi i \frac{n_1}{2} \right)
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) & \\
& &=&
  \sum_{n_2=0}^{N_2-1}
    \left(
      \sum_{n_1=0}^{2-1}
        \left(
          x\left[N_2 n_1 + n_2\right]
        \right)
        \exp\left(-2\pi i \frac{n_1}{2} \right)
    \right)
    \exp\left(-2\pi i \frac{n_2}{2 N_2}\right)
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) & \\
& &=&
  \sum_{n_2=0}^{N_2-1}
    \left(
      x\left[n_2\right] \exp\left( 0 \right) +
      x\left[N_2 + n_2\right] \exp\left( -\pi i \right)
    \right)
    \exp\left(-2\pi i \frac{n_2}{2 N_2}\right)
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) & \\
& &=&
  \sum_{n_2=0}^{N_2-1}
    \left(
      x\left[n_2\right] -
      x\left[N_2 + n_2\right]
    \right)
    \exp\left(-2\pi i \frac{n_2}{2 N_2}\right)
    \exp\left(-2\pi i \frac{n_2 k_2}{N_2} \right) & \\
\end{aligned}
$$

となり、要素数 $N$ の DFT が要素数 $N/2$ の DFT 2 回分に分割することができた。  
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

# FFT による畳み込み

畳み込みの定義は以下である。

$$
\begin{aligned}
& c[m] =
  \sum_{n=0}^{N_a - 1}
    a[n] \times b[m-n] &
  & \left( 0 \le m < N_a + N_b - 1\right) &
\end{aligned}
$$

ここで $N_a, N_b$ はそれぞれ $a, b$ の要素数である。  
また $b[m - n]$ の $m-n$ が負になったり配列の範囲外を指すことがあるが、その場合は $0$ として扱う。

$$
...\textsf{TODO}
$$

この式はコードで説明した方がわかりやすいかもしれない。

```rust
fn convolution(a: &Vec<f64>, b: &Vec<f64>) -> Vec<f64> {
  if a.len() == 0 && b.len() == 0 { return vec![]; }
  let mut c = vec![0f64; a.len() + b.len() - 1];
  for ai in 0..a.len() {
    for bi in 0..b.len() {
      c[ai + bi] += a[ai] * b[bi];
    }
  }
  return c;
}

fn main() {
  let a: Vec<f64> = vec![1.0, 2.0, 3.0];
  let b: Vec<f64> = vec![4.0, 5.0];
  let c = convolution(&a, &b);
  println!("{:?}", c);  // [4.0, 13.0, 22.0, 15.0]
}
```

参考

- [Convolution - Wikipedia](https://en.wikipedia.org/wiki/Convolution)

# Bluestein's FFT

Cooley-Tukey FFT では要素数が $N = N_1 \times N_2$ のように合成数になる配列しか計算量を削減することができなかった。  
Bluestein's FFT は Cooley-Tukey FFT より数倍遅いものの、要素数が素数など任意の配列長を計算量 $O(N \log N)$ でフーリエ変換することができる。

DFT (離散フーリエ変換) の定義は以下である。

$$
\begin{aligned}
& X\left[k\right] =
  \sum_{n=0}^{N-1} x\left[n\right]
  \exp\left(-2\pi i \frac{n k}{N}\right) &
& \left(0 \le k < N\right) &
\end{aligned}
$$

ここで $nk$ を

$$
nk =
  \frac{-(k-n)^2}{2} +
  \frac{n^2}{2} + \frac{k^2}{2}
$$

のように展開し整理すると

$$
\begin{aligned}
& X\left[k\right] &=&
  \sum_{n=0}^{N-1}
    x\left[n\right]
    \exp\left(-2\pi i \frac{n k}{N} \right) & \\
& &=&
  \exp\left(-\pi i \frac{k^2}{N} \right)
  \sum_{n=0}^{N-1}
    x\left[n\right]
    \exp\left(-\pi i \frac{n^2}{N} \right)
  \exp\left(\pi i \frac{(k-n)^2}{N} \right) & \\
\end{aligned}
$$

となる。さらに

$$
\begin{aligned}
& a\left[n\right] &=&
  x\left[n\right]
  \exp\left(-\pi i \frac{n^2}{N} \right) & \\
& b\left[n\right] &=&
  \exp\left(\pi i \frac{n^2}{N} \right) &
\end{aligned}
$$

と置き換えると

$$
\begin{aligned}
& X\left[k\right] =
  b\left[k\right]^*
  \times
  \left(
    \sum_{n=0}^{N-1}
      a\left[n\right]
      b\left[k-n\right]
  \right) &
& \left(b\left[k\right]^*\ \textrm{は}\ b\left[k\right]\ \textrm{の複素共役} \right) &
\end{aligned}
$$

となり、巡回畳み込みと同じ形となる。  
よって、これは Cooley-Tukey FFT を使った畳み込み手法を用いて計算量 $O(N \log N)$ で計算することができる。

**TODO: 以下の文章は書きかけ。もう少しわかりやすく直したい。**

...だた、これだと $a, b$ の配列サイズが合成数ではないケースもあるので、そのままだと Cooley-Tukey FFT を使った畳み込みができない。
そのため、配列の数が $2^N$ などになるよう、配列末尾に $0$ を追加して配列サイズを稼ぐ。

...また、これだと $b\left[k-n\right]$ は $-(N-1) \le k-n < N$ の範囲で利用されるので、マイナス範囲へのアクセスがコードに落とし込みづらい。  
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
& X\left[k\right] =
  b\left[k\right]^*
  \times
  \left(
    \sum_{n=0}^{N-1}
      a\left[n\right]
      b'\left[\left(k-n\right) \mod N_b\right]
  \right) &
& \left(b\left[k\right]^* は b\left[k\right] の複素共役\right) &
\end{aligned}
$$

最後にサンプルコードを示す。

```
TODO
```

参考

- [Chirp Z-transform - Wikipedia](https://en.wikipedia.org/wiki/Chirp_Z-transform#Bluestein's_algorithm)
