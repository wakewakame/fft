FFT やフィルターなどのスクラッチ実装

サンプル実行

```sh
cargo run --release --example bench
```

提供する機能

- [x] 任意の要素数の FFT
- [x] 疑似乱数生成
- [x] 畳み込み
- [x] 双2次フィルタ
- [ ] 音程推定
- [ ] 重畳加算法
- [ ] リバーブ
- [ ] ピッチ変換 (Phase Vocoder か WSOLA のどちらか)
- [ ] JavaScript への移植
