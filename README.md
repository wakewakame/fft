FFT やフィルターなどのスクラッチ実装

やること

- FFT
    - [x] グラフ生成環境の用意
    - [x] 愚直な DFT / iDFT の実装
    - [x] 2 基 Cooley-Tukey FFT の実装
    - [x] 2 基 Stockham FFT の理解・実装
    - [x] テストの作成
    - [x] 畳み込みの実装
    - [x] Bluestein's FFT の理解・実装
    - [x] N 基 Stockham FFT の理解・実装
    - [x] ベンチマークコードの作成
    - [ ] 要素数によってアルゴリズムを切り替える
- 疑似乱数生成
    - [x] xorshift32 の実装
    - [x] f64 の疑似乱数生成
- 双2次フィルタ
    - [ ] ハイパス, ローパス, バンドパスフィルターの実装
- 0BSD でライセンスする
- 備忘録記事を書く
- JavaScript への移植

参考にしているコード/サイト

- FFT
    - http://wwwa.pikara.ne.jp/okojisan/stockham/cooley-tukey.html
    - https://github.com/mreineck/pocketfft
    - https://github.com/ejmahler/RustFFT
