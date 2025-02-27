FFT やフィルターなどのスクラッチ実装

やること

- FFT
    - [x] グラフ生成環境の用意
    - [x] 愚直な DFT / iDFT の実装
    - [x] 2 基 Cooley-Tukey FFT の実装
    - [x] 2 基 Stockham FFT の理解・実装
    - [ ] テストの作成
    - [ ] ベンチマークコードの作成
    - [ ] 畳み込みの実装
    - [ ] N 基 Stockham FFT の理解・実装
    - [ ] Bluestein's FFT の理解・実装
- 疑似乱数生成
    - [x] MT19937 の実装
    - [x] f64 の疑似乱数生成
- 双2次フィルタ
    - [ ] ハイパス, ローパス, バンドパスフィルターの実装
- JavaScript への移植

参考にしているコード/サイト

- FFT
    - http://wwwa.pikara.ne.jp/okojisan/stockham/cooley-tukey.html
    - https://github.com/mreineck/pocketfft
    - https://github.com/ejmahler/RustFFT
- MT19937
    - https://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/mt.html
    - https://cpprefjp.github.io/reference/random/mt19937.html
