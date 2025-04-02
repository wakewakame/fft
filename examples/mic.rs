use cpal::traits::{DeviceTrait, HostTrait, StreamTrait};
use math::complex::*;
use math::fft::*;
use math::plot::*;
use std::sync::mpsc;

fn main() {
    // マイク入力の取得
    let (tx, rx) = mpsc::channel();
    let host = cpal::default_host();
    let device = host.default_input_device().unwrap();
    let config = device.default_input_config().unwrap();
    let stream = device
        .build_input_stream(
            &config.into(),
            move |data: &[f32], _: &cpal::InputCallbackInfo| {
                let samples: Vec<f64> = data.iter().map(|&x| x as f64).collect();
                let _ = tx.send(samples);
            },
            move |err| panic!("error: {}", err),
            None,
        )
        .unwrap();
    stream.play().unwrap();

    // マイク入力の表示
    let mut buffer = ([0f64; 1024], 0usize);
    for samples in rx {
        let mut samples = samples;
        for _ in 0..(1 + samples.len() / buffer.0.len()) {
            if buffer.1 + samples.len() < buffer.0.len() {
                // buffer が埋まるまで待つ
                buffer.0[buffer.1..buffer.1 + samples.len()].copy_from_slice(&samples);
                buffer.1 += samples.len();
            } else {
                let (first, second) = samples.split_at(buffer.0.len() - buffer.1);
                buffer.0[buffer.1..].copy_from_slice(first);
                samples = second.to_vec();
                buffer.1 = 0;

                // 描画
                let fft = fft(&to_complex(&buffer.0.to_vec()), true);
                let wave = fft.iter().map(|&x| x * 100.0).collect();
                plot(&wave, 40, 10, false);
                plot(
                    &buffer
                        .0
                        .iter()
                        .map(|&x| Complex::new(0.5 + 20.0 * x, 0.0))
                        .collect(),
                    40,
                    10,
                    false,
                );
            }
        }
    }
}
