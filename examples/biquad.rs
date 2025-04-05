use cpal::traits::{DeviceTrait, HostTrait, StreamTrait};

fn main() {
    let host = cpal::default_host();
    let device = host.default_output_device().unwrap();
    let config = device.default_output_config().unwrap();
    let rand = std::sync::Arc::new(std::sync::Mutex::new(math::rand::XorShift32::default()));
    let biquad = std::sync::Arc::new(std::sync::Mutex::new(math::biquad::BiQuad::new(
        config.sample_rate().0 as f64,
    )));
    println!("{}", config.sample_rate().0);
    let count = std::sync::Arc::new(std::sync::Mutex::new(0u64));
    let buffer = std::sync::Arc::new(std::sync::Mutex::new(vec![0.0; 1024]));
    let stream = device
        .build_output_stream(
            &config.into(),
            move |data: &mut [f32], _: &cpal::OutputCallbackInfo| {
                for sample in data.iter_mut() {
                    let mut s = (rand.lock().unwrap().f64() * 2.0 - 1.0) as f32;
                    let count = {
                        let mut count = count.lock().unwrap();
                        *count += 1;
                        *count as f64
                    };
                    let cutoff = 660.0 + (count / 44100.0).sin() * 220.0;
                    biquad.lock().unwrap().band_pass(cutoff, 0.1);
                    //biquad.lock().unwrap().low_pass(cutoff, 0.5f64.exp2());
                    //biquad.lock().unwrap().high_pass(cutoff, 0.5f64.exp2());
                    s = biquad.lock().unwrap().process(s as f64) as f32;
                    *sample = s
                }
                let mut buf = buffer.lock().unwrap();
                if buf.len() > data.len() {
                    buf.rotate_left(data.len());
                    buf[data.len()..].copy_from_slice(data);
                } else {
                    let buflen = buf.len();
                    buf.copy_from_slice(&data[data.len() - buflen..]);
                }
                let buf = math::fft::window(
                    &buf.iter()
                        .map(|x| math::complex::Complex::new(*x as f64 * 0.05, 0.0))
                        .collect(),
                    false,
                );
                let fft = math::fft::fft_range(&buf, 0.0, 1320.0 * 1024.0 / 44100.0, 80);
                std::hint::black_box(fft);

                let fft = biquad
                    .lock()
                    .unwrap()
                    .freq_response(0.0, 1320.0, 80)
                    .iter()
                    .map(|x| math::complex::Complex::new(*x * 0.1, 0.0))
                    .collect::<Vec<_>>();
                math::plot::plot(&fft, 80, 20, false);
            },
            move |err| eprintln!("エラーが発生しました: {}", err),
            None,
        )
        .unwrap();
    stream.play().unwrap();
    std::thread::park();
}
