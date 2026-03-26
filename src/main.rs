//! Fluid simulation — live window or video export.
//!
//! Usage:
//!   cargo run --release                        # live window (ESC to quit)
//!   cargo run --release -- --video out.mp4     # render 300 frames to video

use fluid_sim::FluidGrid;
use minifb::{Key, Window, WindowOptions};

// ── Simulation parameters ─────────────────────────────────────────────────────

const GRID_W: u32 = 120;
const GRID_H: u32 = 80;
const CELL_SIZE: f32 = 1.0;
/// Pixels per cell.
const SCALE: usize = 8;
/// Number of Gauss-Seidel pressure iterations per step.
const PRESSURE_ITERS: usize = 80;
/// Horizontal inflow speed on the left wall.
const INFLOW_SPEED: f32 = 15.0;
/// Vertical band (fraction of height) that emits smoke.
const SMOKE_LO: f32 = 0.35;
const SMOKE_HI: f32 = 0.65;
/// Frames to render in video mode.
const VIDEO_FRAMES: usize = 6000;
/// FPS for video mode.
const VIDEO_FPS: usize = 200;

const WIN_W: usize = GRID_W as usize * SCALE;
const WIN_H: usize = GRID_H as usize * SCALE;

// ── Entry point ───────────────────────────────────────────────────────────────

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let video_out = args
        .windows(2)
        .find(|w| w[0] == "--video")
        .map(|w| w[1].clone());

    let mut grid = FluidGrid::new((GRID_W, GRID_H), CELL_SIZE);

    if let Some(path) = video_out {
        run_video(&mut grid, &path);
    } else {
        run_live(&mut grid);
    }
}

// ── Inflow injection (called every step) ─────────────────────────────────────

fn inject_inflow(grid: &mut FluidGrid) {
    let h = grid.cell_count.1;
    let lo = (h as f32 * SMOKE_LO) as u32;
    let hi = (h as f32 * SMOKE_HI) as u32;
    for y in 0..h {
        grid.velocities_x[(0, y)] = INFLOW_SPEED;
        grid.smoke[(0, y)] = if y >= lo && y < hi { 1.0 } else { 0.0 };
    }
}

// ── Rendering ─────────────────────────────────────────────────────────────────
/// Write every cell's colour into `buffer` (0x00RRGGBB, row-major, y=0 at top).
fn render(grid: &FluidGrid, buffer: &mut Vec<u32>) {
    let (w, h) = (grid.cell_count.0 as usize, grid.cell_count.1 as usize);

    for x in 0..w {
        for y in 0..h {
            let smoke = grid.smoke[(x as u32, y as u32)].clamp(0.0, 1.0);
            let (vx, vy) = grid.velocity_at_center(x as u32, y as u32);
            let speed = (vx * vx + vy * vy).sqrt();
            let color = cell_color(smoke, speed);

            // Flip y so that y=0 appears at the bottom of the screen.
            let screen_y = (h - 1 - y) * SCALE;
            let screen_x = x * SCALE;
            for dy in 0..SCALE {
                for dx in 0..SCALE {
                    buffer[(screen_y + dy) * WIN_W + screen_x + dx] = color;
                }
            }
        }
    }
}

/// Smoke density + speed → colour.
///
/// Slow smoke → deep blue-white.
/// Fast smoke → warm orange-white.
/// No smoke   → near-black background with a faint blue tint from velocity.
fn cell_color(density: f32, speed: f32) -> u32 {
    let max_speed = INFLOW_SPEED * 1.5;
    let t = (speed / max_speed).clamp(0.0, 1.0); // 0 = slow, 1 = fast
    let d = density;

    // Smoke colour: lerp blue-white → orange-white based on local speed.
    let sr = lerp(0.35, 1.00, t) * d;
    let sg = lerp(0.55, 0.75, t) * d;
    let sb = lerp(1.00, 0.25, t) * d;

    // Background: faint blue glow proportional to speed (shows flow even without smoke).
    let bg = (t * 0.12).min(1.0);
    let br = bg * 0.1;
    let bg_ = bg * 0.2;
    let bb = bg * 0.6;

    let r = ((sr + br).clamp(0.0, 1.0) * 255.0) as u32;
    let g = ((sg + bg_).clamp(0.0, 1.0) * 255.0) as u32;
    let b = ((sb + bb).clamp(0.0, 1.0) * 255.0) as u32;
    (r << 16) | (g << 8) | b
}

#[inline]
fn lerp(a: f32, b: f32, t: f32) -> f32 {
    a + (b - a) * t
}

// ── Live window mode ──────────────────────────────────────────────────────────

fn run_live(grid: &mut FluidGrid) {
    let mut window = Window::new(
        "Fluid Simulation  |  ESC to quit",
        WIN_W,
        WIN_H,
        WindowOptions::default(),
    )
    .expect("Could not open window");

    window.set_target_fps(60);
    let mut buffer = vec![0u32; WIN_W * WIN_H];
    let mut frame = 0u64;

    while window.is_open() && !window.is_key_down(Key::Escape) {
        inject_inflow(grid);
        grid.step(PRESSURE_ITERS);

        render(grid, &mut buffer);
        window.update_with_buffer(&buffer, WIN_W, WIN_H).unwrap();

        if frame % 60 == 0 {
            println!(
                "frame {} | max divergence: {:.4}",
                frame,
                grid.max_divergence()
            );
        }
        frame += 1;
    }
}

// ── Video export mode ─────────────────────────────────────────────────────────

fn run_video(grid: &mut FluidGrid, output_path: &str) {
    use std::io::Write;
    use std::process::{Command, Stdio};

    println!("Rendering {VIDEO_FRAMES} frames → {output_path}");
    println!("(requires ffmpeg in PATH)");

    let mut ffmpeg = Command::new("ffmpeg")
        .args([
            "-y",
            "-f",
            "rawvideo",
            "-pixel_format",
            "rgb24",
            "-video_size",
            &format!("{WIN_W}x{WIN_H}"),
            "-framerate",
            &VIDEO_FPS.to_string(),
            "-i",
            "pipe:0",
            "-c:v",
            "libx264",
            "-pix_fmt",
            "yuv420p",
            "-crf",
            "18",
            output_path,
        ])
        .stdin(Stdio::piped())
        .spawn()
        .expect("Failed to start ffmpeg — is it installed?");

    let stdin = ffmpeg.stdin.as_mut().unwrap();
    let mut buffer = vec![0u32; WIN_W * WIN_H];

    for i in 0..VIDEO_FRAMES {
        inject_inflow(grid);
        grid.step(PRESSURE_ITERS);
        render(grid, &mut buffer);

        // Convert 0x00RRGGBB → [R, G, B] bytes.
        let bytes: Vec<u8> = buffer
            .iter()
            .flat_map(|&p| {
                [
                    ((p >> 16) & 0xFF) as u8, // R
                    ((p >> 8) & 0xFF) as u8,  // G
                    (p & 0xFF) as u8,         // B
                ]
            })
            .collect();

        stdin.write_all(&bytes).expect("Failed to write frame");

        if (i + 1) % VIDEO_FPS == 0 {
            println!(
                "  {}/{VIDEO_FRAMES} frames  (div={:.4})",
                i + 1,
                grid.max_divergence()
            );
        }
    }

    drop(ffmpeg.stdin.take()); // close stdin → signals EOF to ffmpeg
    ffmpeg.wait().expect("ffmpeg did not finish successfully");
    println!("Done → {output_path}");
}
