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
const INFLOW_SPEED: f32 = 5.0;
/// Frames to render in video mode.
const VIDEO_FRAMES: usize = 3000;
/// FPS for video mode.
const VIDEO_FPS: usize = 60;

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
    }
}

// ── Inflow injection (called every step) ─────────────────────────────────────

fn inject_inflow(grid: &mut FluidGrid) {
    let h = grid.cell_count.1;
    for y in 0..h {
        grid.velocities_x[(0, y)] = INFLOW_SPEED;
    }
}

// ── Rendering ─────────────────────────────────────────────────────────────────
/// Write every cell's colour into `buffer` (0x00RRGGBB, row-major, y=0 at top).
fn render(grid: &FluidGrid, buffer: &mut Vec<u32>) {
    let (w, h) = (grid.cell_count.0 as usize, grid.cell_count.1 as usize);

    for x in 0..w {
        for y in 0..h {
            let (vx, vy) = grid.velocity_at_center(x as u32, y as u32);
            let color = velocity_color(vx, vy);

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

/// Visualisation des vitesses :
/// - Rouge = mouvement vers la droite
/// - Bleu = mouvement vers la gauche
/// - Vert = mouvement vers le haut
/// - Jaune = mouvement vers le bas
/// - Intensité = vitesse (plus rapide = plus brillant)
fn velocity_color(vx: f32, vy: f32) -> u32 {
    let max_speed = INFLOW_SPEED * 1.2;
    let nx = (vx / max_speed).clamp(-1.0, 1.0);
    let ny = (vy / max_speed).clamp(-1.0, 1.0);

    // Composantes directionnelles
    let r = if nx > 0.0 { nx } else { 0.0 }; // Droite = rouge
    let b = if nx < 0.0 { -nx } else { 0.0 }; // Gauche = bleu
    let g = if ny > 0.0 { ny } else { 0.0 }; // Haut = vert
    let y = if ny < 0.0 { -ny } else { 0.0 }; // Bas = jaune (rouge+vert)

    // Mélanger les composantes
    let r_final = (r + y * 0.8).clamp(0.0, 1.0);
    let g_final = (g + y * 0.6).clamp(0.0, 1.0);
    let b_final = b.clamp(0.0, 1.0);

    ((r_final * 255.0) as u32) << 16 | ((g_final * 255.0) as u32) << 8 | ((b_final * 255.0) as u32)
}

#[inline]
fn lerp(a: f32, b: f32, t: f32) -> f32 {
    a + (b - a) * t
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
