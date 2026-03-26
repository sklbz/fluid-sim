//! Fluid simulation — live window or video export.
//!
//! Usage:
//!   cargo run --release                        # live window (ESC to quit)
//!   cargo run --release -- --video out.mp4     # render 300 frames to video
//!
//! Modélise la portance d'un profil NACA 4412 avec angle d'attaque.

use fluid_sim::FluidGrid;
use minifb::{Key, Window, WindowOptions};

// ── Simulation parameters ─────────────────────────────────────────────────────

/// Grille agrandie pour mieux capturer le sillage et la couche limite.
const GRID_W: u32 = 300;
const GRID_H: u32 = 150;
const CELL_SIZE: f32 = 1.0;

/// Pixels par cellule (réduit pour tenir dans une fenêtre raisonnable).
const SCALE: usize = 4;

/// Nombre d'itérations Gauss-Seidel par pas de temps.
const PRESSURE_ITERS: usize = 100;

/// Vitesse d'entrée (écoulement horizontal de gauche à droite).
const INFLOW_SPEED: f32 = 5.0;

/// Corde du profil en cellules.
const CHORD: u32 = 90;

/// Angle d'attaque en degrés (positif = bord d'attaque vers le haut).
const ANGLE_OF_ATTACK_DEG: f32 = 8.0;

const VIDEO_FRAMES: usize = 600;
const VIDEO_FPS: usize = 60;

const WIN_W: usize = GRID_W as usize * SCALE;
const WIN_H: usize = GRID_H as usize * SCALE;

struct ForceStats {
    lift_sum: f32,
    count: usize,
    last_log_time: std::time::Instant,
}

// ── NACA 4412 geometry ────────────────────────────────────────────────────────

/// Épaisseur demi-profil (upper/lower) normalisée du NACA 4412.
///
/// Paramètres NACA 4412 : m=0.04, p=0.4, t=0.12
/// Renvoie (y_upper, y_lower) en unités de corde pour x ∈ [0, 1].
fn naca4412(x_norm: f32) -> (f32, f32) {
    let m = 0.04_f32; // cambrure max
    let p = 0.40_f32; // position de la cambrure max
    let t = 0.12_f32; // épaisseur relative

    // Épaisseur symétrique NACA
    let yt = 5.0
        * t
        * (0.2969 * x_norm.sqrt() - 0.1260 * x_norm - 0.3516 * x_norm.powi(2)
            + 0.2843 * x_norm.powi(3)
            - 0.1015 * x_norm.powi(4));

    // Ligne de cambrure
    let yc = if x_norm < p {
        m / (p * p) * (2.0 * p * x_norm - x_norm * x_norm)
    } else {
        m / ((1.0 - p) * (1.0 - p)) * ((1.0 - 2.0 * p) + 2.0 * p * x_norm - x_norm * x_norm)
    };

    // Pente de la cambrure (pour orienter les normales)
    let dyc_dx = if x_norm < p {
        2.0 * m / (p * p) * (p - x_norm)
    } else {
        2.0 * m / ((1.0 - p) * (1.0 - p)) * (p - x_norm)
    };

    let theta = dyc_dx.atan();

    let xu = x_norm - yt * theta.sin();
    let yu = yc + yt * theta.cos();
    let xl = x_norm + yt * theta.sin();
    let yl = yc - yt * theta.cos();

    // On n'a besoin que des ordonnées (on interpole sur une grille cartésienne)
    let _ = (xu, xl);
    (yu, yl)
}

// ── Airfoil mask generation ───────────────────────────────────────────────────

/// Remplit le masque `solid` avec le profil NACA 4412 tourné de `aoa_deg` degrés.
///
/// Le bord d'attaque est positionné à 1/4 de la grille en x,
/// et centré verticalement.
fn build_airfoil(grid: &mut FluidGrid, aoa_deg: f32) {
    let (w, h) = grid.cell_count;
    let chord = CHORD as f32;

    // Centre du profil dans la grille
    let cx = w as f32 * 0.28; // bord d'attaque à ~28 % → bord de fuite à ~58 %
    let cy = h as f32 * 0.50;

    let aoa_rad = -aoa_deg.to_radians(); // signe : rotation CCW = bord d'attaque vers le haut
    let (sin_a, cos_a) = (aoa_rad.sin(), aoa_rad.cos());

    for x in 0..w {
        for y in 0..h {
            // Position de la cellule dans le repère du profil (non tourné)
            let dx = x as f32 + 0.5 - cx;
            let dy = y as f32 + 0.5 - cy;

            // Rotation inverse pour se ramener dans le repère du profil
            let local_x = dx * cos_a + dy * sin_a;
            let local_y = -dx * sin_a + dy * cos_a;

            // Coordonnée normalisée sur la corde
            let x_norm = local_x / chord;

            if x_norm < 0.0 || x_norm > 1.0 {
                continue;
            }

            let (yu, yl) = naca4412(x_norm);

            // Convertir en cellules
            let y_upper = yu * chord;
            let y_lower = yl * chord;

            if local_y >= y_lower && local_y <= y_upper {
                grid.solid[(x, y)] = true;
            }
        }
    }
}

// ── Entry point ───────────────────────────────────────────────────────────────

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let video_out = args
        .windows(2)
        .find(|w| w[0] == "--video")
        .map(|w| w[1].clone());

    let mut grid = FluidGrid::new((GRID_W, GRID_H), CELL_SIZE);
    build_airfoil(&mut grid, ANGLE_OF_ATTACK_DEG);

    if let Some(path) = video_out {
        run_video(&mut grid, &path);
    } else {
        run_window(&mut grid);
    }
}

// ── Inflow injection ──────────────────────────────────────────────────────────

/// Injecte un écoulement uniforme sur toute la face gauche.
/// On ajoute une légère composante verticale selon l'angle d'attaque
/// (revient à simuler un vent incident incliné).
fn inject_inflow(grid: &mut FluidGrid) {
    let h = grid.cell_count.1;
    let aoa = ANGLE_OF_ATTACK_DEG.to_radians();
    let vx = INFLOW_SPEED * aoa.cos();
    let vy = INFLOW_SPEED * aoa.sin(); // positif = vers le haut

    for y in 0..h {
        grid.velocities_x[(0, y)] = vx;
    }
    // La composante verticale s'applique sur la face y=0→h (face gauche)
    // via les faces v de la colonne x=0.
    for y in 0..=h {
        grid.velocities_y[(0, y)] = vy;
    }
}

// ── Rendering ─────────────────────────────────────────────────────────────────

/// Colorisation :
///   - Gris foncé  = cellule solide (profil)
///   - Rouge       = vitesse vers la droite (surpression extrados)
///   - Bleu        = vitesse vers la gauche
///   - Vert        = vitesse vers le haut (portance)
///   - Jaune       = vitesse vers le bas (sillage)
fn render(grid: &FluidGrid, buffer: &mut Vec<u32>) {
    let (w, h) = (grid.cell_count.0 as usize, grid.cell_count.1 as usize);

    for x in 0..w {
        for y in 0..h {
            let color = if grid.solid[(x as u32, y as u32)] {
                0x2C2C2A // gris anthracite = solide
            } else {
                let (vx, vy) = grid.velocity_at_center(x as u32, y as u32);
                velocity_color(vx, vy)
            };

            // y=0 est en bas de la physique, mais en haut du buffer écran.
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

fn velocity_color(vx: f32, vy: f32) -> u32 {
    let max_speed = INFLOW_SPEED * 1.5;
    let nx = (vx / max_speed).clamp(-1.0, 1.0);
    let ny = (vy / max_speed).clamp(-1.0, 1.0);

    let r = if nx > 0.0 { nx } else { 0.0 };
    let b = if nx < 0.0 { -nx } else { 0.0 };
    let g = if ny > 0.0 { ny } else { 0.0 };
    let y_comp = if ny < 0.0 { -ny } else { 0.0 };

    let r_final = (r + y_comp * 0.8).clamp(0.0, 1.0);
    let g_final = (g + y_comp * 0.6).clamp(0.0, 1.0);
    let b_final = b.clamp(0.0, 1.0);

    ((r_final * 255.0) as u32) << 16 | ((g_final * 255.0) as u32) << 8 | (b_final * 255.0) as u32
}

// ── Live window mode ──────────────────────────────────────────────────────────

fn run_window(grid: &mut FluidGrid) {
    let mut window = Window::new(
        "Fluid sim — NACA 4412  (ESC = quit)",
        WIN_W,
        WIN_H,
        WindowOptions::default(),
    )
    .expect("Failed to create window");

    let mut buffer = vec![0u32; WIN_W * WIN_H];
    let mut frame = 0usize;

    while window.is_open() && !window.is_key_down(Key::Escape) {
        inject_inflow(grid);
        grid.step(PRESSURE_ITERS);
        render(grid, &mut buffer);
        window.update_with_buffer(&buffer, WIN_W, WIN_H).unwrap();

        frame += 1;
        if frame % 60 == 0 {
            println!("frame {frame}  div={:.4}", grid.max_divergence());
        }
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

        let bytes: Vec<u8> = buffer
            .iter()
            .flat_map(|&p| {
                [
                    ((p >> 16) & 0xFF) as u8,
                    ((p >> 8) & 0xFF) as u8,
                    (p & 0xFF) as u8,
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

    drop(ffmpeg.stdin.take());
    ffmpeg.wait().expect("ffmpeg did not finish successfully");
    println!("Done → {output_path}");
}
