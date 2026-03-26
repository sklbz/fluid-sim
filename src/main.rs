use fluid_sim::FluidGrid;
use std::fs::write;
use std::path::Path;
use std::{thread::sleep, time::Duration};

fn main() {
    let mut grid = FluidGrid::new((10, 10), 1.0);
    grid.velocities_x.randomize(10.0..=10.0);
    grid.velocities_y.randomize(10.0..=10.0);

    // display_divergence(&grid);
    // display_pressure(&grid);

    let path = Path::new("/home/sklbz/code/fluid-sim/frames/-1");
    match write(path, dump_json(&grid)) {
        Ok(_) => println!("Successfully wrote to file"),
        Err(e) => println!("Failed to write to file: {}", e),
    };

    grid.advect_velocities();

    let path = Path::new("/home/sklbz/code/fluid-sim/frames/-1advection");
    match write(path, dump_json(&grid)) {
        Ok(_) => println!("Successfully wrote to file"),
        Err(e) => println!("Failed to write to file: {}", e),
    };

    grid.gauss_seidel();
    let path = Path::new("/home/sklbz/code/fluid-sim/frames/-1gaussseidel");
    match write(path, dump_json(&grid)) {
        Ok(_) => println!("Successfully wrote to file"),
        Err(e) => println!("Failed to write to file: {}", e),
    };

    return;
    for i in 0..100 {
        grid.advect_velocities();
        grid.gauss_seidel();

        let path = Path::new("/home/sklbz/code/fluid-sim/frames/").join(i.to_string());
        match write(path, dump_json(&grid)) {
            Ok(_) => println!("Successfully wrote to file"),
            Err(e) => println!("Failed to write to file: {}", e),
        };
    }

    return;
    loop {
        println!("\x1b[2J\x1b[H");

        // grid.advect_velocities();
        grid.gauss_seidel();

        display_divergence(&grid);
        display_pressure(&grid);
        println!("{}", dump_json(&grid));

        sleep(Duration::from_millis(5000));
    }

    for _ in 0..1 {
        for _ in 0..40 {
            println!("\x1b[2J\x1b[H");

            display_divergence(&grid);
            display_pressure(&grid);
            // display_velocities(&grid);

            grid.pressure_solve();
            sleep(Duration::from_millis(10));
        }
        grid.update_velocities();
        sleep(Duration::from_millis(20));
    }
}

fn dump_json(grid: &FluidGrid) -> String {
    let cx = grid.cell_count.0 as usize;
    let cy = grid.cell_count.1 as usize;

    let fmt = |v: &[f32]| {
        let s: Vec<String> = v.iter().map(|x| format!("{:.4}", x)).collect();
        format!("[{}]", s.join(","))
    };

    let pressure: Vec<String> = (0..cx)
        .map(|x| {
            let col: Vec<f32> = (0..cy)
                .map(|y| grid.pressure_map[(x as u32, y as u32)])
                .collect();
            fmt(&col)
        })
        .collect();

    let vx: Vec<String> = (0..=cx)
        .map(|x| {
            let col: Vec<f32> = (0..cy)
                .map(|y| grid.velocities_x[(x as u32, y as u32)])
                .collect();
            fmt(&col)
        })
        .collect();

    let vy: Vec<String> = (0..cx)
        .map(|x| {
            let col: Vec<f32> = (0..=cy)
                .map(|y| grid.velocities_y[(x as u32, y as u32)])
                .collect();
            fmt(&col)
        })
        .collect();

    format!(
        r#"{{"cell_count":[{},{}],"cell_size":{:.4},"pressure":[{}],"velocities_x":[{}],"velocities_y":[{}]}}"#,
        cx,
        cy,
        grid.cell_size,
        pressure.join(","),
        vx.join(","),
        vy.join(",")
    )
}
fn display_divergence(grid: &FluidGrid) {
    println!("Divergence:");
    for x in 0..grid.cell_count.0 {
        for y in 0..grid.cell_count.1 {
            let div = grid.divergence_at_cell(x, y);
            if div.abs() > 0.5 {
                if div > 0.0 {
                    print!("\x1b[31m");
                } else {
                    print!("\x1b[34m");
                }
            } else {
                print!("\x1b[32m");
            }
            print!("{:^5.1}", if div.abs() < 1e-2 { 0.0 } else { div });
            print!("\x1b[0m");
        }
        println!();
    }
}

fn display_pressure(grid: &FluidGrid) {
    println!("Pressure:");
    for x in 0..grid.cell_count.0 {
        for y in 0..grid.cell_count.1 {
            print!("{:^5.1}", grid.get_pressure(x, y));
        }
        println!();
    }
}

/// Display velocities
fn display_velocities(grid: &FluidGrid) {
    println!("Velocities X:");
    for x in 0..grid.cell_count.0 {
        for y in 0..grid.cell_count.1 {
            print!("{:^5.1}", grid.velocities_x[(x, y)]);
        }
        println!();
    }
    println!("Velocities Y:");
    for x in 0..grid.cell_count.0 {
        for y in 0..grid.cell_count.1 {
            print!("{:^5.1}", grid.velocities_y[(x, y)]);
        }
        println!();
    }
}
