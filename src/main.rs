use fluid_sim::FluidGrid;
use std::{thread::sleep, time::Duration};

fn main() {
    let mut grid = FluidGrid::new((10, 10), 1.0);
    grid.velocities_x.randomize(-10.0..=10.0);
    grid.velocities_y.randomize(-10.0..=10.0);

    display_divergence(&grid);
    display_pressure(&grid);

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
    display_divergence(&grid);
    display_pressure(&grid);
    // display_velocities(&grid);
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
