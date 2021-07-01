extern crate criterion_plot;

use criterion_plot::prelude::*;
use itertools_num::linspace;
use std::path::Path;


pub fn plot_vector(y_values: Vec<f32>, dataname: &'static str, filename: &'static str, log: bool) {
    let x_values = linspace::<f32>(0.0, y_values.len() as f32, y_values.len()).collect::<Vec<_>>();

    // Make a new Figure to plot our vector:
    let mut f = Figure::new();

    // Configure settings for the output of the plot:
    f.set(Font("Helvetica"));
    f.set(FontSize(16.0));
    f.set(Output(Path::new(filename)));
    f.set(Size(1000, 400));

    // If log, set y axis to log mode:
    if log {
        f.configure(Axis::BottomX, |a| a
            .set(Scale::Logarithmic)
            //.set(Range::Limits(20.0, 44100.0/2.0))
        );
        f.configure(Axis::LeftY, |a| a
            .set(Scale::Logarithmic)
            //.set(Range::Limits(1000.0, 100000000.0))
        );
    }

    // Configure the key for the plot
    f.configure(Key, |k| {
        k.set(Boxed::Yes)
            .set(Position::Inside(Vertical::Top, Horizontal::Left))
    });

    // Plot the vector (in memory):
    f.plot(
        Lines {
            x: x_values,
            y: y_values,
        },
        |l| {
            l.set(Color::Rgb(255, 0, 0))
                .set(Label(dataname))
                .set(LineType::Solid)
        }
    );
    

    // Spit out the plot to a .svg file:
    f.draw()
        .ok()
        .and_then(|gnuplot| {
            gnuplot.wait_with_output()
                .ok()
                .and_then(|p| String::from_utf8(p.stderr).ok())
        }).expect("ERROR occurred while plotting");
}