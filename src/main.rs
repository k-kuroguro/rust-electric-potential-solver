use solver::Solver;
use std::{
   io::Write,
   process::{Command, Stdio},
};

mod solver;

const DELTA: f64 = 0.1;
const EPS: f64 = 1e-3;
const L_X: f64 = 30.;
const L_Y: f64 = 25.;
const R: f64 = 0.4;
const D: f64 = 3.;
const V: f64 = 20.;

fn plot(data_file: &str, out_file: &str, l_x: f64, l_y: f64, v: f64) {
   let mut gnuplot = Command::new("gnuplot")
      .arg("-p")
      .stdin(Stdio::piped())
      .spawn()
      .unwrap();
   let stdin = gnuplot.stdin.as_mut().unwrap();
   writeln!(
      stdin,
      "
      set noxtics
      set noytics
      set nokey

      set contour
      set view map
      unset surface
      set cntrlabel onecolor
      set cntrparam levels increment {},1,{}

      set terminal pngcairo size {},{}
      set output '{}'

      splot [{}:{}] [{}:{}] '{}' w l lw 2 lc rgb '#000000'
      ",
      -(v / 2.) + 1.,
      v / 2. - 1.,
      l_x * 80.,
      l_y * 80.,
      out_file,
      -l_x / 2.,
      l_x / 2.,
      -l_y / 2.,
      l_y / 2.,
      data_file
   )
   .unwrap();
   gnuplot.wait().unwrap();
}

fn main() {
   let stop: usize = 1;
   let data_file = "data.dat";
   let out_file = "equipotential_lines.png";
   let mut solver = Solver::new(DELTA, EPS, L_X, L_Y, R, D, V);
   for i in 0..=stop {
      solver.solve();
      if i != stop {
         solver = Solver::create_double(&solver);
      }
   }
   solver.dump(data_file);
   plot(data_file, out_file, L_X, L_Y, V);
}
