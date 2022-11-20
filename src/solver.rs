use std::{f64::consts::PI, fs, io::Write};

pub struct Solver {
   delta: f64,
   eps: f64,
   l_x: f64,
   l_y: f64,
   r: f64,
   d: f64,
   v: f64,
   n_x: usize,
   n_y: usize,
   omega: f64,
   pixels: Vec<f64>,
   masks: Vec<bool>,
}

impl Solver {
   pub fn new(delta: f64, eps: f64, l_x: f64, l_y: f64, r: f64, d: f64, v: f64) -> Solver {
      let n_x: usize = (l_x / delta) as usize;
      let n_y: usize = (l_y / delta) as usize;

      let mu = 0.5 * ((PI / n_x as f64).cos() + (PI / n_y as f64).cos());
      let omega = 2. / (1. + (1. - mu.powf(2.)).sqrt());

      let mut solver = Solver {
         delta,
         eps,
         l_x,
         l_y,
         r,
         d,
         v,
         n_x,
         n_y,
         omega,
         pixels: vec![0f64; (n_x + 1) * (n_y + 1)],
         masks: vec![false; (n_x + 1) * (n_y + 1)],
      };
      solver.set_boundary();
      solver
   }

   pub fn create_double(solver: &Solver) -> Solver {
      let mut new_solver = Solver::new(
         solver.delta / 2.,
         solver.eps,
         solver.l_x,
         solver.l_y,
         solver.r,
         solver.d,
         solver.v,
      );
      for i in 0..=new_solver.n_x {
         let i2 = i / 2;
         for j in 0..=new_solver.n_y {
            if new_solver.get_mask(i, j) {
               continue;
            }
            let j2 = j / 2;
            if i % 2 == 1 {
               if j % 2 == 1 {
                  new_solver.set_pixel(
                     i,
                     j,
                     (solver.get_pixel(i2, j2)
                        + solver.get_pixel(i2 + 1, j2)
                        + solver.get_pixel(i2, j2 + 1)
                        + solver.get_pixel(i2 + 1, j2 + 1))
                        / 4.,
                  );
               } else {
                  new_solver.set_pixel(
                     i,
                     j,
                     (solver.get_pixel(i2, j2) + solver.get_pixel(i2 + 1, j2)) / 2.,
                  );
               }
            } else {
               if j % 2 == 1 {
                  new_solver.set_pixel(
                     i,
                     j,
                     (solver.get_pixel(i2, j2) + solver.get_pixel(i2, j2 + 1)) / 2.,
                  );
               } else {
                  new_solver.set_pixel(i, j, solver.get_pixel(i2, j2));
               }
            }
         }
      }
      drop(solver);
      new_solver
   }

   pub fn solve(&mut self) {
      loop {
         let error = self.step();
         println!("error = {:?}", error);
         if error < self.eps {
            break;
         }
      }
   }

   pub fn dump(&self, path: &str) {
      let mut file = fs::File::create(path).unwrap();
      for i in 0..=self.n_x {
         for j in 0..=self.n_y {
            let x = i as f64 * self.delta - self.l_x / 2.;
            let y = j as f64 * self.delta - self.l_y / 2.;
            file
               .write_all(
                  String::from(
                     x.to_string()
                        + " "
                        + &y.to_string()
                        + " "
                        + &self.get_pixel(i, j).to_string()
                        + "\n",
                  )
                  .as_bytes(),
               )
               .unwrap();
         }
         file.write_all(String::from("\n").as_bytes()).unwrap();
      }
   }

   fn set_boundary(&mut self) {
      let r_square = self.r.powf(2.);
      for i in 0..=self.n_x {
         let x = i as f64 * self.delta - self.l_x / 2.;
         let right_x_square = (x - self.d).powf(2.);
         let left_x_square = (x + self.d).powf(2.);
         for j in 0..=self.n_y {
            let y = j as f64 * self.delta - self.l_y / 2.;
            let y_square = y.powf(2.);

            if right_x_square + y_square <= r_square {
               self.set_mask(i, j, true);
               self.set_pixel(i, j, -(self.v / 2.));
            }
            if left_x_square + y_square <= r_square {
               self.set_mask(i, j, true);
               self.set_pixel(i, j, self.v / 2.);
            }
         }
      }
   }

   fn set_pixel(&mut self, i: usize, j: usize, v: f64) {
      self.pixels[i + (self.n_x + 1) * j] = v;
   }

   fn get_pixel(&self, i: usize, j: usize) -> f64 {
      self.pixels[i + (self.n_x + 1) * j]
   }

   fn set_mask(&mut self, i: usize, j: usize, v: bool) {
      self.masks[i + (self.n_x + 1) * j] = v;
   }

   fn get_mask(&self, i: usize, j: usize) -> bool {
      self.masks[i + (self.n_x + 1) * j]
   }

   fn calc_new_value(&self, top: f64, right: f64, bottom: f64, left: f64) -> f64 {
      0.25 * (top + right + bottom + left)
   }

   fn half_step(&mut self, is_even: bool) -> f64 {
      let mut error: f64 = 0.;
      for j in 0..=self.n_y {
         for i in (((is_even ^ (j % 2 != 0)) as usize)..=self.n_x).step_by(2) {
            if self.get_mask(i, j) {
               continue;
            }
            let value = self.get_pixel(i, j);
            let new_value = self.calc_new_value(
               self.get_pixel(i, self.add_y_index(j, 1)),
               self.get_pixel(self.add_x_index(i, 1), j),
               self.get_pixel(i, self.sub_y_index(j, 1)),
               self.get_pixel(self.sub_x_index(i, 1), j),
            );
            self.set_pixel(i, j, (1. - self.omega) * value + self.omega * new_value);
            error = (new_value - value).max(error);
         }
      }
      error
   }

   fn step(&mut self) -> f64 {
      self.half_step(true).max(self.half_step(false))
   }

   fn is_boundary(&self, i: usize, j: usize) -> bool {
      i == 0 || j == 0 || i == self.n_x || j == self.n_y
   }

   fn add_x_index(&self, i: usize, delta_i: usize) -> usize {
      i.saturating_add(delta_i).min(self.n_x)
   }

   fn add_y_index(&self, j: usize, delta_j: usize) -> usize {
      j.saturating_add(delta_j).min(self.n_y)
   }

   fn sub_x_index(&self, i: usize, delta_i: usize) -> usize {
      i.saturating_sub(delta_i).min(self.n_x)
   }

   fn sub_y_index(&self, j: usize, delta_j: usize) -> usize {
      j.saturating_sub(delta_j).min(self.n_y)
   }
}
