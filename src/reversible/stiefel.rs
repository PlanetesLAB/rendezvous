use super::stumpff;

pub fn gs(_gs: &mut [f64], _beta: f64, _x: f64) {}

pub fn gs3(gs: &mut [f64], beta: f64, x: f64) {
    let x2 = x * x;
    stumpff::cs3(gs, beta * x2);
    gs[1] *= x;
    gs[2] *= x2;
    gs[3] *= x2 * x;
}
