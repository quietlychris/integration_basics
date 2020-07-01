#![feature(proc_macro_hygiene)]
use inline_python::{python, Context};

fn main() {
    // Prints the results of an accleration under Euler's method. Experiment with the time-step
    // explicit_euler(0.01);

    // Simulates the movement of spring-mass-damper using either explicit or semi-implicit Euler's method
    // The explicit method will improperly diverge. Changes in the physical parameters must be altered below
    // Specify "explicit" or "semi-implicit" as the integration mode
    // spring_mass_damper("explicit");

    // Simulates the same SMD system using a 4th-order Runge-Kutta method
    let (time, position) = rk4();

    // The basics of returning the time/position vectors from each method has been implemented for RK4
    // If you're interested in overlaying the plots, it should be fairly trivial to alter the spring_mass_damper()
    // function to return them, then build a graph using the basics of the plot() function below. 
}

fn explicit_euler(dt: f64) {
    let mut t = 0.;

    let mut velocity = 0.;
    let mut position = 0.;
    let mut force = 10.;
    let mut mass = 1.;

    while t <= 10. {
        position = position + velocity * dt;
        velocity = velocity + (force / mass) * dt;
        println!(
            "t={:.2}:    position = {}    velocity = {}",
            t, position, velocity
        );
        t += dt;
    }
}

fn spring_mass_damper(mode: &str) {
    let mut position_vec = Vec::new();
    let mut time_vec = Vec::new();

    let mut t = 0.;
    let dt = 0.01;
    let mut velocity = 0.;
    let mut position = 1000.;
    let mut mass = 1.;
    let k = 15.; // Spring constant
    let b = 0.1; // Damping coefficient

    while t <= 100. {
        let force = (-k * position) - (b * velocity);
        match mode {
            "explicit" => {
                position = position + velocity * dt;
                velocity = velocity + (force / mass) * dt;
            }
            "semi-implicit" => {
                velocity = velocity + (force / mass) * dt;
                position = position + velocity * dt;
            }
            _ => panic!("Please specify 'explicit' or 'implicit' as the mode"),
        };
        println!(
            "t={:.2}:    position = {}    velocity = {}",
            t, position, velocity
        );

        position_vec.push(position); // Adds position and time to vectors for graphing
        time_vec.push(t);
        t += dt;
    }

    let ctx = Context::new();
    ctx.run(python! {

        import numpy as np;
        import matplotlib.pyplot as plt;

        plt.plot('time_vec,'position_vec,"b-");

        plt.show();
    });
}

/// Implementing RK4

#[derive(Debug, Clone, Copy)]
struct State {
    x: f64, // position
    v: f64, // velocity
}

// It's nice to have convenience methods, so we'll make one that can build a simple
// struct with a default values set to zero
impl State {
    fn default() -> Self {
        State { x: 0., v: 0. }
    }
}

struct Derivative {
    dx: f64, // dx/dt = velocity
    dv: f64, // dv/dt = acceleration
}

impl Derivative {
    fn default() -> Self {
        Derivative { dx: 0., dv: 0. }
    }
}

fn evaluate(initial: &State, t: f64, dt: f64, d: &Derivative) -> Derivative {
    let mut state = State::default(); // We instantiate a struct

    state.x = initial.x + d.dx * dt;
    state.v = initial.v + d.dv * dt;

    let mut output = Derivative::default();
    output.dx = state.v;
    output.dv = acceleration(&state, t + dt);
    output
}

fn acceleration(state: &State, t: f64) -> f64 {
    let k = 15.;
    let b = 0.1;
    (-k * state.x) - (b * state.v)
}

fn integrate(state: &mut State, t: f64, dt: f64) {
    let mut a = Derivative::default();
    let mut b = Derivative::default();
    let mut c = Derivative::default();
    let mut d = Derivative::default();

    a = evaluate(state, t, 0., &Derivative::default());
    b = evaluate(state, t, dt * 0.5, &a);
    c = evaluate(state, t, dt * 0.5, &b);
    d = evaluate(state, t, dt, &c);

    let dxdt = (1. / 6.) * (a.dx + 2. * (&b.dx + &c.dx) + d.dx);
    let dvdt = (1. / 6.) * (a.dv + 2. * (&b.dv + &c.dv) + d.dv);

    state.x = state.x + dxdt * dt;
    state.v = state.v + dvdt * dt;
}

fn rk4() -> (Vec<f64>,Vec<f64>) {
    let mut initial = State { x: 1000., v: 0. };
    let mut position_vec = Vec::new();
    let mut time_vec = Vec::new();

    let mut t = 0.;
    let dt = 0.01;

    while t < 100. {
        integrate(&mut initial, t, dt);
        assert!(initial.x.is_normal() == true);

        position_vec.push(initial.x); // Adds position and time to vectors for graphing
        time_vec.push(t);
        t += dt;
    }

    plot(&time_vec,&position_vec,"r-");

    (time_vec, position_vec)

}

fn plot(x: &Vec<f64>, y: &Vec<f64>,line_type: &str) {
    let ctx = Context::new();
    ctx.run(python! {

        import numpy as np;
        import matplotlib.pyplot as plt;

        plt.plot('x,'y,'line_type);

        plt.show();
    });
}