<div id="top"></div>

# Bending of Bernoulli Beams - Numerical Analysis (**B3NA**)

Python implementation of the finite element method (FEM) to solve the bending equation for an elastic Bernoulli beam. For the sake of simplicity, we do not use Sobolev spaces and work with piecewise differentiable functions instead. The repository also includes tests and experiments to validate the correct functionality of the package and to illustrate how to get started with it. 

> This project is part of the course *Project in Numerical Anlysis* of the MSc in Computer Simulations for Science and Engineering, lectured by M. Karow at TU Berlin (SoSe 2022)

## Built With 

* [NumPy](https://numpy.org/)
* [SciPy](https://scipy.org/)
* [Matplotlib](https://matplotlib.org/)
* [Seaborn](https://seaborn.pydata.org/)

<p align="right">(<a href="#top">back to top</a>)</p>


## Authors (alphabetically)

F. de Andres - deandresvertferran@gmail.com

F. Sindy - f.sindy@student.tudelft.nl

V. Yelnyk - volodymyr.yelnyk@gmail.com

<p align="right">(<a href="#top">back to top</a>)</p>


## Case of application and usage (problem statement) 

### Pure bending of a single beam, both statically and dynamically: 

![imagen](https://user-images.githubusercontent.com/92535468/166642289-8f03800e-aa82-49ea-a4b5-e03723cb9f97.png)

### Vibration of frameworks of beams: 

*TO BE INCLUDED*

<p align="right">(<a href="#top">back to top</a>)</p>


## Usage and install

0. Read the documentation and be aware of the [license terms](https://github.com/ferdean/beam-num-analysis/blob/main/LICENSE)
1. Git-clone this repo and change directory 
    
    <pre>git clone https://github.com/ferdean/beam-num-analysis</pre>
    <pre>cd beam-num-analysis</pre>
    
2. Install modules using pip.
   <pre>pip install -r requirements.txt</pre>

<p align="right">(<a href="#top">back to top</a>)</p>

## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>


## Roadmap 

- [ ] Single beam, static case
    - [x] Implement static deformation of cantilever beam
    - [ ] Implement different BC (e. g. beam supported at both ends)
    - [x] Develop analytic solutions to design tests and validation cases
    - [x] Clean code, collect everything in packages/libraries and objects
- [ ] Single beam, dynamic case
    - [x] Implement time-dependent solver (and repeat all points of static case)
    - [x] Test homogeneous case with a deformed initial condition
    - [x] Test non-homogeneous cases
    - [x] Develop visualization library + GUI supporting animations (Newmark-method)
    - [ ] Make forces and BCs interactive in the GUI
- [ ] Eigenvalue analysis
- [ ] Vibration of frameworks of beams
    - [ ] Implementation of a 2D solver for frameworks
    - [ ] Implementation of a 3D solver for framekorks
- [ ] Documentation and cleaning
    - [ ] Oral presentation (10-20 min video)
    - [ ] Written report (around 20 pages)
    - [ ] Update README.md and repository documentation. 
 
 <p align="right">(<a href="#top">back to top</a>)</p>

