Regression model:
<p align="center">y_k = b_0 + b_1 t_k + x_k,</p>
in particular, noise x_k could fit an AR(1) model:
<p align="center">x_k = a x_{k - 1} + sqrt(1 - a^2) e_k,</p>
where e_k is a white noise.
<p>
Results:
<p>
[1. Exact values of coverage factor K for normal e_k](exactK/K.pdf)
[2. Values of coverage factor K for e_k having two-side power (TSP) distribution (obtained by Monte Carlo simulation)](AR1Noise/K.pdf)
[3. Values of coverage factor K when x_k is 1/f^alpha noise (obtained by Monte Carlo simulation)](coloredNoise/K.pdf)
