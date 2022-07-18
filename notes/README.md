# Project notes

<!-- <details><summary> -->
<h2> 18 July 2022 </h2>
<!-- </summary> -->

The goal of this code is to __infer how the probability of transmission depends on the proximity and duration of contact.__ 

We're starting with a straightforward toy model that specifies the instantaneous hazard (force) of infection, $\lambda$. This hazard can take a few forms: 

- A step function: 

$$
	\lambda(d) = 
	\begin{cases} 
		k & d \leq d* \\ 
		0 & d > d*
	\end{cases} 
$$ 

<p align="center">
<img src="2022-07-18-stepfun.png" style="width:50%">
</p>

- A power-law decay with distance: 

$$ 
	\lambda(d) = \frac{k}{1 + d^\alpha}
$$ 

- An exponential decay with distance: 

$$ 
	\lambda(d) = k e^{-\phi d}
$$ 

The questions become: 

- Given some observations (locations over time, timing of infection), how precisely can we determine the kernel parameters ($k, d* , \alpha, \phi$)? 
- Under which circumstances can we distinguish between these models? 
- When does it matter to be able to distinguish between models? 
- What sorts of experiments do we need to run (sample size, frequency of observation, precision with which we need to know epidemiologic links) to measure the infection kernel with sufficient accuracy to inform interventions (_e.g.,_ isolation period, gathering size restrictions, general risk communication)? 

<!-- </details> -->

</body>