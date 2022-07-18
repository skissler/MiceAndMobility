# Project notes

## 18 July 2022 

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

- A power-law decay with distance: 

$$ 
	\lambda(d) = \frac{k}{1 + d^\alpha}
$$ 

- An exponential decay with distance: 

$$ 
	\lambda(d) = k e^{-\phi d}
$$ 