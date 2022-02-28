# global-optimization
Unknown function approximation, given input-output measurements, using Genetic Algorithm 

Preliminary implementation of a [Genetic Algorithm](https://www.geeksforgeeks.org/genetic-algorithms/), utilized to approximate an unknown function, given 
input-output measurements. Flexibility is provided to obtain different results and bias the approximation towards a desired direction, through selection of
multiple parameters such as, the size of the population, the k more accurate chromosomes to pick, the number of generations, and the tester-set of the input 
values.

**note**

For optimal approximation, the threshold (desired MSE) must be kept low enough, anywhere in the range 0.0001-0.0005. This however, translates to excessive computational
requirements, and a pretty poor performance.

