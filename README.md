# global-optimization
Unknown function approximation, given input-output measurements, using Genetic Algorithm 

Preliminary implementation of a [Genetic Algorithm](https://www.geeksforgeeks.org/genetic-algorithms/), utilized to approximate an unknown function, given 
input-output measurements. Flexibility is provided to obtain different results and bias the approximation towards a desired direction, through selection of
multiple parameters such as, the size of the population, the k more accurate chromosomes to pick, the number of generations, and the tester-set of the input 
values.

Version-1
```
1. Less values given as input
2. Computational requirements are kept on a reasonable level
3. Both the number of generations and the size of population can reach high values, without punishing effect on performance
4. Less accurate approximation of the function most of the times
5. Result is biased towards the area of the specific values, given as input
```

Version-2
```
1. Whole grid given as input, containing every possible value inside the range
2. Computational requirements increase drastically
3. Convention must be made about the selection of parameters, as increase of all at the same time leads to very poor performance
4. More accurate approximation most of the time
5. Result is less biased towards the area of the input
```
