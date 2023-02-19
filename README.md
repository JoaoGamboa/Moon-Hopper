# Trust-propelled Moon Hopper: Design and Trajectory optimization

## Motivation

The study focus on the sizing and trajectory assessment of a thurst-propelled lunar hopper. The system has an hybrid mobility system:(i) a wheeled mobility to traverse relatively smooth terrain and a (ii) propulsion system to perform assisted propelled hops. 
- The code in "Mass_Estimation.py" runs a small script to estimate the final mass of the hopper based on (i) the mean propellant consumption, (ii) the expeted coverable distance, (iii)the one-hop distance.
- The code in "Initial_Guess" provides a first analytical initial guess to the optimization routine.
- The code in "Optimization" provides an optimization routine based on the GEKKO python library for the hopper trajectory.

You are free re-use the code and change it. However, if you are using it for your project, please cite:
```
  @article{gamboa2022sizing,
    title={Sizing of a Propelled-Hopping System on the Moon},
    author={Gamboa, Joao and Rimani, Jasmine and and Lizy-Destrez, Stéphanie},
    publisher={Proceedings of the 73rd International Astronautical Congress, IAC},
    year={2022}
  }
```

## Licence 
This work is distributed under the MIT License.

## Status
The repository would be update with a new code on control optimization in the next few moths.
