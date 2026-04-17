**Overview:** This project develops a multi-period portfolio optimization model for a charitable foundation aiming to grow its 
assets from $40M to $100M over 10 years while continuing annual spending and receiving donations. The model determines optimal
investment decisions over time under uncertainty, balancing growth objectives with risk and donor-driven constraints.

**Methodology:**
1. Problem Setup
     - Initial wealth: $40M
     - Investment horizon: 10 years
     - Annual spending and donations incorporated into cash flows
     - Investment universe: 7 pre-defined fund mixes across asset classes
2. Return Modeling
    - Modeled asset returns as jointly lognormal distributions
    - Estimated parameters:
    - Expected returns (μ)
    - Volatility (σ)
    - Correlation structure
    - Generated scenarios using Monte Carlo simulation (5,000 paths)
3. Optimization Approach - Formulated as a dynamic programming problem using Bellman recursion
    - At each time step:
         1. Chose the optimal fund allocation
         2. Maximized expected future utility of wealth
    - Constraints:
        1. Wealth evolves with spending and donations
        2. Investment decisions limited to predefined fund mixes
4. Utility Function - Modeled nonlinear, piecewise utility to reflect client preferences:
       a. Strong preference for reaching $100M target
       b. Increasing reward for exceeding $110M
       c. Penalties for underperformance
6. Simulation & Evaluation - Forward-simulated optimal policy over 5,000 paths
      - Key outputs:
          - Expected terminal wealth
          - Probability of reaching $100M and $110M
          - Optimal fund allocation over time
       
**Conclusion:** In conclusion, this project demonstrates how dynamic programming can be used to determine optimal investment 
strategies under uncertainty for long-term financial goals. By incorporating stochastic returns, spending, donations, and 
a nonlinear utility function, the model adapts portfolio decisions over time to maximize expected outcomes while managing 
risk. The results highlight the importance of flexible, state-dependent strategies in improving the likelihood of achieving 
target wealth levels. 
